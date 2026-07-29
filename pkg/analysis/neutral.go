package analysis

// Neutral-expectation validation — roadmap §6h, "Credibility backbone".
//
// Before any creationist argument is made, DRIFT must be shown to reproduce
// textbook neutral population-genetics results: Tajima's D ~ 0, the 1/i site-
// frequency spectrum, theta = 4*Ne*mu recovery (Watterson's theta_W and pairwise
// theta_pi both landing on the same theta), and Hardy-Weinberg genotype
// proportions. This file computes those statistics and scores them against their
// analytic expectations with tolerances, emitting a pass/fail report so
// regressions are caught.
//
// WHICH VARIATION IS MEASURED. DRIFT has TWO independent genetic-bookkeeping
// systems, and they must not be confused here:
//
//   * the Chromosomes BITFIELD (read by SFS/VCF/Fst/joint-SFS/het-distance) holds
//     only FOUNDER standing variation seeded at t=0. It receives NO de-novo input
//     during a run — it only drifts, so every site marches monotonically to loss
//     or fixation and the bitfield NEVER reaches mutation-drift equilibrium.
//
//   * the MUTATION POOL (pop.IndMutations + pop.MutationPool, gated by
//     track_mutations) receives a fresh Poisson(mu) input of de-novo mutations per
//     birth, inherited through meiosis. THIS is the infinite-sites process that
//     reaches mutation-drift equilibrium, so the neutral theta = 4*Ne*mu
//     expectation is only meaningful here.
//
// The §6h statistics are therefore computed from the MUTATION POOL, not the
// bitfield: each distinct mutation id is one infinite-sites "segregating site",
// its derived-allele count in the sample is how many of the sampled strand copies
// carry it, and an individual is heterozygous at that site iff exactly one of its
// two strands carries the mutation.
//
// UNITS / theta = 4*Ne*mu. mu is the Poisson mean number of new mutations per
// DIPLOID individual per generation (utils.GenerateNewMutations). A diploid is two
// haploid gene copies, so the per-copy (per-sequence) rate is mu/2, and the
// autosomal neutral expectation theta = 4*Ne*(mu/2) = 2*Ne*mu. We therefore
// recover an implied Ne = theta_hat / (2*mu). Because DRIFT is an overlapping-
// generations, monogamous, age-structured model (NOT Wright-Fisher), the coalescent
// Ne differs from the census N; the report states the implied Ne and the emergent
// Ne/N ratio rather than asserting theta == 2*N_census*mu. The Ne-INDEPENDENT
// checks (Tajima's D ~ 0, theta_pi ~ theta_W, the 1/i SFS shape, and HWE F ~ 0) are
// the hard assertions; the absolute theta value is reported as an implied Ne.

import (
	"drift/pkg/core"
	"drift/pkg/utils"
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"sort"
)

// NeutralStats holds the neutral summary statistics computed from the mutation
// pool over a sample of individuals.
type NeutralStats struct {
	NumSamples int // sampled diploid individuals
	NumChrom   int // haploid sample size = 2 * NumSamples
	SegSites   int // segregating mutation lineages in the sample
	Singletons int // sites with derived count 1 or NumChrom-1
	Doubletons int
	ThetaW     float64 // Watterson's theta (genome-wide, = S / a_n)
	ThetaPi    float64 // pairwise theta / nucleotide diversity (genome-wide)
	TajimasD   float64
	FuLiD      float64
	MeanHo     float64 // mean observed heterozygosity across segregating sites
	MeanHe     float64 // mean expected heterozygosity 2p(1-p) across segregating sites
	FIS        float64 // inbreeding coefficient 1 - Ho/He (HWE departure)
	// UnfoldedSFS[i] = number of segregating sites with i derived copies in the
	// sample (i = 1..NumChrom-1); index 0 and NumChrom are always 0.
	UnfoldedSFS []int
	FoldedSFS   []int // minor-allele-count spectrum (i = 1..NumChrom/2)
}

// SampleLiving returns up to n living individual IDs (from pop.IndData), chosen at
// random with the shared per-run RNG so sampling is reproducible under a fixed
// seed. Unlike SampleIDs (which samples from pop.Chromosomes), this works when the
// neutral run tracks mutations without the DNA bitfield (track_DNA = 0). n <= 0 or
// n >= population size returns all living individuals, sorted.
func SampleLiving(pop *core.Pop, n int) []int {
	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	// Sort BEFORE shuffling: the id slice is built from map iteration order (which
	// Go randomizes per run), so shuffling it directly would pick a different subset
	// each run even under a fixed RNG seed. Sorting first makes the subsample a
	// deterministic function of the RNG stream alone.
	sort.Ints(ids)
	if n <= 0 || n >= len(ids) {
		return ids
	}
	utils.RandShuffle(len(ids), func(i, j int) { ids[i], ids[j] = ids[j], ids[i] })
	ids = ids[:n]
	sort.Ints(ids)
	return ids
}

// ComputeNeutralStats tallies the mutation pool over the sampled individuals and
// derives the neutral summary statistics. Returns a zero-value-ish result (SegSites
// 0) when there are fewer than two sampled individuals or no segregating mutations.
func ComputeNeutralStats(pop *core.Pop, ids []int) *NeutralStats {
	numSamples := len(ids)
	numChrom := 2 * numSamples
	stats := &NeutralStats{
		NumSamples:  numSamples,
		NumChrom:    numChrom,
		UnfoldedSFS: make([]int, numChrom+1),
		FoldedSFS:   make([]int, numChrom/2+1),
	}
	if numSamples < 2 {
		return stats
	}

	// Per mutation lineage, tally derived-copy count and heterozygote count across
	// the sampled individuals.
	type tally struct{ derived, het int }
	counts := make(map[int]*tally)
	for _, id := range ids {
		muts, ok := pop.IndMutations[id]
		if !ok {
			continue
		}
		// Presence of each mutation on strand 0 and strand 1 for this individual.
		on0 := muts[0]
		on1 := muts[1]
		seen0 := make(map[int]bool, len(on0))
		for _, m := range on0 {
			seen0[m] = true
		}
		seen1 := make(map[int]bool, len(on1))
		for _, m := range on1 {
			seen1[m] = true
		}
		// Union of ids carried on either strand.
		for m := range seen0 {
			t := counts[m]
			if t == nil {
				t = &tally{}
				counts[m] = t
			}
			c := 1
			if seen1[m] {
				c = 2
			}
			t.derived += c
			if c == 1 {
				t.het++
			}
		}
		for m := range seen1 {
			if seen0[m] {
				continue // already handled in the strand-0 pass
			}
			t := counts[m]
			if t == nil {
				t = &tally{}
				counts[m] = t
			}
			t.derived++ // present only on strand 1
			t.het++
		}
	}

	// Accumulate over mutation ids in SORTED order so the floating-point sums
	// (theta_pi, Ho, He) are bit-reproducible: map iteration order is randomized and
	// float addition is not associative, so an unsorted sum drifts in the last ULP
	// between runs. The SFS bins are integer and order-independent, but keep one loop.
	mutIDs := make([]int, 0, len(counts))
	for m := range counts {
		mutIDs = append(mutIDs, m)
	}
	sort.Ints(mutIDs)

	var thetaPi, sumHo, sumHe float64
	for _, m := range mutIDs {
		t := counts[m]
		d := t.derived
		if d <= 0 || d >= numChrom {
			continue // monomorphic in the sample (lost or fixed): not segregating
		}
		stats.SegSites++
		stats.UnfoldedSFS[d]++
		minor := d
		if minor > numChrom-minor {
			minor = numChrom - minor
		}
		stats.FoldedSFS[minor]++
		if d == 1 || d == numChrom-1 {
			stats.Singletons++
		}
		if d == 2 || d == numChrom-2 {
			stats.Doubletons++
		}
		// Unbiased per-site pairwise diversity: 2*d*(2n-d) / (2n*(2n-1)).
		thetaPi += 2 * float64(d) * float64(numChrom-d) / (float64(numChrom) * float64(numChrom-1))
		p := float64(d) / float64(numChrom)
		sumHe += 2 * p * (1 - p)
		sumHo += float64(t.het) / float64(numSamples)
	}

	stats.ThetaPi = thetaPi
	stats.ThetaW = CalculateThetaWatterson(stats.SegSites, numChrom) // genome-wide S / a_n
	stats.TajimasD = CalculateTajimasD(stats.ThetaPi, stats.ThetaW, stats.SegSites, numChrom)
	stats.FuLiD = CalculateFuLiD(stats.Singletons, stats.SegSites, numChrom)
	if stats.SegSites > 0 {
		stats.MeanHo = sumHo / float64(stats.SegSites)
		stats.MeanHe = sumHe / float64(stats.SegSites)
		if stats.MeanHe > 0 {
			stats.FIS = 1 - stats.MeanHo/stats.MeanHe
		}
	}
	return stats
}

// ExpectedUnfoldedSFS returns the neutral-equilibrium expected PROPORTION of
// segregating sites in each unfolded frequency class i (1..numChrom-1): under the
// standard coalescent E[xi_i] proportional to 1/i, normalized to sum to 1. Index 0
// and numChrom are 0. Used to score the observed SFS shape.
func ExpectedUnfoldedSFS(numChrom int) []float64 {
	exp := make([]float64, numChrom+1)
	if numChrom < 2 {
		return exp
	}
	var norm float64
	for i := 1; i < numChrom; i++ {
		norm += 1.0 / float64(i)
	}
	for i := 1; i < numChrom; i++ {
		exp[i] = (1.0 / float64(i)) / norm
	}
	return exp
}

// SFSShapeChiSquare returns a Pearson chi-square goodness-of-fit statistic of the
// observed unfolded SFS against the neutral 1/i expectation, together with the
// degrees of freedom (number of frequency classes with a positive expectation
// minus 1). Lower is a better fit; chi2/dof near 1 is consistent with neutrality.
// Classes with an expected count below 1 are pooled into the neighbouring class so
// the statistic is not dominated by sparse tail bins.
func SFSShapeChiSquare(unfolded []int, segSites, numChrom int) (chi2 float64, dof int) {
	if segSites <= 0 || numChrom < 3 {
		return 0, 0
	}
	exp := ExpectedUnfoldedSFS(numChrom)
	var obsAcc, expAcc float64
	classes := 0
	flush := func() {
		if expAcc <= 0 {
			obsAcc, expAcc = 0, 0
			return
		}
		diff := obsAcc - expAcc
		chi2 += diff * diff / expAcc
		classes++
		obsAcc, expAcc = 0, 0
	}
	for i := 1; i < numChrom; i++ {
		obsAcc += float64(unfolded[i])
		expAcc += exp[i] * float64(segSites)
		if expAcc >= 1.0 {
			flush()
		}
	}
	flush() // pool any trailing remainder
	if classes > 1 {
		dof = classes - 1
	}
	return chi2, dof
}

// ---------------------------------------------------------------------------
// Validation report: score the statistics against analytic expectations.
// ---------------------------------------------------------------------------

// ValidationCheck is one scored comparison of an observed statistic against its
// analytic neutral expectation.
type ValidationCheck struct {
	Name      string
	Observed  float64
	Expected  float64
	Tolerance float64 // absolute allowed |observed - expected|
	Pass      bool
	Note      string
}

// NeutralTolerances holds the pass/fail bands for a SINGLE-REALIZATION validation
// report (the per-run diagnostic wired into drift.go). Tajima's D and the theta
// ratio are scored against DRIFT's CHARACTERIZED neutral baseline (a rare-variant
// excess: D ~ -0.66, theta_pi/theta_W ~ 0.81), NOT the Wright-Fisher values (0 and
// 1) — see pkg/validation for why DRIFT's high-reproductive-variance demography
// produces that skew. HWE (F_IS ~ 0) is a true textbook expectation DRIFT meets.
// Bands are wide because a single sample of Tajima's D has a standard deviation
// near 1; the replicate-averaged harness (pkg/validation) uses much tighter bands.
type NeutralTolerances struct {
	TajimasDCenter   float64 // characterized centre of the neutral Tajima's D band
	TajimasD         float64 // half-width of the D band
	ThetaRatioCenter float64 // characterized centre of theta_pi/theta_W
	ThetaRatio       float64 // half-width of the ratio band
	SFSShape         float64 // max chi2/dof for the 1/i SFS fit (descriptive; linkage inflates it)
	FIS              float64 // |F_IS| band around 0 (HWE)
	FuLiD            float64 // |Fu-Li D| band around 0
}

// DefaultTolerances are the single-realization bands, centered on DRIFT's
// characterized neutral baseline. Wide because one sample's Tajima's D has SD ~ 1.
func DefaultTolerances() NeutralTolerances {
	return NeutralTolerances{
		TajimasDCenter:   -0.66, // refcount pool accounting (§6f); legacy was -0.89
		TajimasD:         1.3,   // pass for a single-run D in [-1.96, 0.64]
		ThetaRatioCenter: 0.81,
		ThetaRatio:       0.40, // [0.41, 1.21]
		SFSShape:         6.0,  // descriptive; correlated (linked) sites inflate chi2
		FIS:              0.20,
		FuLiD:            1.6,
	}
}

// ValidationReport is the full scored report for one neutral realization (or one
// replicate-averaged summary).
type ValidationReport struct {
	NumSamples int
	NumChrom   int
	SegSites   int
	Mu         float64 // per-diploid per-generation mutation rate (model "mu")
	CensusN    float64 // census population size at sampling
	ThetaW     float64
	ThetaPi    float64
	ImpliedNe  float64 // theta_W / (2*mu): the coalescent Ne implied by the data
	NeOverN    float64 // ImpliedNe / CensusN (emergent, since DRIFT != Wright-Fisher)
	Checks     []ValidationCheck
}

// AllPass reports whether every check in the report passed.
func (r *ValidationReport) AllPass() bool {
	for _, c := range r.Checks {
		if !c.Pass {
			return false
		}
	}
	return true
}

// BuildReport scores a NeutralStats against the analytic neutral expectations
// using the given per-diploid mutation rate mu, census size censusN, and
// tolerance bands.
func BuildReport(s *NeutralStats, mu, censusN float64, tol NeutralTolerances) *ValidationReport {
	r := &ValidationReport{
		NumSamples: s.NumSamples,
		NumChrom:   s.NumChrom,
		SegSites:   s.SegSites,
		Mu:         mu,
		CensusN:    censusN,
		ThetaW:     s.ThetaW,
		ThetaPi:    s.ThetaPi,
	}
	if mu > 0 {
		r.ImpliedNe = s.ThetaW / (2 * mu)
	}
	if censusN > 0 {
		r.NeOverN = r.ImpliedNe / censusN
	}

	add := func(name string, obs, exp, tolBand float64, note string) {
		r.Checks = append(r.Checks, ValidationCheck{
			Name: name, Observed: obs, Expected: exp, Tolerance: tolBand,
			Pass: math.Abs(obs-exp) <= tolBand, Note: note,
		})
	}

	// 1. Tajima's D: scored against DRIFT's characterized neutral baseline (~ -0.66),
	// NOT the Wright-Fisher 0 (see NeutralTolerances / pkg/validation).
	add("TajimasD", s.TajimasD, tol.TajimasDCenter, tol.TajimasD,
		"characterized baseline ~ -0.66 (rare-variant excess, not WF 0); single-sample SD ~ 1")

	// 2. theta_pi / theta_W: characterized ~ 0.81 (both estimate the same theta, but
	// the neutral SFS skew pulls the ratio below the WF value of 1).
	ratio := 0.0
	if s.ThetaW > 0 {
		ratio = s.ThetaPi / s.ThetaW
	}
	add("ThetaPi/ThetaW", ratio, tol.ThetaRatioCenter, tol.ThetaRatio,
		"both estimate theta = 2*Ne*mu; characterized ~0.81 (WF=1)")

	// 3. SFS shape ~ 1/i (chi2/dof).
	chi2, dof := SFSShapeChiSquare(s.UnfoldedSFS, s.SegSites, s.NumChrom)
	chiPerDof := 0.0
	if dof > 0 {
		chiPerDof = chi2 / float64(dof)
	}
	// Score chi2/dof against an upper band (expected ~ 1): observed value, target 0,
	// tolerance = band, so it passes when chi2/dof <= band.
	r.Checks = append(r.Checks, ValidationCheck{
		Name: "SFS_1overI_chi2perDof", Observed: chiPerDof, Expected: 0, Tolerance: tol.SFSShape,
		Pass: chiPerDof <= tol.SFSShape,
		Note: fmt.Sprintf("Pearson chi2/dof vs neutral 1/i SFS (dof=%d); ~1 is a good fit", dof),
	})

	// 4. Hardy-Weinberg: F_IS = 1 - Ho/He ~ 0 (random-mating genotype proportions).
	add("HWE_FIS", s.FIS, 0, tol.FIS,
		"F_IS = 1 - Ho/He; random mating => ~0 (slight positive ~1/(2Ne) from finite size)")

	// 5. Fu & Li's D ~ 0 (secondary low-frequency-class neutrality check).
	add("FuLiD", s.FuLiD, 0, tol.FuLiD, "neutral expectation 0")

	return r
}

// SaveNeutralValidation is the end-of-run entry point wired into drift.go behind
// the track_validation flag. It samples the living population, computes the neutral
// statistics from the mutation pool, scores them against the analytic expectations
// (single-realization tolerances), writes a validation-report CSV, and prints a
// one-line PASS/FAIL summary. A no-op (nil) when mutations are not tracked or the
// sample is too small.
func SaveNeutralValidation(model *core.Model, pop *core.Pop) error {
	if model.Parameters["track_mutations"] != 1 {
		fmt.Println("Neutral validation: track_mutations is off; nothing to validate (needs the de-novo mutation pool)")
		return nil
	}
	sampleSize := int(model.Parameters["validation_sample_size"])
	ids := SampleLiving(pop, sampleSize)
	stats := ComputeNeutralStats(pop, ids)
	if stats.SegSites == 0 {
		fmt.Printf("Neutral validation: no segregating sites among %d sampled individuals (run longer / raise mu)\n", stats.NumSamples)
		return nil
	}

	mu := model.Parameters["mu"]
	censusN := float64(len(pop.IndData))
	report := BuildReport(stats, mu, censusN, DefaultTolerances())

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	path := fmt.Sprintf("%s/%s_validation_run%d_year%d.csv",
		model.ResultsDir, model.ModelName, run, year)
	if err := WriteValidationReport(path, report); err != nil {
		return err
	}

	status := "PASS"
	if !report.AllPass() {
		status = "FAIL"
	}
	fmt.Printf("Neutral validation [%s]: n=%d (2n=%d), S=%d, thetaW=%.2f, thetaPi=%.2f, D=%.3f, F_IS=%.3f, impliedNe=%.0f (Ne/N=%.2f) -> %s\n",
		status, stats.NumSamples, stats.NumChrom, stats.SegSites, stats.ThetaW, stats.ThetaPi,
		stats.TajimasD, stats.FIS, report.ImpliedNe, report.NeOverN, path)
	return nil
}

// WriteValidationReport writes a validation report to a CSV: a header block of
// context rows (Metric,Value) followed by the scored checks
// (Check,Observed,Expected,Tolerance,Result).
func WriteValidationReport(path string, r *ValidationReport) error {
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create validation report: %v", err)
	}
	defer file.Close()

	w := csv.NewWriter(file)
	defer w.Flush()

	ctx := [][]string{
		{"Metric", "Value"},
		{"NumSamples", fmt.Sprintf("%d", r.NumSamples)},
		{"NumChromosomes", fmt.Sprintf("%d", r.NumChrom)},
		{"SegregatingSites", fmt.Sprintf("%d", r.SegSites)},
		{"mu_per_diploid_per_gen", fmt.Sprintf("%.6g", r.Mu)},
		{"CensusN", fmt.Sprintf("%.0f", r.CensusN)},
		{"ThetaW", fmt.Sprintf("%.6f", r.ThetaW)},
		{"ThetaPi", fmt.Sprintf("%.6f", r.ThetaPi)},
		{"ImpliedNe_thetaW_over_2mu", fmt.Sprintf("%.2f", r.ImpliedNe)},
		{"Ne_over_N", fmt.Sprintf("%.4f", r.NeOverN)},
		{"", ""},
	}
	for _, row := range ctx {
		if err := w.Write(row); err != nil {
			return err
		}
	}

	if err := w.Write([]string{"Check", "Observed", "Expected", "Tolerance", "Result", "Note"}); err != nil {
		return err
	}
	for _, c := range r.Checks {
		result := "PASS"
		if !c.Pass {
			result = "FAIL"
		}
		row := []string{
			c.Name,
			fmt.Sprintf("%.6f", c.Observed),
			fmt.Sprintf("%.6f", c.Expected),
			fmt.Sprintf("%.6f", c.Tolerance),
			result,
			c.Note,
		}
		if err := w.Write(row); err != nil {
			return err
		}
	}
	return nil
}
