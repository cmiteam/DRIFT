package analysis

// Molecular-clock / dating layer (roadmap §6e, "Coalescence -> calendar dates under
// explicit, swappable mutation-rate / generation-time assumptions"). The famous
// human dates (Y-Adam / Mt-Eve ~200 kya, out-of-Africa ~50-70 kya) are molecular-
// clock products: an observed genetic divergence D is turned into a time by
//     T_generations = D / (2 * mu)          (mu = per-lineage-per-generation rate)
//     T_years       = T_generations * generation_time.
// Both mu and generation_time are ASSUMED, and both are unsettled in the field
// itself — the direct pedigree rate (~1.2e-8/site/gen) is about half the
// phylogenetic rate, and estimates of generation time range ~22-33 years. Halving
// the assumed mu doubles every date. This layer makes that dependence quantitative
// instead of rhetorical.
//
// DRIFT'S LEVERAGE. DRIFT knows all three quantities a real analyst must assume: the
// true mu it used, the true generation_time, and — through the recorded genealogy —
// the true coalescence dates (Y-Adam / Mt-Eve). So it can date the SAME divergence
// under a grid of assumptions and show exactly how far the inferred date swings from
// the self-consistent (true-mu, true-g) value, with the true genealogical dates
// reported alongside for scale.
//
// SIGNAL (chosen with Rob, over Y-locus-specific dating): the observed divergence is
// the AUTOSOMAL mean pairwise difference theta_pi from the de-novo mutation pool
// (reusing the §6h neutral machinery — the mutation pool is the equilibrium object;
// the founder bitfield is not). Under neutrality E[theta_pi] = 2*mu*T_pair, so the
// molecular pairwise-coalescence time is T_pair = theta_pi / (2*mu) generations. The
// true genealogical Y-Adam / Mt-Eve years-back (FindYAdam / FindMtEve) are reported
// as the ground-truth dates DRIFT actually knows.
//
// WIRING. track_dating (default off) gates it; requires track_mutations (the
// molecular signal). Read-only, end-of-run, no RNG on the default full-sample path
// (dating_sample_size 0). Emits <model>_dating_* (point estimate + true anchors +
// swing) and <model>_dating_sensitivity_* (the mu x g matrix).

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// Default sensitivity grids. Mu factors bracket the pedigree-vs-phylogenetic
// factor-of-two in each direction; generation times span the human range.
var defaultMuFactors = []float64{0.5, 1, 2}
var defaultGenTimes = []float64{20, 25, 29, 35}

// DatingCell is one (assumed mu, assumed generation time) entry of the sensitivity
// matrix: the calendar date the analyst would infer from the observed divergence.
type DatingCell struct {
	MuFactor  float64 // assumed mu as a multiple of the true mu
	Mu        float64 // assumed per-diploid-per-generation rate (TrueMu * MuFactor)
	GenTime   float64 // assumed generation time (years)
	TMRCAGen  float64 // theta_pi / (2*Mu): pairwise coalescence time in generations
	DateYears float64 // TMRCAGen * GenTime
}

// DatingResult is the output of a dating pass.
type DatingResult struct {
	NumSamples int
	SegSites   int
	ThetaPi    float64 // observed autosomal mean pairwise difference (the divergence)

	TrueMu         float64 // the mu DRIFT actually used
	TrueGenTime    float64 // the generation_time DRIFT actually used
	PointTMRCAGen  float64 // theta_pi / (2*TrueMu): self-consistent coalescence time
	PointDateYears float64 // PointTMRCAGen * TrueGenTime: the self-consistent date

	YAdamYearsBack int // true genealogical years back to Y-Adam (-1 if not found)
	MtEveYearsBack int // true genealogical years back to Mt-Eve (-1 if not found)

	Cells       []DatingCell
	MinDate     float64
	MaxDate     float64
	SwingFactor float64 // MaxDate / MinDate across the grid (the headline)
}

// parseFloatList parses a ';'/','/whitespace-separated list of floats, falling back
// to def when the string is empty or yields no valid numbers. Positive values only
// (a mutation-rate factor or generation time is positive); non-positive or malformed
// tokens are skipped.
func parseFloatList(s string, def []float64) []float64 {
	fields := strings.FieldsFunc(s, func(r rune) bool {
		return r == ';' || r == ',' || r == ' ' || r == '\t'
	})
	out := make([]float64, 0, len(fields))
	for _, f := range fields {
		v, err := strconv.ParseFloat(strings.TrimSpace(f), 64)
		if err == nil && v > 0 {
			out = append(out, v)
		}
	}
	if len(out) == 0 {
		return def
	}
	return out
}

// ComputeDating dates the sampled population's autosomal coalescence under the true
// (mu, generation_time) and over the configured sensitivity grid, and attaches the
// true genealogical Y-Adam / Mt-Eve dates. Returns nil when mutations are not tracked
// or there is no divergence to date.
func ComputeDating(model *core.Model, pop *core.Pop, ids []int) *DatingResult {
	stats := ComputeNeutralStats(pop, ids)
	if stats.SegSites == 0 || stats.ThetaPi <= 0 {
		return nil
	}

	trueMu := model.Parameters["mu"]
	trueG := model.Parameters["generation_time"]
	if trueG <= 0 {
		trueG = 1
	}

	res := &DatingResult{
		NumSamples:     stats.NumSamples,
		SegSites:       stats.SegSites,
		ThetaPi:        stats.ThetaPi,
		TrueMu:         trueMu,
		TrueGenTime:    trueG,
		YAdamYearsBack: -1,
		MtEveYearsBack: -1,
	}
	if trueMu > 0 {
		res.PointTMRCAGen = stats.ThetaPi / (2 * trueMu)
		res.PointDateYears = res.PointTMRCAGen * trueG
	}

	// True genealogical anchors (the dates DRIFT actually knows). These use the
	// recorded paternal/maternal databases; when genealogy is not tracked they
	// report "not found" and stay at -1.
	if ya := FindYAdam(model, pop); ya.FoundYAdam {
		res.YAdamYearsBack = ya.GenerationsBack // computed as year - birthYear, i.e. calendar years
	}
	if me := FindMtEve(model, pop); me.FoundMtEve {
		res.MtEveYearsBack = me.GenerationsBack
	}

	// Sensitivity matrix over assumed mu x generation time.
	muFactors := parseFloatList(model.StringParam("dating_mu_factors", ""), defaultMuFactors)
	genTimes := parseFloatList(model.StringParam("dating_gen_times", ""), defaultGenTimes)
	first := true
	for _, mf := range muFactors {
		mu := trueMu * mf
		if mu <= 0 {
			continue
		}
		tmrcaGen := stats.ThetaPi / (2 * mu)
		for _, g := range genTimes {
			date := tmrcaGen * g
			res.Cells = append(res.Cells, DatingCell{
				MuFactor: mf, Mu: mu, GenTime: g, TMRCAGen: tmrcaGen, DateYears: date,
			})
			if first || date < res.MinDate {
				res.MinDate = date
			}
			if first || date > res.MaxDate {
				res.MaxDate = date
			}
			first = false
		}
	}
	if res.MinDate > 0 {
		res.SwingFactor = res.MaxDate / res.MinDate
	}
	return res
}

// SaveDating writes the point-estimate/anchor summary and the sensitivity matrix to
// two CSVs in the model's results directory and prints a one-line summary. A no-op
// (nil) when the result is nil (nothing to date).
func SaveDating(model *core.Model, result *DatingResult) error {
	if result == nil {
		fmt.Println("Dating (§6e): no divergence to date (needs track_mutations and segregating mutations)")
		return nil
	}

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]

	// --- Point estimate + true anchors + swing summary ---
	sumPath := fmt.Sprintf("%s/%s_dating_run%d_year%d.csv", model.ResultsDir, model.ModelName, run, year)
	sf, err := os.Create(sumPath)
	if err != nil {
		return fmt.Errorf("failed to create dating summary file: %v", err)
	}
	defer sf.Close()
	sw := csv.NewWriter(sf)
	defer sw.Flush()

	if err := sw.Write([]string{"Quantity", "Value", "Units", "Note"}); err != nil {
		return err
	}
	rows := [][]string{
		{"ObservedThetaPi", fmt.Sprintf("%.4f", result.ThetaPi), "mutations", "autosomal mean pairwise difference (mutation pool)"},
		{"SegSites", fmt.Sprintf("%d", result.SegSites), "sites", "segregating mutation lineages in the sample"},
		{"NumSamples", fmt.Sprintf("%d", result.NumSamples), "diploids", "sampled individuals"},
		{"TrueMu", fmt.Sprintf("%.6g", result.TrueMu), "per diploid/gen", "the mutation rate DRIFT used"},
		{"TrueGenerationTime", fmt.Sprintf("%.4g", result.TrueGenTime), "years", "the generation time DRIFT used"},
		{"PointTMRCA_generations", fmt.Sprintf("%.2f", result.PointTMRCAGen), "generations", "theta_pi / (2*mu) at true mu"},
		{"PointTMRCA_years", fmt.Sprintf("%.1f", result.PointDateYears), "years", "self-consistent molecular date (true mu, true g)"},
		{"TrueYAdam_yearsBack", yearsBackStr(result.YAdamYearsBack), "years", "genealogical ground truth (paternal MRCA)"},
		{"TrueMtEve_yearsBack", yearsBackStr(result.MtEveYearsBack), "years", "genealogical ground truth (maternal MRCA)"},
		{"MinDate_years", fmt.Sprintf("%.1f", result.MinDate), "years", "youngest date across the mu x g grid"},
		{"MaxDate_years", fmt.Sprintf("%.1f", result.MaxDate), "years", "oldest date across the mu x g grid"},
		{"SwingFactor", fmt.Sprintf("%.2f", result.SwingFactor), "x", "MaxDate / MinDate: how far assumptions move the date"},
	}
	for _, r := range rows {
		if err := sw.Write(r); err != nil {
			return err
		}
	}

	// --- Sensitivity matrix ---
	matPath := fmt.Sprintf("%s/%s_dating_sensitivity_run%d_year%d.csv", model.ResultsDir, model.ModelName, run, year)
	mf, err := os.Create(matPath)
	if err != nil {
		return fmt.Errorf("failed to create dating sensitivity file: %v", err)
	}
	defer mf.Close()
	mw := csv.NewWriter(mf)
	defer mw.Flush()

	if err := mw.Write([]string{"MuFactor", "AssumedMu", "GenerationTime", "TMRCA_generations", "Date_years"}); err != nil {
		return err
	}
	for _, c := range result.Cells {
		row := []string{
			fmt.Sprintf("%.4g", c.MuFactor),
			fmt.Sprintf("%.6g", c.Mu),
			fmt.Sprintf("%.4g", c.GenTime),
			fmt.Sprintf("%.2f", c.TMRCAGen),
			fmt.Sprintf("%.1f", c.DateYears),
		}
		if err := mw.Write(row); err != nil {
			return err
		}
	}

	fmt.Printf("Dating (§6e): theta_pi=%.3f -> point TMRCA %.0f y (true mu/g); grid spans %.0f-%.0f y (%.2fx swing); true Y-Adam=%s, Mt-Eve=%s -> %s\n",
		result.ThetaPi, result.PointDateYears, result.MinDate, result.MaxDate, result.SwingFactor,
		yearsBackStr(result.YAdamYearsBack), yearsBackStr(result.MtEveYearsBack), sumPath)
	return nil
}

// yearsBackStr formats a genealogical years-back value, or "n/a" when not found (-1).
func yearsBackStr(v int) string {
	if v < 0 {
		return "n/a"
	}
	return strconv.Itoa(v)
}
