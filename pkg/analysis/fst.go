package analysis

// Fst — population differentiation (roadmap §6b, "Differentiation & admixture
// statistics"). Fst quantifies how much of the total genetic variation is due to
// differences BETWEEN demes rather than within them; the out-of-Africa story is
// fundamentally an argument about continental differentiation, so Fst is a
// non-negotiable part of the modern OoA toolkit.
//
// POPULATIONS IN DRIFT. DRIFT is individual-based and (until now) single-
// population. "Populations" for Fst are the discrete DEMES introduced with the
// island-model infrastructure (num_demes / deme_migration_rate): every individual
// carries a Deme label (individual.Deme), founders are partitioned across demes,
// the label is inherited maternally, and mating is partitioned by deme. With
// num_demes = 1 there is a single deme and Fst is undefined (returned as an empty
// result) — you must run with num_demes >= 2 to get differentiation.
//
// ESTIMATORS. Two standard estimators are computed, both as a RATIO OF SUMS across
// sites (Bhatia et al. 2013 "ratio of averages" — the numerator and denominator are
// each summed over sites and then divided, which is far less biased than averaging
// per-site ratios):
//
//   - Hudson's Fst (pairwise). Per site, for demes i,j with sample allele
//     frequencies p_i,p_j and haploid sample sizes n_i,n_j:
//         num = (p_i - p_j)^2 - p_i(1-p_i)/(n_i-1) - p_j(1-p_j)/(n_j-1)
//         den = p_i(1-p_j) + p_j(1-p_i)
//     Fst_Hudson = Σ num / Σ den.
//
//   - Weir & Cockerham 1984 θ (pairwise AND global over all demes). The standard
//     a/(a+b+c) variance-components estimator using diploid sample sizes, allele
//     frequencies, and OBSERVED heterozygosities (a heterozygote is an individual
//     whose two phased strand copies differ at the site). θ = Σa / Σ(a+b+c).
//
// A site contributes only if it is polymorphic across the eligible demes. Demes
// with fewer than 2 sampled diploids are dropped (the estimators are undefined for
// them) with a logged note.
//
// Wired into drift.go behind the track_Fst parameter (default off); requires
// track_DNA so chromosomes exist, and num_demes >= 2 so there is structure to
// measure. Reuses SampleIDs, buildContigs and bitAt (vcf.go).

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"encoding/csv"
	"fmt"
	"os"
	"sort"
)

// FstPair holds both estimators for one unordered pair of demes.
type FstPair struct {
	DemeA, DemeB int
	Hudson       float64 // Hudson's Fst
	WC           float64 // Weir & Cockerham θ (r=2)
	NumSites     int     // polymorphic sites used for this pair
}

// FstResult is the output of a differentiation pass.
type FstResult struct {
	Demes     []int       // eligible demes (>= 2 sampled diploids), sorted
	DemeSizes map[int]int // sampled diploids per eligible deme
	Pairs     []FstPair   // one entry per unordered deme pair
	GlobalWC  float64     // Weir & Cockerham θ across all eligible demes
	NumSites  int         // sites polymorphic across the eligible demes
}

// demeMembers groups the sampled ids that carry two strand copies by their Deme
// label. The returned slices are in the (already sorted) input id order.
func demeMembers(pop *core.Pop, ids []int) map[int][]int {
	byDeme := make(map[int][]int)
	for _, id := range ids {
		c, ok := pop.Chromosomes[id]
		if !ok || len(c) < 2 || c[0] == nil || c[1] == nil {
			continue
		}
		d := pop.IndData[id][individual.Deme]
		byDeme[d] = append(byDeme[d], id)
	}
	return byDeme
}

// siteCounts holds, for one site and one deme, the derived-allele count and the
// heterozygote count over that deme's sampled individuals. The diploid sample size
// is constant per deme (all members carry chromosomes), so it lives in FstResult.
type siteCounts struct {
	derived int // number of 1-bits across both strand copies
	het     int // number of individuals whose two copies differ
}

// countSite tallies derived alleles and heterozygotes for one deme's members at a
// single genome bit.
func countSite(pop *core.Pop, members []int, bit int) siteCounts {
	var sc siteCounts
	for _, id := range members {
		c := pop.Chromosomes[id]
		a0 := bitAt(c[0], bit)
		a1 := bitAt(c[1], bit)
		if a0 {
			sc.derived++
		}
		if a1 {
			sc.derived++
		}
		if a0 != a1 {
			sc.het++
		}
	}
	return sc
}

// ComputeFst runs the differentiation pass over the given sampled individual IDs
// (see SampleIDs), grouping them into demes by their Deme label. Returns an empty
// result (NumSites 0, no pairs) when there are fewer than two eligible demes.
func ComputeFst(model *core.Model, pop *core.Pop, ids []int) *FstResult {
	return computeFstFromPartition(model, pop, demeMembers(pop, ids), "deme")
}

// computeFstFromPartition is the estimator core shared by the deme-label Fst
// (ComputeFst) and the geographic-region Fst (ComputeGeographicFst, geographic_fst.go).
// It takes an already-built partition of sampled ids into groups (deme labels or
// geographic-region ids) and computes Hudson's and Weir & Cockerham's estimators over
// them. groupLabel is only used in the informational log lines ("deme"/"region"). The
// FstResult's Demes/Pairs carry whatever integer keys the partition used. Returns an
// empty result (NumSites 0, no pairs) when there are fewer than two eligible groups.
func computeFstFromPartition(model *core.Model, pop *core.Pop, byDeme map[int][]int, groupLabel string) *FstResult {
	// Eligible groups have at least 2 sampled diploids (both estimators are
	// undefined below that). Keep them sorted for deterministic output.
	demes := make([]int, 0, len(byDeme))
	for d, members := range byDeme {
		if len(members) >= 2 {
			demes = append(demes, d)
		} else if len(members) > 0 {
			fmt.Printf("Fst: %s %d has only %d sampled diploid(s); excluded (need >= 2)\n", groupLabel, d, len(members))
		}
	}
	sort.Ints(demes)

	res := &FstResult{Demes: demes, DemeSizes: map[int]int{}}
	for _, d := range demes {
		res.DemeSizes[d] = len(byDeme[d])
	}
	if len(demes) < 2 {
		fmt.Printf("Fst: %d eligible %s(s); need >= 2 (run long enough for each to persist)\n", len(demes), groupLabel)
		return res
	}

	// Pairwise accumulators (ratio-of-sums), keyed by canonical (low,high) pair.
	type pairAcc struct {
		hudNum, hudDen float64
		wcA, wcB, wcC  float64
		sites          int
	}
	pairs := map[[2]int]*pairAcc{}
	for i := 0; i < len(demes); i++ {
		for j := i + 1; j < len(demes); j++ {
			pairs[[2]int{demes[i], demes[j]}] = &pairAcc{}
		}
	}
	// Global Weir & Cockerham accumulators (across all eligible demes).
	var gA, gB, gC float64

	// Walk every genome position via the chromosome arms (same coordinate scan as
	// the VCF export), tallying per-deme counts once per site.
	contigs := buildContigs(model)
	counts := make(map[int]siteCounts, len(demes)) // reused per site
	for _, ct := range contigs {
		for _, arm := range ct.arms {
			start, length := arm[0], arm[1]
			for bit := start; bit < start+length; bit++ {
				totalDerived, totalAlleles := 0, 0
				for _, d := range demes {
					sc := countSite(pop, byDeme[d], bit)
					counts[d] = sc
					totalDerived += sc.derived
					totalAlleles += 2 * len(byDeme[d])
				}
				// Skip sites monomorphic across the eligible demes.
				if totalDerived == 0 || totalDerived == totalAlleles {
					continue
				}
				res.NumSites++

				// Pairwise estimators.
				for i := 0; i < len(demes); i++ {
					for j := i + 1; j < len(demes); j++ {
						di, dj := demes[i], demes[j]
						acc := pairs[[2]int{di, dj}]
						num, den := hudsonSite(counts[di], len(byDeme[di]), counts[dj], len(byDeme[dj]))
						acc.hudNum += num
						acc.hudDen += den
						a, b, c := wcSite([]siteCounts{counts[di], counts[dj]}, []int{len(byDeme[di]), len(byDeme[dj])})
						acc.wcA += a
						acc.wcB += b
						acc.wcC += c
						acc.sites++
					}
				}

				// Global Weir & Cockerham across all eligible demes.
				scs := make([]siteCounts, len(demes))
				ns := make([]int, len(demes))
				for k, d := range demes {
					scs[k] = counts[d]
					ns[k] = len(byDeme[d])
				}
				a, b, c := wcSite(scs, ns)
				gA += a
				gB += b
				gC += c
			}
		}
	}

	// Finalize pairwise ratios.
	for i := 0; i < len(demes); i++ {
		for j := i + 1; j < len(demes); j++ {
			di, dj := demes[i], demes[j]
			acc := pairs[[2]int{di, dj}]
			res.Pairs = append(res.Pairs, FstPair{
				DemeA:    di,
				DemeB:    dj,
				Hudson:   ratio(acc.hudNum, acc.hudDen),
				WC:       ratio(acc.wcA, acc.wcA+acc.wcB+acc.wcC),
				NumSites: acc.sites,
			})
		}
	}
	res.GlobalWC = ratio(gA, gA+gB+gC)
	return res
}

// hudsonSite returns the per-site Hudson numerator and denominator for two demes.
// n1,n2 are DIPLOID sample sizes; the haploid sizes used in the estimator are 2n.
func hudsonSite(a siteCounts, n1 int, b siteCounts, n2 int) (num, den float64) {
	hn1 := float64(2 * n1) // haploid sample sizes
	hn2 := float64(2 * n2)
	p1 := float64(a.derived) / hn1
	p2 := float64(b.derived) / hn2
	diff := p1 - p2
	num = diff*diff - p1*(1-p1)/(hn1-1) - p2*(1-p2)/(hn2-1)
	den = p1*(1-p2) + p2*(1-p1)
	return num, den
}

// wcSite returns the per-site Weir & Cockerham (1984) a, b, c variance components
// for r demes, given each deme's site counts and DIPLOID sample size. The caller
// accumulates a,b,c across sites and forms θ = Σa / Σ(a+b+c). Demes are assumed
// to have >= 2 diploids each (guaranteed by the eligibility filter).
func wcSite(scs []siteCounts, ns []int) (a, b, c float64) {
	r := float64(len(scs))

	var sumN, sumN2 float64
	for _, n := range ns {
		sumN += float64(n)
		sumN2 += float64(n) * float64(n)
	}
	nBar := sumN / r
	nC := (sumN - sumN2/sumN) / (r - 1) // sample-size correction

	// Weighted mean allele frequency and mean observed heterozygosity.
	var pBar, hBar float64
	for i := range scs {
		ni := float64(ns[i])
		pi := float64(scs[i].derived) / (2 * ni)
		hi := float64(scs[i].het) / ni
		pBar += ni * pi
		hBar += ni * hi
	}
	pBar /= sumN
	hBar /= sumN

	// Sample variance of allele frequency across demes.
	var s2 float64
	for i := range scs {
		ni := float64(ns[i])
		pi := float64(scs[i].derived) / (2 * ni)
		d := pi - pBar
		s2 += ni * d * d
	}
	s2 /= (r - 1) * nBar

	pq := pBar * (1 - pBar)
	a = (nBar / nC) * (s2 - (1/(nBar-1))*(pq-((r-1)/r)*s2-hBar/4))
	b = (nBar / (nBar - 1)) * (pq - ((r-1)/r)*s2 - ((2*nBar-1)/(4*nBar))*hBar)
	c = hBar / 2
	return a, b, c
}

// ratio returns num/den, or 0 when the denominator is ~0 (an undefined Fst, e.g.
// no variation), so the output stays finite.
func ratio(num, den float64) float64 {
	if den == 0 {
		return 0
	}
	return num / den
}

// SaveFst writes the pairwise Fst table (both estimators) plus a global-θ row to a
// CSV in the model's results directory and prints a one-line summary, matching the
// other analyses' end-of-run reporting. A no-op (nil error) when the result has no
// pairs (fewer than two demes).
func SaveFst(model *core.Model, result *FstResult) error {
	if result == nil || len(result.Pairs) == 0 {
		return nil
	}

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	path := fmt.Sprintf("%s/%s_fst_run%d_year%d.csv",
		model.ResultsDir, model.ModelName, run, year)

	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create Fst file: %v", err)
	}
	defer file.Close()

	w := csv.NewWriter(file)
	defer w.Flush()

	if err := w.Write([]string{"DemeA", "DemeB", "Fst_Hudson", "Fst_WC", "NumSites"}); err != nil {
		return err
	}
	for _, p := range result.Pairs {
		row := []string{
			fmt.Sprintf("%d", p.DemeA),
			fmt.Sprintf("%d", p.DemeB),
			fmt.Sprintf("%.6f", p.Hudson),
			fmt.Sprintf("%.6f", p.WC),
			fmt.Sprintf("%d", p.NumSites),
		}
		if err := w.Write(row); err != nil {
			return err
		}
	}
	// Global Weir & Cockerham θ across all eligible demes (DemeB = -1 sentinel).
	if err := w.Write([]string{
		"ALL", "-1",
		"", // Hudson is pairwise-only
		fmt.Sprintf("%.6f", result.GlobalWC),
		fmt.Sprintf("%d", result.NumSites),
	}); err != nil {
		return err
	}

	fmt.Printf("Fst: %d demes, %d pairs, %d polymorphic sites (global W&C θ = %.4f) → %s\n",
		len(result.Demes), len(result.Pairs), result.NumSites, result.GlobalWC, path)
	return nil
}
