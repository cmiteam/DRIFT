package analysis

// Heterozygosity-vs-distance-from-origin (roadmap §6c, second bullet). The flagship
// out-of-Africa signature is the near-linear DECLINE of within-population
// heterozygosity with distance from an African origin (Ramachandran et al. 2005
// PNAS; Prugnolle et al. 2005 Curr Biol). It is produced by a SERIAL FOUNDER
// EFFECT: each colonization step founds the next population from a small sample of
// the previous one, shedding a little diversity at every step, so heterozygosity
// falls roughly linearly with the number of founder events — and geographic
// distance is just a proxy for that count.
//
// This analysis measures the gradient DRIFT actually produces, so a canonical
// out-of-Africa scenario and a single-origin (Babel) sequential-dispersal scenario
// — run in the SAME engine under the SAME fixed seed — can be compared head to
// head. It computes, per deme:
//
//   - Expected heterozygosity He (Nei's unbiased gene diversity): the mean over
//     sites of (2n/(2n-1))·2p(1-p), with p the deme's sample derived-allele
//     frequency and n its diploid sample size. This is the quantity the published
//     gradients are drawn in.
//   - Observed heterozygosity Ho: the mean over sites of the fraction of sampled
//     individuals whose two phased strand copies differ (reusing the fst.go het
//     tally).
//
// To keep demes on a common footing, both are averaged over the SAME site set —
// the sites polymorphic across the pooled eligible-deme sample (the Fst NumSites
// convention). Each deme is tagged with its serial-founder RANK (distance from the
// origin), obtained from the demographic scheduler when it exposes DemeRanks (the
// demography.Scenario does); with no scheduler the deme label is used as a rank
// proxy and HasRanks is false. Finally an ordinary-least-squares line He ~ rank is
// fitted across demes: its SLOPE is the head-to-head number — how fast diversity
// decays per founder step — with R² for how linear the decline is.
//
// Wired into drift.go behind track_het_distance (default off); requires track_DNA
// (chromosomes) and num_demes/a scenario that yields >= 2 eligible demes. Reuses
// demeMembers/countSite (fst.go) and buildContigs/bitAt (vcf.go).

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"os"
	"sort"
)

// demeRanker is the optional capability a demographic scheduler advertises to
// supply serial-founder ranks. Kept local (rather than widening
// core.DemographyScheduler) so ordinary schedulers and test mocks need not
// implement it. pkg/demography.Scenario satisfies it.
type demeRanker interface {
	DemeRanks() map[int]int
}

// DemeHet is one deme's heterozygosity summary at its distance-from-origin rank.
type DemeHet struct {
	Deme            int
	Rank            int     // serial-founder steps from origin (or deme label if unknown)
	SampledDiploids int     // eligible diploids sampled from this deme
	ExpHet          float64 // Nei unbiased gene diversity, mean over sites
	ObsHet          float64 // observed heterozygosity, mean over sites
	NumSites        int     // sites used (pooled-polymorphic)
}

// HetDistanceResult is the output of a heterozygosity-gradient pass.
type HetDistanceResult struct {
	Demes     []DemeHet // one per eligible deme, sorted by (rank, deme)
	NumSites  int       // pooled-polymorphic sites (shared denominator)
	HasRanks  bool      // true when real serial-founder ranks were available
	Slope     float64   // OLS slope of He on rank (the gradient; <0 = decline)
	Intercept float64   // OLS intercept (He at the origin)
	R2        float64   // coefficient of determination of the He~rank fit
}

// ComputeHetDistance runs the gradient pass over the given sampled individual IDs
// (see SampleIDs), grouping them into demes by their Deme label. Returns an empty
// result (NumSites 0) when fewer than two demes have >= 2 sampled diploids.
func ComputeHetDistance(model *core.Model, pop *core.Pop, ids []int) *HetDistanceResult {
	byDeme := demeMembers(pop, ids)

	// Eligible demes have >= 2 sampled diploids (the unbiased estimators divide by
	// 2n-1 and need at least one heterozygosity observation). Keep sorted.
	demes := make([]int, 0, len(byDeme))
	for d, members := range byDeme {
		if len(members) >= 2 {
			demes = append(demes, d)
		} else if len(members) > 0 {
			fmt.Printf("HetDistance: deme %d has only %d sampled diploid(s); excluded (need >= 2)\n", d, len(members))
		}
	}
	sort.Ints(demes)

	res := &HetDistanceResult{}
	if len(demes) < 2 {
		fmt.Printf("HetDistance: %d eligible deme(s); need >= 2 (set num_demes / a scenario and run long enough for demes to persist)\n", len(demes))
		return res
	}

	// Serial-founder ranks from the scheduler when it can supply them; otherwise
	// fall back to the deme label as a rank proxy.
	ranks, hasRanks := demeRanks(model, demes)
	res.HasRanks = hasRanks

	// Per-deme accumulators over the shared (pooled-polymorphic) site set.
	sumExp := make(map[int]float64, len(demes))
	sumObs := make(map[int]float64, len(demes))

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
				// Skip sites monomorphic across the pooled eligible-deme sample.
				if totalDerived == 0 || totalDerived == totalAlleles {
					continue
				}
				res.NumSites++
				for _, d := range demes {
					n := len(byDeme[d])
					sc := counts[d]
					twoN := float64(2 * n)
					p := float64(sc.derived) / twoN
					// Nei unbiased gene diversity for a biallelic locus.
					sumExp[d] += (twoN / (twoN - 1)) * 2 * p * (1 - p)
					sumObs[d] += float64(sc.het) / float64(n)
				}
			}
		}
	}

	for _, d := range demes {
		dh := DemeHet{
			Deme:            d,
			Rank:            ranks[d],
			SampledDiploids: len(byDeme[d]),
			NumSites:        res.NumSites,
		}
		if res.NumSites > 0 {
			dh.ExpHet = sumExp[d] / float64(res.NumSites)
			dh.ObsHet = sumObs[d] / float64(res.NumSites)
		}
		res.Demes = append(res.Demes, dh)
	}
	// Sort by rank, then deme, so the gradient reads in dispersal order.
	sort.Slice(res.Demes, func(i, j int) bool {
		if res.Demes[i].Rank != res.Demes[j].Rank {
			return res.Demes[i].Rank < res.Demes[j].Rank
		}
		return res.Demes[i].Deme < res.Demes[j].Deme
	})

	res.Slope, res.Intercept, res.R2 = fitLine(res.Demes)
	return res
}

// demeRanks returns the serial-founder rank of each eligible deme and whether real
// ranks were available. When the scheduler cannot supply them, each deme's label
// stands in as its rank (correct for scenarios that number demes in dispersal
// order, and harmless otherwise — HasRanks flags the fallback).
func demeRanks(model *core.Model, demes []int) (map[int]int, bool) {
	out := make(map[int]int, len(demes))
	if r, ok := model.DemographyScheduler.(demeRanker); ok && model.DemographyScheduler != nil {
		all := r.DemeRanks()
		missing := false
		for _, d := range demes {
			if rk, present := all[d]; present {
				out[d] = rk
			} else {
				out[d] = d // deme not in the scenario tree; fall back to its label
				missing = true
			}
		}
		return out, !missing
	}
	for _, d := range demes {
		out[d] = d
	}
	return out, false
}

// fitLine returns the ordinary-least-squares slope, intercept and R² of ExpHet on
// Rank across the demes. With fewer than two distinct ranks the fit is undefined
// and all three are 0.
func fitLine(demes []DemeHet) (slope, intercept, r2 float64) {
	n := float64(len(demes))
	if n < 2 {
		return 0, 0, 0
	}
	var sx, sy, sxx, sxy, syy float64
	for _, d := range demes {
		x := float64(d.Rank)
		y := d.ExpHet
		sx += x
		sy += y
		sxx += x * x
		sxy += x * y
		syy += y * y
	}
	denom := n*sxx - sx*sx
	if denom == 0 { // all ranks identical → no gradient defined
		return 0, 0, 0
	}
	slope = (n*sxy - sx*sy) / denom
	intercept = (sy - slope*sx) / n
	ssTot := syy - sy*sy/n
	if ssTot > 0 {
		ssRes := 0.0
		for _, d := range demes {
			pred := intercept + slope*float64(d.Rank)
			e := d.ExpHet - pred
			ssRes += e * e
		}
		r2 = 1 - ssRes/ssTot
	}
	return slope, intercept, r2
}

// SaveHetDistance writes two CSVs in the model's results directory: the per-deme
// heterozygosity-vs-rank table (the gradient) and a one-row gradient summary (the
// head-to-head number — slope/intercept/R²). It also prints a one-line summary,
// matching the other analyses' end-of-run reporting. A no-op (nil) when the result
// has fewer than two demes.
func SaveHetDistance(model *core.Model, result *HetDistanceResult) error {
	if result == nil || len(result.Demes) < 2 {
		return nil
	}

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]

	tablePath := fmt.Sprintf("%s/%s_hetdistance_run%d_year%d.csv",
		model.ResultsDir, model.ModelName, run, year)
	if err := writeHetTable(tablePath, result); err != nil {
		return err
	}

	gradPath := fmt.Sprintf("%s/%s_hetgradient_run%d_year%d.csv",
		model.ResultsDir, model.ModelName, run, year)
	if err := writeHetGradient(gradPath, result); err != nil {
		return err
	}

	note := ""
	if !result.HasRanks {
		note = " (deme-label ranks: no scheduler)"
	}
	fmt.Printf("HetDistance: %d demes, %d polymorphic sites, He~rank slope = %.5f (R² = %.3f)%s → %s\n",
		len(result.Demes), result.NumSites, result.Slope, result.R2, note, gradPath)
	return nil
}

func writeHetTable(path string, result *HetDistanceResult) error {
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create HetDistance file: %v", err)
	}
	defer file.Close()
	w := csv.NewWriter(file)
	defer w.Flush()

	if err := w.Write([]string{"Deme", "Rank", "SampledDiploids", "ExpHet", "ObsHet", "NumSites"}); err != nil {
		return err
	}
	for _, d := range result.Demes {
		if err := w.Write([]string{
			fmt.Sprintf("%d", d.Deme),
			fmt.Sprintf("%d", d.Rank),
			fmt.Sprintf("%d", d.SampledDiploids),
			fmt.Sprintf("%.6f", d.ExpHet),
			fmt.Sprintf("%.6f", d.ObsHet),
			fmt.Sprintf("%d", d.NumSites),
		}); err != nil {
			return err
		}
	}
	return nil
}

func writeHetGradient(path string, result *HetDistanceResult) error {
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create HetGradient file: %v", err)
	}
	defer file.Close()
	w := csv.NewWriter(file)
	defer w.Flush()

	if err := w.Write([]string{"NumDemes", "Slope", "Intercept", "R2", "NumSites", "HasRanks"}); err != nil {
		return err
	}
	return w.Write([]string{
		fmt.Sprintf("%d", len(result.Demes)),
		fmt.Sprintf("%.6f", result.Slope),
		fmt.Sprintf("%.6f", result.Intercept),
		fmt.Sprintf("%.6f", result.R2),
		fmt.Sprintf("%d", result.NumSites),
		fmt.Sprintf("%t", result.HasRanks),
	})
}
