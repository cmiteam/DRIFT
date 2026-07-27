package analysis

// Ne-through-time — effective population size across the run (roadmap §6d,
// "Effective population size & the bottleneck question"). Creationist history
// posits recent severe bottlenecks (Flood, Babel); the field's recent-Ne methods
// (PSMC/MSMC/GONE) are actively debated. Because DRIFT is forward-time and knows
// the true census trajectory, it can show what these estimators recover from a
// known truth — and, in particular, whether a temporal / LD estimate detects a
// bottleneck the way the census does.
//
// TWO ESTIMATORS (this file = the temporal core; ne_ld.go = the caveated LD leg):
//
//   - TEMPORAL (Waples 1989, plan II). The rigorous, calibration-free core. It
//     needs no recombination map (DRIFT has none yet — §6a) and makes NO
//     equilibrium assumption, which matters because the Chromosomes bitfield read
//     here holds only FOUNDER standing variation drifting toward loss/fixation
//     (§6h) — a non-equilibrium object. The temporal method simply measures the
//     realized drift between two time samples: for a locus with sample derived-
//     allele frequencies x (earlier) and y (later),
//         Fc = Σ (x-y)^2 / Σ [ (x+y)/2 - x*y ]        (ratio of sums over loci)
//         Ne = t / ( 2 * ( Fc - 1/n0 - 1/nt ) )
//     where t is the number of generations between the samples and n0, nt are the
//     HAPLOID sample sizes (= 2 * sampled diploids). The ratio-of-sums form (as in
//     fst.go) is less biased than a mean of per-locus ratios. Loci are those
//     segregating in the EARLIER sample and tracked forward (a now-fixed locus
//     correctly contributes y ∈ {0,1}). If the sampling-corrected Fc is <= 0 the
//     drift signal is below sampling noise and Ne is reported as +Inf (unresolved).
//
//   - LD (ne_ld.go). GONE-style recent Ne from r^2, present but explicitly
//     approximate: its absolute calibration needs a recombination map (§6a), so
//     the bit-distance -> recombination-fraction conversion is a single documented
//     knob (ne_cM_per_bit), off by default.
//
// The per-window series contrasts directly with the §6h long-term coalescent Ne
// (ImpliedNe = theta_W / (2*mu), emergent Ne/N ≈ 0.17): the temporal estimator sees
// the RECENT trajectory (including bottlenecks) that the theta-based value averages
// over.
//
// WIRING. track_Ne (default off) gates the whole thing. Capture fires in the main
// run loop at the ne_interval cadence (drift.go); each capture is a deterministic
// full-population scan of the bitfield — it draws NO RNG (ne_sample_size > 0 to
// subsample is opt-in and does draw) — so a track_Ne run's simulation is
// byte-identical to a track_Ne-off run. SaveNeTimeSeries writes the series at
// end-of-run. Reuses buildContigs/bitAt (vcf.go) and individual.Deme (§6b).

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"math"
	"sort"
)

// CaptureNe records one Ne-through-time window (§6d). It compares the current
// allele frequencies against the previous snapshot (temporal estimator), optionally
// adds the LD recent-Ne companion, stores an NeSnapshot in pop.NeHistory, then
// replaces the previous snapshot with the current one. Draws no RNG unless
// ne_sample_size > 0 (subsampling). A no-op the first time it is called (there is no
// earlier sample yet) except to seed the previous snapshot.
func CaptureNe(model *core.Model, pop *core.Pop) {
	if pop.NeHistory == nil {
		pop.NeHistory = map[int]*core.NeSnapshot{}
	}
	year := model.FreeParameters["year"]

	// Which individuals to include. Default (0) = all with tracked chromosomes,
	// deterministic and RNG-free. A positive ne_sample_size subsamples (opt-in;
	// draws RNG, so it perturbs the stream for the rest of the run).
	sampleSize := int(model.Parameters["ne_sample_size"])
	ids := SampleIDs(pop, sampleSize)

	// Group current members: -1 = pooled, deme ids when the island model is active.
	members := groupMembers(model, pop, ids)

	// Temporal estimate against the previous snapshot (if any).
	snap := &core.NeSnapshot{
		Year:    year,
		CensusN: len(pop.IndData),
	}
	if prev := pop.NePrevSnap; prev != nil {
		genTime := model.Parameters["generation_time"]
		if genTime <= 0 {
			genTime = 1 // fall back to years if generation time is unset
		}
		tGen := float64(year-prev.Year) / genTime
		snap.IntervalGen = tGen

		// Pooled temporal Ne.
		if pooled, ok := members[-1]; ok {
			ne, nsites := temporalNe(pop, prev, -1, pooled, tGen)
			snap.TemporalNe = ne
			snap.TemporalNumSites = nsites
		}
		// Per-deme temporal Ne (island model only).
		if len(members) > 1 {
			snap.DemeNe = map[int]float64{}
			demes := make([]int, 0, len(members))
			for g := range members {
				if g >= 0 {
					demes = append(demes, g)
				}
			}
			sort.Ints(demes)
			for _, g := range demes {
				ne, _ := temporalNe(pop, prev, g, members[g], tGen)
				snap.DemeNe[g] = ne
			}
		}

		// LD-based recent-Ne companion (approximate; disabled unless ne_cM_per_bit > 0).
		snap.LDNe, snap.LDNumPairs = ldRecentNe(model, pop, members[-1])
	}
	pop.NeHistory[year] = snap

	// The current window becomes the reference for the next capture.
	pop.NePrevSnap = snapshotFreqs(model, pop, members, year)
}

// groupMembers partitions the sampled ids into estimation groups: key -1 is the
// pooled whole-population sample; keys >= 0 are demes (only when num_demes > 1).
// Only individuals carrying two strand copies are included. Members preserve the
// sorted input id order (SampleIDs returns sorted ids), so all downstream sums are
// order-deterministic.
func groupMembers(model *core.Model, pop *core.Pop, ids []int) map[int][]int {
	numDemes := int(model.Parameters["num_demes"])
	groups := map[int][]int{}
	for _, id := range ids {
		c, ok := pop.Chromosomes[id]
		if !ok || len(c) < 2 || c[0] == nil || c[1] == nil {
			continue
		}
		groups[-1] = append(groups[-1], id)
		if numDemes > 1 {
			d := pop.IndData[id][individual.Deme]
			groups[d] = append(groups[d], id)
		}
	}
	return groups
}

// snapshotFreqs builds a FreqSnapshot of the segregating-site derived-allele
// frequencies for each group at the given year. A single scan over genome positions
// tallies every group at once. Only sites polymorphic within a group are stored (a
// monomorphic site carries no temporal information); groups with fewer than two
// diploids are skipped. No RNG.
func snapshotFreqs(model *core.Model, pop *core.Pop, members map[int][]int, year int) *core.FreqSnapshot {
	snap := &core.FreqSnapshot{
		Year:     year,
		Groups:   map[int]map[int]float64{},
		HaploidN: map[int]int{},
	}
	// Sorted group keys for deterministic construction.
	groups := make([]int, 0, len(members))
	for g, m := range members {
		if len(m) >= 2 {
			groups = append(groups, g)
			snap.Groups[g] = map[int]float64{}
			snap.HaploidN[g] = 2 * len(m)
		}
	}
	sort.Ints(groups)

	contigs := buildContigs(model)
	for _, ct := range contigs {
		for _, arm := range ct.arms {
			start, length := arm[0], arm[1]
			for bit := start; bit < start+length; bit++ {
				for _, g := range groups {
					d := derivedCount(pop, members[g], bit)
					hn := snap.HaploidN[g]
					if d == 0 || d == hn {
						continue // monomorphic within this group
					}
					snap.Groups[g][bit] = float64(d) / float64(hn)
				}
			}
		}
	}
	return snap
}

// temporalNe estimates Ne for one group over the interval from prev to the current
// population using the Waples (1989) plan-II temporal method. It tracks the sites
// that were segregating in prev (recomputing each one's current frequency, so a
// now-fixed site correctly contributes 0 or 1) and forms the ratio-of-sums Fc. It
// returns +Inf when the sampling-corrected Fc is non-positive (drift below sampling
// noise ⇒ unresolved). Sites are visited in sorted order for reproducible summation.
func temporalNe(pop *core.Pop, prev *core.FreqSnapshot, group int, curMembers []int, tGen float64) (float64, int) {
	prevFreqs := prev.Groups[group]
	n0 := prev.HaploidN[group]
	nt := 2 * len(curMembers)
	if len(prevFreqs) == 0 || n0 < 2 || nt < 2 || tGen <= 0 {
		return math.Inf(1), 0
	}

	sites := make([]int, 0, len(prevFreqs))
	for s := range prevFreqs {
		sites = append(sites, s)
	}
	sort.Ints(sites)

	var num, den float64
	used := 0
	for _, s := range sites {
		x := prevFreqs[s]
		y := float64(derivedCount(pop, curMembers, s)) / float64(nt)
		z := (x+y)/2 - x*y // expected heterozygosity-like denominator
		if z <= 0 {
			continue // both samples fixed for the same allele: no information
		}
		diff := x - y
		num += diff * diff
		den += z
		used++
	}
	if den <= 0 || used == 0 {
		return math.Inf(1), used
	}
	fc := num / den
	corrected := fc - 1.0/float64(n0) - 1.0/float64(nt)
	if corrected <= 0 {
		return math.Inf(1), used // drift below sampling noise
	}
	return tGen / (2 * corrected), used
}

// derivedCount returns the number of derived (1) alleles across both strand copies
// of every member at genome bit s.
func derivedCount(pop *core.Pop, members []int, s int) int {
	d := 0
	for _, id := range members {
		c := pop.Chromosomes[id]
		if bitAt(c[0], s) {
			d++
		}
		if bitAt(c[1], s) {
			d++
		}
	}
	return d
}
