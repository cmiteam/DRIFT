package analysis

// Genetic load, DFE, and mutational-meltdown tracking — roadmap §6f, goal 1
// (creationist / genetic-entropy thesis).
//
// This is a READ-ONLY characterization layer over the de-novo mutation pool
// (pop.IndMutations + pop.MutationPool, gated by track_mutations — the same
// equilibrium object §6h/§6d/§6e read; the Chromosomes bitfield holds only
// decaying founder standing variation). It consumes NO RNG: every capture is a
// deterministic full-population scan, so a track_load run is byte-identical to a
// track_load-off run and strictly-neutral runs still pass `drift -validate`.
//
// WHAT IT MEASURES, and WHY these statistics. An empirical probe of the shipped
// engine (documented in TODO.md §6f) established three facts this layer is built
// to expose:
//
//   1. DRIFT has a DRIFT BARRIER at |s| ~ 1/(2*Ne). With the emergent Ne ~ 0.19*N
//      (§6h), mutations weaker than ~1/(2*0.19*N) are effectively neutral:
//      selection_mode has zero effect on them and they accumulate as drift alone
//      dictates. The ENTIRE shipped Weibull DFE (|s| ~ 5e-8) and even "moderate"
//      |s| ~ 1e-3 sit below this barrier — the nearly-neutral regime the genetic-
//      entropy argument lives in. The summary reports the barrier and the fraction
//      of the realized deleterious DFE that falls below it.
//
//   2. DRIFT does NOT spontaneously melt down to extinction at realistic params,
//      because its birth-then-cull-to-K demography regenerates the census each
//      year regardless of ABSOLUTE mean fitness (selection is soft/relative). Mean
//      fitness declines but census N does not, so there is no small-N -> more-drift
//      -> more-load feedback. The load time series therefore reaches a fluctuating
//      quasi-balance rather than a runaway; the summary states the verdict
//      (rising / stationary / declining) rather than assuming meltdown.
//
//   3. The decisive object is FIXED vs SEGREGATING load. Deleterious variants
//      below the barrier mostly stay segregating and rare; the ratchet only bites
//      (fixed deleterious accumulate) when Ne is small. NumFixedDel is that
//      counter; the DFE file contrasts the realized DFE of segregating vs fixed
//      variants — the depletion of large-|s| deleterious from the fixed column is
//      selection's fingerprint.
//
//   4. POOL GARBAGE-COLLECTION UNDER-COUNTS LOAD (found while building this layer).
//      DRIFT reference-counts pop.MutationPool: death.go decrements Mutation.Count
//      per dead copy and deletes the lineage at Count<=0, but InheritMutations
//      increments Count only on strand-1 inheritance (strand-0 bumps the global
//      pop.MutationCount instead). So Count systematically under-counts, and the
//      lineages with the most inheritance/death events — the COMMON and FIXED ones —
//      are garbage-collected FIRST, while still carried. Two consequences: (a) the
//      engine's own realized fitness under-counts load, because
//      CountFitnessAndMutations skips ids missing from the pool — a GC'd deleterious
//      mutation simply stops being expressed; and (b) FIXED load / the Muller's
//      ratchet are structurally invisible in the pool (in the §6f smoke run
//      ~98% of carried mutation ids were already GC'd, and NumFixed among pool-
//      resident lineages was 0). This layer therefore reports pool-RESIDENT load
//      only, exactly what selection acts on — the numbers are self-consistent with
//      the engine, and the GC is itself a (further) reason DRIFT does not melt down:
//      accumulated load is silently discarded. The opt-in
//      mutation_count_model="refcount" (utils.InheritMutations) corrects the Count
//      accounting so common/fixed lineages persist: with it on, THIS SAME layer
//      reports the true accumulating load and the fixed-deleterious ratchet
//      (NumFixedDel) that legacy hides (in the §6f smoke run mean fitness falls
//      1.0 -> 0.69 with 10 fixed deleterious by year 2000 under refcount, vs a
//      spurious flat 0.997 under legacy). Default "legacy" stays byte-identical;
//      refcount is opt-in (it changes the realized pool).
//
// Composition. The trajectory IS mean(CombineFitness), so it reflects
// fitness_model (additive / multiplicative / synergistic), epistasis_coefficient,
// and per-locus dominance exactly. selection_mode is purely observational (named
// in the summary header so the below-barrier mode-independence is visible).
// Everything partitions by mutation Class (stamped on each Mutation) and by Deme
// (§6b/§6c). The load layer pairs with the opt-in gamma DFE (dfe_model=gamma) so a
// literature-anchored input DFE can be fed and its realized fate measured here.

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"sort"
)

// emergentNeOverN is DRIFT's characterized neutral Ne/N ratio (§6h): the coalescent
// effective size is an emergent ~0.19 of the census size because the engine is
// overlapping-generations, monogamous, and high-reproductive-variance (NOT Wright-
// Fisher). Used only to report the operative selection (drift) barrier 1/(2*Ne);
// it is a documented reference value, not a live measurement.
const emergentNeOverN = 0.19

// CaptureLoad snapshots the living population's genetic load for the current year
// into pop.LoadHistory (roadmap §6f). Called from the run loop at the load_interval
// cadence when track_load is set. A deterministic full-population scan — consumes
// no RNG. Lazily initializes LoadHistory (nil until the first capture, like
// NeHistory/SFSHistory).
func CaptureLoad(model *core.Model, pop *core.Pop) {
	year := model.FreeParameters["year"]
	snap := computeLoadSnapshot(model, pop, year)
	if pop.LoadHistory == nil {
		pop.LoadHistory = make(map[int]*core.LoadSnapshot)
	}
	pop.LoadHistory[year] = snap
}

// computeLoadSnapshot scans the living population and the mutation pool once and
// returns the load statistics for the given year. It is the shared core of both
// CaptureLoad (time-series windows) and the end-of-run summary.
func computeLoadSnapshot(model *core.Model, pop *core.Pop, year int) *core.LoadSnapshot {
	living := sortedLivingIDs(pop)
	N := len(living)
	snap := &core.LoadSnapshot{Year: year, CensusN: N}
	if N == 0 {
		return snap
	}
	twoN := 2 * N

	copies := make(map[int]int)
	demeLoadSum := make(map[int]float64)
	demeCount := make(map[int]int)
	classLoadSum := make(map[int]float64)
	var fitSum, copySum, delCopySum float64

	for _, id := range living {
		muts := pop.IndMutations[id]
		for _, s := range []int{0, 1} {
			for _, m := range muts[s] {
				// Count only pool-RESIDENT mutations. DRIFT reference-counts the pool
				// (death.go decrements Mutation.Count per dead copy and deletes at 0),
				// but InheritMutations only increments Count on strand-1 inheritance
				// (strand-0 bumps the global MutationCount instead), so Count under-
				// counts and some still-carried mutations are garbage-collected early.
				// CountFitnessAndMutations already skips such orphaned ids, so selection
				// never sees them; counting them here (they'd default to Effect 0) would
				// desync the load/DFE view from the fitness the engine actually applies.
				if _, ok := pop.MutationPool[m]; !ok {
					continue
				}
				copies[m]++
				copySum++
				if pop.MutationPool[m].Effect < 0 {
					delCopySum++
				}
			}
		}
		nCopies, fit := utils.CountFitnessAndMutations(pop, id, model)
		// CountFitnessAndMutations returns 0.0 (not 1.0) for an individual with no
		// mutation copies (its len(strands)==0 early return) — e.g. founders before
		// any de-novo input. For LOAD purposes a mutation-free individual carries no
		// load, so its relative fitness is 1.0. (A child that inherited an empty
		// strand map already returns 1.0 via CombineFitness([]), so this only rescues
		// the never-mutated founders and keeps the year-0 window from reading load=1.)
		if nCopies == 0 {
			fit = 1.0
		}
		fitSum += fit
		deme := pop.IndData[id][individual.Deme]
		demeLoadSum[deme] += 1 - fit
		demeCount[deme]++
		accumulateClassLoad(pop, id, classLoadSum)
	}

	snap.MeanFitness = fitSum / float64(N)
	snap.TotalLoad = 1 - snap.MeanFitness
	snap.MeanMutPerInd = copySum / float64(N)
	snap.MeanDelPerInd = delCopySum / float64(N)

	// Fixed vs segregating classification and the fixed-only load. Sorted mutation
	// ids so the float sum in CombineFitness is reproducible under a fixed seed.
	ids := sortedKeys(copies)
	fixedContribs := make([]float64, 0)
	for _, m := range ids {
		c := copies[m]
		if c <= 0 || c > twoN {
			continue
		}
		eff := pop.MutationPool[m].Effect
		if c == twoN {
			snap.NumFixed++
			fixedContribs = append(fixedContribs, eff) // homozygous derived in all: full effect
			switch {
			case eff < 0:
				snap.NumFixedDel++
			case eff > 0:
				snap.NumFixedBen++
			}
		} else {
			snap.NumSeg++
		}
	}
	snap.FixedLoad = 1 - utils.CombineFitness(fixedContribs, model)
	snap.SegLoad = snap.TotalLoad - snap.FixedLoad

	if len(demeCount) > 1 {
		snap.DemeLoad = make(map[int]float64, len(demeCount))
		for d, n := range demeCount {
			snap.DemeLoad[d] = demeLoadSum[d] / float64(n)
		}
	}
	// Per-class additive load (per individual). Reported when a real class spectrum
	// is in play (>1 class, or a single non-point class index).
	if len(classLoadSum) > 1 || (len(classLoadSum) == 1 && !hasOnlyClassZero(classLoadSum)) {
		snap.ClassLoad = make(map[int]float64, len(classLoadSum))
		for cl, v := range classLoadSum {
			snap.ClassLoad[cl] = v / float64(N)
		}
	}
	return snap
}

// accumulateClassLoad adds individual ind's per-locus additive LOAD (= -fitness
// contribution; positive for deleterious) into classLoadSum keyed by mutation
// class. Mirrors the per-locus zygosity scoring in CountFitnessAndMutations (hom =
// full effect, het = h*effect) but groups the contributions by class. Sorted
// mutation ids keep the float sum reproducible. This is the additive decomposition
// (sums to the total additive load); under multiplicative/synergistic fitness it
// is a diagnostic that does not exactly equal TotalLoad.
func accumulateClassLoad(pop *core.Pop, ind int, classLoadSum map[int]float64) {
	strands := pop.IndMutations[ind]
	copies := make(map[int]int)
	for _, s := range []int{0, 1} {
		for _, m := range strands[s] {
			copies[m]++
		}
	}
	for _, m := range sortedKeys(copies) {
		mut, ok := pop.MutationPool[m]
		if !ok {
			continue
		}
		var contrib float64
		if copies[m] >= 2 {
			contrib = mut.Effect
		} else {
			contrib = float64(mut.Dominance) / 100.0 * mut.Effect
		}
		classLoadSum[mut.Class] += -contrib
	}
}

// ---------------------------------------------------------------------------
// DFE: input vs realized (segregating vs fixed), binned by signed log10 effect.
// ---------------------------------------------------------------------------

// dfeKey buckets a mutation by class, sign (deleterious/neutral/beneficial), and
// order of magnitude of |effect| (floor(log10|effect|); unused for neutral).
type dfeKey struct {
	class int
	sign  int // -1 deleterious, 0 neutral, +1 beneficial
	exp   int // floor(log10(|effect|)); 0 for neutral
}

func dfeKeyFor(mut core.Mutation) dfeKey {
	if mut.Effect == 0 {
		return dfeKey{mut.Class, 0, 0}
	}
	a := mut.Effect
	sign := 1
	if a < 0 {
		sign = -1
		a = -a
	}
	return dfeKey{mut.Class, sign, int(math.Floor(math.Log10(a)))}
}

// computeDFE tallies the realized DFE of the circulating mutation pool among the
// living, split into segregating (0<f<1) and fixed (f==1) lineages, per (class,
// sign, magnitude) bucket. Counts only pool-resident mutations (see the GC note in
// computeLoadSnapshot). The segregating-vs-fixed contrast is selection's
// fingerprint: under effective selection large-|s| deleterious lineages are
// depleted from the fixed column relative to the segregating one. (The historical
// INPUT DFE is not reported: DRIFT garbage-collects lost lineages from the pool, so
// it cannot be reconstructed reliably here.)
func computeDFE(pop *core.Pop) (seg, fixed map[dfeKey]int, keys []dfeKey) {
	seg = map[dfeKey]int{}
	fixed = map[dfeKey]int{}

	living := sortedLivingIDs(pop)
	twoN := 2 * len(living)
	copies := livingPoolCopies(pop, living)
	for _, m := range sortedKeys(copies) {
		c := copies[m]
		if c <= 0 || c > twoN {
			continue
		}
		k := dfeKeyFor(pop.MutationPool[m])
		if c == twoN {
			fixed[k]++
		} else {
			seg[k]++
		}
	}

	seen := map[dfeKey]bool{}
	for k := range seg {
		seen[k] = true
	}
	for k := range fixed {
		seen[k] = true
	}
	for k := range seen {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		a, b := keys[i], keys[j]
		if a.class != b.class {
			return a.class < b.class
		}
		if a.sign != b.sign {
			return a.sign < b.sign
		}
		return a.exp < b.exp
	})
	return seg, fixed, keys
}

// ---------------------------------------------------------------------------
// Output.
// ---------------------------------------------------------------------------

// SaveLoadTimeSeries writes the three §6f artifacts: the per-window load time
// series, the realized segregating-vs-fixed DFE, and a characterization summary
// (drift barrier + meltdown verdict). Wired into drift.go behind track_load; requires
// track_mutations (the de-novo pool). Ensures the final year is represented in the
// series even if the run ended off the capture cadence.
func SaveLoadTimeSeries(model *core.Model, pop *core.Pop) error {
	if model.Parameters["track_mutations"] != 1 {
		fmt.Println("Genetic-load tracking: track_mutations is off; nothing to measure (needs the de-novo mutation pool)")
		return nil
	}
	if pop.LoadHistory == nil {
		pop.LoadHistory = make(map[int]*core.LoadSnapshot)
	}
	// Guarantee an end-of-run window.
	finalYear := model.FreeParameters["year"]
	if _, ok := pop.LoadHistory[finalYear]; !ok {
		pop.LoadHistory[finalYear] = computeLoadSnapshot(model, pop, finalYear)
	}

	run := model.FreeParameters["run"]
	if err := writeLoadSeries(model, pop, run); err != nil {
		return err
	}
	if err := writeDFE(model, pop, run); err != nil {
		return err
	}
	return writeLoadSummary(model, pop, run)
}

func writeLoadSeries(model *core.Model, pop *core.Pop, run int) error {
	path := fmt.Sprintf("%s/%s_load_timeseries_run%d.csv", model.ResultsDir, model.ModelName, run)
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create load time series: %v", err)
	}
	defer f.Close()
	w := csv.NewWriter(f)
	defer w.Flush()

	if err := w.Write([]string{
		"Year", "CensusN", "MeanFitness", "TotalLoad", "SegLoad", "FixedLoad",
		"NumSeg", "NumFixed", "NumFixedDel", "NumFixedBen", "MeanMutPerInd", "MeanDelPerInd",
	}); err != nil {
		return err
	}
	for _, y := range sortedSnapshotYears(pop.LoadHistory) {
		s := pop.LoadHistory[y]
		if err := w.Write([]string{
			fmt.Sprintf("%d", s.Year), fmt.Sprintf("%d", s.CensusN),
			fmt.Sprintf("%.8f", s.MeanFitness), fmt.Sprintf("%.6g", s.TotalLoad),
			fmt.Sprintf("%.6g", s.SegLoad), fmt.Sprintf("%.6g", s.FixedLoad),
			fmt.Sprintf("%d", s.NumSeg), fmt.Sprintf("%d", s.NumFixed),
			fmt.Sprintf("%d", s.NumFixedDel), fmt.Sprintf("%d", s.NumFixedBen),
			fmt.Sprintf("%.3f", s.MeanMutPerInd), fmt.Sprintf("%.3f", s.MeanDelPerInd),
		}); err != nil {
			return err
		}
	}
	return nil
}

func writeDFE(model *core.Model, pop *core.Pop, run int) error {
	path := fmt.Sprintf("%s/%s_dfe_run%d.csv", model.ResultsDir, model.ModelName, run)
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create DFE file: %v", err)
	}
	defer f.Close()
	w := csv.NewWriter(f)
	defer w.Flush()

	seg, fixed, keys := computeDFE(pop)
	if err := w.Write([]string{
		"Class", "Sign", "EffectLo", "EffectHi", "SegCount", "FixedCount",
	}); err != nil {
		return err
	}
	signName := map[int]string{-1: "deleterious", 0: "neutral", 1: "beneficial"}
	for _, k := range keys {
		lo, hi := "0", "0"
		if k.sign != 0 {
			lo = fmt.Sprintf("%g", math.Pow(10, float64(k.exp)))
			hi = fmt.Sprintf("%g", math.Pow(10, float64(k.exp+1)))
		}
		if err := w.Write([]string{
			fmt.Sprintf("%d", k.class), signName[k.sign], lo, hi,
			fmt.Sprintf("%d", seg[k]), fmt.Sprintf("%d", fixed[k]),
		}); err != nil {
			return err
		}
	}
	return nil
}

func writeLoadSummary(model *core.Model, pop *core.Pop, run int) error {
	path := fmt.Sprintf("%s/%s_load_summary_run%d.csv", model.ResultsDir, model.ModelName, run)
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create load summary: %v", err)
	}
	defer f.Close()
	w := csv.NewWriter(f)
	defer w.Flush()

	final := computeLoadSnapshot(model, pop, model.FreeParameters["year"])
	N := final.CensusN
	// Operative selection (drift) barrier at the emergent Ne ~ 0.19*N (§6h).
	barrier := 0.0
	if N > 0 {
		ne := emergentNeOverN * float64(N)
		if ne > 0 {
			barrier = 1.0 / (2.0 * ne)
		}
	}
	belowBarrier, delTotal := deleteriousBelowBarrier(pop, barrier)
	belowFrac := 0.0
	if delTotal > 0 {
		belowFrac = float64(belowBarrier) / float64(delTotal)
	}
	verdict := meltdownVerdict(pop.LoadHistory)

	rows := [][]string{
		{"Metric", "Value"},
		{"selection_mode", model.StringParam("selection_mode", "fecundity")},
		{"fitness_model", model.StringParam("fitness_model", "additive")},
		{"dfe_model", model.StringParam("dfe_model", "weibull")},
		{"mu", fmt.Sprintf("%.6g", model.Parameters["mu"])},
		{"f_neutral", fmt.Sprintf("%.6g", model.Parameters["f_neutral"])},
		{"CensusN", fmt.Sprintf("%d", N)},
		{"MeanFitness", fmt.Sprintf("%.8f", final.MeanFitness)},
		{"TotalLoad", fmt.Sprintf("%.6g", final.TotalLoad)},
		{"SegLoad", fmt.Sprintf("%.6g", final.SegLoad)},
		{"FixedLoad", fmt.Sprintf("%.6g", final.FixedLoad)},
		{"NumSegregating", fmt.Sprintf("%d", final.NumSeg)},
		{"NumFixed", fmt.Sprintf("%d", final.NumFixed)},
		{"NumFixedDeleterious", fmt.Sprintf("%d", final.NumFixedDel)},
		{"NumFixedBeneficial", fmt.Sprintf("%d", final.NumFixedBen)},
		{"MeanMutPerInd", fmt.Sprintf("%.3f", final.MeanMutPerInd)},
		{"MeanDelPerInd", fmt.Sprintf("%.3f", final.MeanDelPerInd)},
		{"EmergentNe", fmt.Sprintf("%.1f", emergentNeOverN*float64(N))},
		{"DriftBarrier_s", fmt.Sprintf("%.6g", barrier)},
		{"DelBelowBarrierFrac", fmt.Sprintf("%.4f", belowFrac)},
		{"MeltdownVerdict", verdict},
	}
	for _, r := range rows {
		if err := w.Write(r); err != nil {
			return err
		}
	}
	return nil
}

// deleteriousBelowBarrier counts, among the living population's segregating+fixed
// deleterious lineages, how many have |effect| below the drift barrier (effectively
// neutral: selection cannot see them). Returns (belowBarrier, totalDeleterious).
func deleteriousBelowBarrier(pop *core.Pop, barrier float64) (below, total int) {
	living := sortedLivingIDs(pop)
	twoN := 2 * len(living)
	copies := livingPoolCopies(pop, living)
	for _, m := range sortedKeys(copies) {
		c := copies[m]
		if c <= 0 || c > twoN {
			continue
		}
		eff := pop.MutationPool[m].Effect
		if eff >= 0 {
			continue
		}
		total++
		if -eff < barrier {
			below++
		}
	}
	return below, total
}

// meltdownVerdict compares mean total load in the first vs the last third of the
// captured windows: rising (accumulating toward meltdown), declining (selection
// purging), or stationary (mutation-selection-drift quasi-balance). Needs >= 3
// windows; otherwise "insufficient-windows".
func meltdownVerdict(hist map[int]*core.LoadSnapshot) string {
	years := sortedSnapshotYears(hist)
	if len(years) < 3 {
		return "insufficient-windows"
	}
	third := len(years) / 3
	if third < 1 {
		third = 1
	}
	var firstSum, lastSum float64
	for i := 0; i < third; i++ {
		firstSum += hist[years[i]].TotalLoad
	}
	for i := len(years) - third; i < len(years); i++ {
		lastSum += hist[years[i]].TotalLoad
	}
	first := firstSum / float64(third)
	last := lastSum / float64(third)
	switch {
	case last > first*1.05:
		return "rising (load accumulating)"
	case last < first*0.95:
		return "declining (selection purging)"
	default:
		return "stationary (quasi-balance)"
	}
}

// ---------------------------------------------------------------------------
// Small deterministic helpers.
// ---------------------------------------------------------------------------

func sortedLivingIDs(pop *core.Pop) []int {
	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	sort.Ints(ids)
	return ids
}

func sortedKeys(m map[int]int) []int {
	ks := make([]int, 0, len(m))
	for k := range m {
		ks = append(ks, k)
	}
	sort.Ints(ks)
	return ks
}

// livingPoolCopies returns the number of carried copies of each pool-RESIDENT
// mutation id across the living population. Orphaned ids (garbage-collected from
// MutationPool while still carried; see computeLoadSnapshot) are skipped so the
// view matches the fitness the engine actually applies.
func livingPoolCopies(pop *core.Pop, living []int) map[int]int {
	copies := make(map[int]int)
	for _, id := range living {
		muts := pop.IndMutations[id]
		for _, s := range []int{0, 1} {
			for _, m := range muts[s] {
				if _, ok := pop.MutationPool[m]; ok {
					copies[m]++
				}
			}
		}
	}
	return copies
}

func sortedSnapshotYears(hist map[int]*core.LoadSnapshot) []int {
	ys := make([]int, 0, len(hist))
	for y := range hist {
		ys = append(ys, y)
	}
	sort.Ints(ys)
	return ys
}

func hasOnlyClassZero(m map[int]float64) bool {
	for k := range m {
		if k != 0 {
			return false
		}
	}
	return true
}
