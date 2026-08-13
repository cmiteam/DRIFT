package simulation

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"testing"
)

// Integration tests for diploid pedigree retention (TMR4A.md W1) driven through the
// REAL birthStandard / RIP paths, not a hand-rolled replay — the point is that the two
// hooks are placed correctly and that the flag is genuinely inert when off.

// pedModel builds the smallest model birthStandard will run: everyone mature, nobody in
// menopause, no birth spacing, birth_prob 1 (RandIntn(1) is always 0), no wandering, and
// no DNA/mutation tracking, so the only RNG draws per birth are the pregnancy roll, the
// fitness roll, and the child's sex.
func pedModel(trackPedigree float64) *core.Model {
	return &core.Model{
		Parameters: map[string]float64{
			"maturity":           0,
			"menopause":          1,
			"spacing":            0,
			"birth_prob":         1,
			"mu_scale_factor":    1,
			"lifespan_drop":      1,
			"min_lifespan":       100,
			"lifespan":           100,
			"wander":             0,
			"track_DNA":          0,
			"track_mutations":    0,
			"track_coalescence":  0,
			"track_pedigree":     trackPedigree,
			"recomb_obligate_cM": 50,
		},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"year": 0, "indID": 0, "genome_bits": 64},
		ChromosomeArms: map[int]map[int][]int{},
	}
}

// pedPop seeds nPairs married couples (even ids male, odd ids female) and, when the flag
// is on, gives each founder a pedigree node exactly as config.InitializePop does.
func pedPop(model *core.Model, nPairs int) *core.Pop {
	pop := &core.Pop{
		IndData:      map[int][]int{},
		Chromosomes:  map[int][][]uint64{},
		Centromeres:  map[int][2][]uint64{},
		IndMutations: map[int]map[int][]int{},
		MutationPool: map[int]core.Mutation{},
		MutationHist: map[int]int{},
		Tracking:     map[string]int{},
		MaleDB:       map[int]core.Ancestor{},
		FemaleDB:     map[int]core.Ancestor{},
	}
	for i := 0; i < 2*nPairs; i++ {
		d := individual.MakeIndData()
		d[individual.Sex] = i % 2 // even = male, odd = female
		d[individual.BirthYear] = -20
		d[individual.Lifespan] = 100
		d[individual.Fitness] = 1
		d[individual.NumBirths] = 0
		d[individual.LastBirthYear] = -100
		d[individual.MarriageState] = i ^ 1 // pair 0-1, 2-3, ...
		pop.IndData[i] = d
		if model.Parameters["track_pedigree"] == 1 {
			pop.PedEnsure(i, d[individual.BirthYear], d[individual.Sex])
		}
	}
	model.FreeParameters["indID"] = 2 * nPairs
	return pop
}

// runPedSim advances `years` years of births, then kills every individual whose id is
// divisible by killMod, exercising the prune cascade against real RIP calls.
func runPedSim(t *testing.T, model *core.Model, pop *core.Pop, years, killMod int) {
	t.Helper()
	for y := 1; y <= years; y++ {
		model.FreeParameters["year"] = y
		birthStandard(model, pop)
		if y%3 == 0 {
			// Deterministic culling: sorted so the kill order is reproducible.
			var victims []int
			for id := range pop.IndData {
				if id%killMod == 0 {
					victims = append(victims, id)
				}
			}
			sortInts(victims)
			for _, v := range victims {
				RIP(v, pop, model)
			}
		}
	}
}

func sortInts(a []int) {
	for i := 1; i < len(a); i++ {
		for j := i; j > 0 && a[j] < a[j-1]; j-- {
			a[j], a[j-1] = a[j-1], a[j]
		}
	}
}

// THE GATE (TMR4A.md §2): every item must leave a default run byte-identical. The
// pedigree path draws no RNG and writes no existing field, so an identical seed must
// leave the RNG in the identical state and produce the identical population.
func TestPedigreeIsByteIdenticalWhenOff(t *testing.T) {
	run := func(track float64) ([]byte, map[int][]int, int) {
		utils.SeedRNG(20260813)
		model := pedModel(track)
		pop := pedPop(model, 6)
		runPedSim(t, model, pop, 12, 4)
		return utils.RNGState(), pop.IndData, model.FreeParameters["indID"]
	}
	offState, offInds, offNext := run(0)
	onState, onInds, onNext := run(1)

	if string(offState) != string(onState) {
		t.Error("track_pedigree consumed RNG: the flag must not perturb the stream")
	}
	if offNext != onNext {
		t.Errorf("different number of births: off=%d on=%d", offNext, onNext)
	}
	if len(offInds) != len(onInds) {
		t.Fatalf("different survivors: off=%d on=%d", len(offInds), len(onInds))
	}
	for id, off := range offInds {
		on, ok := onInds[id]
		if !ok {
			t.Fatalf("individual %d present only with the flag off", id)
		}
		for f := range off {
			if off[f] != on[f] {
				t.Fatalf("individual %d field %d: off=%d on=%d", id, f, off[f], on[f])
			}
		}
	}
}

// With the flag off nothing is recorded at all — no silent allocation on the hot path.
func TestPedigreeNotRecordedWhenOff(t *testing.T) {
	utils.SeedRNG(7)
	model := pedModel(0)
	pop := pedPop(model, 4)
	runPedSim(t, model, pop, 6, 5)
	if pop.Pedigree != nil {
		t.Errorf("pedigree map allocated with track_pedigree off (%d nodes)", len(pop.Pedigree))
	}
}

// Both parents of every living individual must be reachable — including a son's mother
// and a daughter's father, which MaleDB/FemaleDB do not record. This is the actual gap
// W1 closes, so assert it against the existing chains rather than in the abstract.
func TestPedigreeRecordsWhatCoalescenceChainsCannot(t *testing.T) {
	utils.SeedRNG(99)
	model := pedModel(1)
	pop := pedPop(model, 6)
	runPedSim(t, model, pop, 15, 4)

	checked := 0
	for id, d := range pop.IndData {
		n := pop.Pedigree[id]
		if n == nil {
			t.Fatalf("living individual %d has no pedigree node", id)
		}
		if n.Dad == core.PedFounder && n.Mom == core.PedFounder {
			continue // a founder: no parents by construction
		}
		checked++
		if n.Dad != d[individual.Dad] || n.Mom != d[individual.Mom] {
			t.Errorf("individual %d: pedigree says %d/%d, IndData says %d/%d",
				id, n.Dad, n.Mom, d[individual.Dad], d[individual.Mom])
		}
		if n.BirthYear != d[individual.BirthYear] {
			t.Errorf("individual %d: pedigree birth year %d, IndData %d",
				id, n.BirthYear, d[individual.BirthYear])
		}
		// The cross-sex parent is exactly what the single-sex chains drop.
		crossSexParent := n.Mom
		if d[individual.Sex] == 1 {
			crossSexParent = n.Dad
		}
		if crossSexParent == core.PedFounder {
			t.Errorf("individual %d has no recorded cross-sex parent", id)
		}
	}
	if checked == 0 {
		t.Fatal("no non-founders were born; the test proved nothing")
	}
}

// Every ancestor of a living individual must survive its own death — the retention
// property the walker depends on — and nothing else may.
func TestPedigreeRetainsAncestorsOfTheLiving(t *testing.T) {
	utils.SeedRNG(4242)
	model := pedModel(1)
	pop := pedPop(model, 6)
	runPedSim(t, model, pop, 18, 3)

	// Closure of the living over both parent edges: everything the walker can reach.
	reachable := map[int]bool{}
	var stack []int
	for id := range pop.IndData {
		stack = append(stack, id)
	}
	for len(stack) > 0 {
		cur := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		if reachable[cur] || cur == core.PedFounder {
			continue
		}
		reachable[cur] = true
		n, ok := pop.Pedigree[cur]
		if !ok {
			t.Fatalf("ancestor %d of a living individual was pruned", cur)
		}
		stack = append(stack, n.Dad, n.Mom)
	}
	for id := range pop.Pedigree {
		if !reachable[id] {
			t.Errorf("node %d is retained but unreachable from any living individual", id)
		}
	}
	if len(reachable) <= len(pop.IndData) {
		t.Fatal("no dead ancestors were retained; the retention rule was not exercised")
	}
	t.Logf("living=%d retained pedigree nodes=%d", len(pop.IndData), len(pop.Pedigree))
}

// Extinction empties the pedigree completely: no leak survives a total wipe.
func TestPedigreeEmptiesOnExtinction(t *testing.T) {
	utils.SeedRNG(11)
	model := pedModel(1)
	pop := pedPop(model, 4)
	runPedSim(t, model, pop, 10, 5)

	var all []int
	for id := range pop.IndData {
		all = append(all, id)
	}
	sortInts(all)
	for _, id := range all {
		RIP(id, pop, model)
	}
	if len(pop.Pedigree) != 0 {
		t.Errorf("refcount leak: %d nodes retained after total extinction", len(pop.Pedigree))
	}
}
