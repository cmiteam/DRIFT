package simulation

import (
	"drift/pkg/analysis"
	"drift/pkg/utils"
	"sort"
	"testing"
)

// The locus-sweep escape hatch (TMR4A.md §5, required by W8).
//
// The plan for scaling past what fits in memory is to re-run with a DIFFERENT focal
// locus set under the same seed and rely on the genealogy being byte-identical, sweeping
// loci in passes instead of holding a whole ARG at once. §5 is explicit that "W8 must
// actually test this, or the whole strategy rests on an assumption" — so this drives the
// real birth/death path twice with different arg_loci and compares the walker's answers
// at the loci the two runs share.

// sweepRun runs a fixed mini-simulation with the given focal locus spec and returns the
// walker's per-locus results, keyed by genome-bit position, plus the final RNG state.
func sweepRun(t *testing.T, argLoci string) (map[int]analysis.TMRKALocus, []byte, int) {
	t.Helper()
	utils.SeedRNG(90210)
	model := pedModel(1)
	model.Parameters["track_DNA"] = 1
	model.Parameters["track_arg"] = 1
	model.Parameters["tmrka_k"] = 4
	model.Parameters["generation_time"] = 25
	model.StringParams["arg_loci"] = argLoci
	model.ChromosomeArms = map[int]map[int][]int{
		1: {0: {0, 64}, 1: {64, 64}},
		2: {0: {128, 64}, 1: {192, 64}},
	}
	model.FreeParameters["genome_bits"] = 256

	pop := pedPop(model, 8)
	runPedSim(t, model, pop, 24, 5)

	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	sort.Ints(ids)

	res := analysis.ComputeTMRKA(model, pop, ids)
	out := map[int]analysis.TMRKALocus{}
	if res != nil {
		for _, l := range res.Loci {
			out[l.Position] = l
		}
	}
	return out, utils.RNGState(), len(ids)
}

func TestTMRKALocusSweepIsReproducible(t *testing.T) {
	full, fullState, fullN := sweepRun(t, "0-255:16")
	if len(full) != 16 {
		t.Fatalf("expected 16 focal loci, got %d", len(full))
	}

	// A second pass over a DIFFERENT subset of loci, same seed. The genealogy is a
	// property of the run, not of which loci were watched, so every shared locus must
	// come back identical — trajectory included.
	part, partState, partN := sweepRun(t, "32,96,240")
	if len(part) != 3 {
		t.Fatalf("expected 3 focal loci, got %d", len(part))
	}
	if string(fullState) != string(partState) || fullN != partN {
		t.Fatal("the locus set perturbed the run itself; sweeping loci in passes is unsound")
	}

	for pos, got := range part {
		want, ok := full[pos]
		if !ok {
			t.Fatalf("locus %d missing from the full pass", pos)
		}
		if got.Outcome != want.Outcome || got.TMRKAYear != want.TMRKAYear ||
			got.HasTMRKA != want.HasTMRKA || got.Coalescences != want.Coalescences ||
			got.FinalLineages != want.FinalLineages {
			t.Errorf("locus %d differs between passes:\n  sweep: %+v\n  full:  %+v", pos, got, want)
		}
		if len(got.Trajectory) != len(want.Trajectory) {
			t.Fatalf("locus %d: trajectory lengths %d vs %d", pos, len(got.Trajectory), len(want.Trajectory))
		}
		for i := range got.Trajectory {
			if got.Trajectory[i] != want.Trajectory[i] {
				t.Fatalf("locus %d trajectory step %d: %+v vs %+v",
					pos, i, got.Trajectory[i], want.Trajectory[i])
			}
		}
	}
}

// The same run walked twice must give the same answer — the walk itself must not depend
// on Go's randomized map iteration order.
func TestTMRKAWalkIsDeterministic(t *testing.T) {
	utils.SeedRNG(1234)
	model := pedModel(1)
	model.Parameters["track_DNA"] = 1
	model.Parameters["track_arg"] = 1
	model.StringParams["arg_loci"] = "0-255:32"
	model.ChromosomeArms = map[int]map[int][]int{1: {0: {0, 128}, 1: {128, 128}}}
	model.FreeParameters["genome_bits"] = 256

	pop := pedPop(model, 8)
	runPedSim(t, model, pop, 20, 4)

	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	sort.Ints(ids)

	a := analysis.ComputeTMRKA(model, pop, ids)
	b := analysis.ComputeTMRKA(model, pop, ids)
	if a == nil || b == nil {
		t.Fatal("walk produced nothing")
	}
	if len(a.Loci) != len(b.Loci) {
		t.Fatalf("locus counts differ: %d vs %d", len(a.Loci), len(b.Loci))
	}
	for i := range a.Loci {
		if a.Loci[i].Position != b.Loci[i].Position ||
			a.Loci[i].Outcome != b.Loci[i].Outcome ||
			a.Loci[i].TMRKAYear != b.Loci[i].TMRKAYear ||
			a.Loci[i].FinalLineages != b.Loci[i].FinalLineages {
			t.Fatalf("locus %d: %+v vs %+v", a.Loci[i].Position, a.Loci[i], b.Loci[i])
		}
		for j := range a.Loci[i].Trajectory {
			if a.Loci[i].Trajectory[j] != b.Loci[i].Trajectory[j] {
				t.Fatalf("locus %d trajectory step %d differs between identical walks",
					a.Loci[i].Position, j)
			}
		}
	}
}

// Sanity anchor: over a run with a founding population of F diploids, the number of
// distinct founder haplotypes surviving at a locus can never exceed 2F, must be at
// least 1, and must fall over time as drift removes lineages. This is the headline
// falsifiable quantity (TMR4A.md §0), so its bounds are worth pinning.
func TestTMRKAFounderHaplotypeBounds(t *testing.T) {
	utils.SeedRNG(555)
	const founderPairs = 8
	model := pedModel(1)
	model.Parameters["track_DNA"] = 1
	model.Parameters["track_arg"] = 1
	model.StringParams["arg_loci"] = "0-255:16"
	model.ChromosomeArms = map[int]map[int][]int{1: {0: {0, 128}, 1: {128, 128}}}
	model.FreeParameters["genome_bits"] = 256

	pop := pedPop(model, founderPairs)
	runPedSim(t, model, pop, 30, 4)

	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	sort.Ints(ids)

	res := analysis.ComputeTMRKA(model, pop, ids)
	if res == nil || len(res.Loci) == 0 {
		t.Fatal("walk produced nothing")
	}
	maxFounderStrands := 2 * 2 * founderPairs // 2 strands × the founding individuals
	total := 0
	for _, l := range res.Loci {
		if l.FinalLineages < 1 {
			t.Errorf("locus %d: %d surviving founder haplotypes", l.Position, l.FinalLineages)
		}
		if l.FinalLineages > maxFounderStrands {
			t.Errorf("locus %d: %d surviving founder haplotypes exceeds the %d the founders could transmit",
				l.Position, l.FinalLineages, maxFounderStrands)
		}
		if l.FinalLineages > l.InitialLineages {
			t.Errorf("locus %d: lineages grew backwards in time", l.Position)
		}
		total += l.FinalLineages
	}
	mean := float64(total) / float64(len(res.Loci))
	if mean >= float64(maxFounderStrands) {
		t.Errorf("mean surviving founder haplotypes %.2f — drift removed none over 30 years", mean)
	}
	t.Logf("mean founder haplotypes surviving: %.2f of a possible %d (n=%d sampled, %d loci)",
		mean, maxFounderStrands, len(ids), len(res.Loci))
}
