package simulation

import (
	"drift/pkg/analysis"
	"drift/pkg/utils"
	"sort"
	"testing"
)

// W3 founder-allele labels driven through the real birth/death path (TMR4A.md W3).
// The claim under test is the one the study rests on: with fewer created alleles than
// founder strands, the number of distinct CREATED alleles surviving at a locus is
// strictly below the number of distinct founder haplotypes — identity by state and
// identity by descent come apart, and the walker reports both rather than conflating them.

// labelledRun runs a fixed mini-simulation, optionally labelling founders with A created
// alleles, and returns the walker's loci plus the final RNG state.
func labelledRun(t *testing.T, alleles int) ([]analysis.TMRKALocus, []byte) {
	t.Helper()
	utils.SeedRNG(31415)
	const founderPairs = 8 // 16 founders ⇒ 32 founder strands
	model := pedModel(1)
	model.Parameters["track_DNA"] = 1
	model.Parameters["track_arg"] = 1
	model.Parameters["tmrka_k"] = 4
	model.StringParams["arg_loci"] = "0-255:16"
	model.ChromosomeArms = map[int]map[int][]int{1: {0: {0, 128}, 1: {128, 128}}}
	model.FreeParameters["genome_bits"] = 256

	pop := pedPop(model, founderPairs)
	if alleles > 0 {
		ids := make([]int, 0, len(pop.IndData))
		for id := range pop.IndData {
			ids = append(ids, id)
		}
		pop.AssignFounderLabels(ids, alleles)
	}
	runPedSim(t, model, pop, 26, 4)

	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	sort.Ints(ids)

	res := analysis.ComputeTMRKA(model, pop, ids)
	if res == nil {
		t.Fatal("walk produced nothing")
	}
	return res.Loci, utils.RNGState()
}

// A = 4 over 32 founder strands: eight founder haplotypes per created allele, so wherever
// more than one founder haplotype survives they usually share a created allele.
func TestTMRKACreatedAllelesUnderShipsFounderStrands(t *testing.T) {
	loci, _ := labelledRun(t, 4)
	strictlyFewer := 0
	for _, l := range loci {
		if l.CreatedAllelesSurviving < 1 {
			t.Fatalf("locus %d: %d created alleles", l.Position, l.CreatedAllelesSurviving)
		}
		if l.CreatedAllelesSurviving > l.FounderStrandsSurviving {
			t.Errorf("locus %d: %d created alleles from %d founder haplotypes — impossible",
				l.Position, l.CreatedAllelesSurviving, l.FounderStrandsSurviving)
		}
		if l.CreatedAllelesSurviving > 4 {
			t.Errorf("locus %d: %d created alleles but only 4 were ever created",
				l.Position, l.CreatedAllelesSurviving)
		}
		if l.UnlabelledSources != 0 {
			t.Errorf("locus %d: %d unlabelled sources in a fully labelled run",
				l.Position, l.UnlabelledSources)
		}
		if l.CreatedAllelesSurviving < l.FounderStrandsSurviving {
			strictlyFewer++
		}
	}
	if strictlyFewer == 0 {
		t.Error("created-allele count never fell below the founder-haplotype count; " +
			"identity by state and identity by descent never came apart, so the test proved nothing")
	}
	t.Logf("%d/%d loci have strictly fewer created alleles than founder haplotypes",
		strictlyFewer, len(loci))
}

// A = 1 is the null arm: every founder strand carries the same created allele, so NO pair
// is separated by created divergence and every difference in the run is mutational.
func TestTMRKASingleCreatedAlleleLeavesNoCreatedPairs(t *testing.T) {
	loci, _ := labelledRun(t, 1)
	for _, l := range loci {
		if l.CreatedAllelesSurviving != 1 {
			t.Errorf("locus %d: %d created alleles, want 1", l.Position, l.CreatedAllelesSurviving)
		}
		if l.CreatedPairs != 0 {
			t.Errorf("locus %d: %d created pairs under a single created allele",
				l.Position, l.CreatedPairs)
		}
	}
}

// The §2 gate: labelling draws no RNG and changes no by-descent quantity. Only the
// created columns move.
func TestTMRKALabelsAreByteIdenticalAndInert(t *testing.T) {
	off, offState := labelledRun(t, 0)
	on, onState := labelledRun(t, 4)

	if string(offState) != string(onState) {
		t.Error("founder-allele labelling perturbed the RNG stream")
	}
	if len(off) != len(on) {
		t.Fatalf("locus counts differ: %d vs %d", len(off), len(on))
	}
	for i := range off {
		if off[i].Position != on[i].Position ||
			off[i].Outcome != on[i].Outcome ||
			off[i].TMRKAYear != on[i].TMRKAYear ||
			off[i].FinalLineages != on[i].FinalLineages ||
			off[i].FounderStrandsSurviving != on[i].FounderStrandsSurviving ||
			off[i].CoalescedPairs != on[i].CoalescedPairs {
			t.Fatalf("locus %d: labelling changed a by-descent quantity\n  off: %+v\n  on:  %+v",
				off[i].Position, off[i], on[i])
		}
		if off[i].CreatedAllelesSurviving != -1 || off[i].CreatedPairs != -1 {
			t.Errorf("locus %d: created columns reported without labels", off[i].Position)
		}
		if on[i].CreatedAllelesSurviving < 1 || on[i].CreatedPairs < 0 {
			t.Errorf("locus %d: created columns missing with labels on", on[i].Position)
		}
	}
}
