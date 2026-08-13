package simulation

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"testing"
)

// Strand provenance (TMR4A.md W2) validated against the REAL meiosis it claims to
// describe. The whole module rests on one inversion — meiosis takes parental copy 0
// where the mask bit is SET — so the decisive test is not "does the bitset round-trip"
// but "does the recorded strand predict which parental allele actually arrived".

// argParent builds a parent whose two strands are trivially distinguishable: strand 0 is
// all zeros, strand 1 is all ones. After meiosis, the child's bit AT a position is then
// exactly the index of the strand it came from.
func argParent(words int) [][]uint64 {
	c0 := make([]uint64, words)
	c1 := make([]uint64, words)
	for i := range c1 {
		c1[i] = ^uint64(0)
	}
	return [][]uint64{c0, c1}
}

// THE LOAD-BEARING TEST. For every focal locus, the recorded strand must equal the bit
// meiosis actually deposited in the child. Run over many independent masks so both
// crossover configurations and both starting homologs are exercised.
func TestARGRecordedStrandMatchesMeiosisResult(t *testing.T) {
	model := linkageModel("arm") // one chromosome, p [0,100), q [100,200), genome_bits 200
	model.Parameters["track_arg"] = 1
	model.StringParams["arg_loci"] = "0,1,50,99,100,101,150,199"
	loci := core.EnsureARGLoci(model)
	if len(loci) != 8 {
		t.Fatalf("focal loci parsed as %v", loci)
	}
	words := (model.FreeParameters["genome_bits"] + 63) / 64

	utils.SeedRNG(31337)
	const trials = 400
	seenStrand := [2]int{}
	for tr := 0; tr < trials; tr++ {
		pop := &core.Pop{
			IndData: map[int][]int{1: individual.MakeIndData(), 2: individual.MakeIndData()},
			Chromosomes: map[int][][]uint64{
				1: argParent(words),
				2: {make([]uint64, words), make([]uint64, words)},
			},
		}
		maskDad, _ := createMask(model, 0)
		maskMom, _ := createMask(model, 1)

		pop.ARGCapture(2, loci, maskDad, maskMom)
		meiosis(pop, maskDad, 1, 2, 0) // child strand 0 from the "dad" gamete

		for i, pos := range loci {
			recorded, ok := pop.ARGStrand(2, 0, i, len(loci))
			if !ok {
				t.Fatalf("trial %d locus %d: no record", tr, pos)
			}
			actual := int((pop.Chromosomes[2][0][pos/64] >> (uint(pos) % 64)) & 1)
			if recorded != actual {
				t.Fatalf("trial %d locus %d: recorded strand %d but meiosis delivered strand %d",
					tr, pos, recorded, actual)
			}
			seenStrand[recorded]++
		}
	}
	// Guard against a vacuous pass: if every draw had come from one strand the equality
	// above would hold trivially.
	if seenStrand[0] == 0 || seenStrand[1] == 0 {
		t.Fatalf("only one strand was ever sampled (%v); the test proved nothing", seenStrand)
	}
}

// The same check under the §6a real-map recombination path, which builds masks by a
// completely different route (crossovers placed in cM space, 50/50 starting homolog).
func TestARGRecordedStrandMatchesMeiosisUnderRecombMap(t *testing.T) {
	model := linkageModel("arm")
	model.StringParams["recombination_model"] = "map"
	model.ArmCM = map[int]map[int]float64{1: {0: 60, 1: 60}}
	model.Parameters["recomb_obligate_cM"] = 50
	model.Parameters["track_arg"] = 1
	model.StringParams["arg_loci"] = "0-199:7"
	loci := core.EnsureARGLoci(model)
	words := (model.FreeParameters["genome_bits"] + 63) / 64

	utils.SeedRNG(5150)
	seenStrand := [2]int{}
	for tr := 0; tr < 200; tr++ {
		pop := &core.Pop{
			IndData: map[int][]int{1: individual.MakeIndData(), 2: individual.MakeIndData()},
			Chromosomes: map[int][][]uint64{
				1: argParent(words),
				2: {make([]uint64, words), make([]uint64, words)},
			},
		}
		mask, _ := createMask(model, 0)
		pop.ARGCapture(2, loci, mask, mask)
		meiosis(pop, mask, 1, 2, 0)
		for i, pos := range loci {
			recorded, _ := pop.ARGStrand(2, 0, i, len(loci))
			actual := int((pop.Chromosomes[2][0][pos/64] >> (uint(pos) % 64)) & 1)
			if recorded != actual {
				t.Fatalf("map mode, trial %d locus %d: recorded %d, delivered %d", tr, pos, recorded, actual)
			}
			seenStrand[recorded]++
		}
	}
	if seenStrand[0] == 0 || seenStrand[1] == 0 {
		t.Fatalf("only one strand was ever sampled (%v)", seenStrand)
	}
}

// Adjacent focal loci must usually agree (they are linked) but not always (recombination
// separates them). A recorded provenance that were constant along the chromosome, or
// independent between neighbours, would both pass the equality test above on a
// degenerate mask; this pins the actual linkage structure.
func TestARGProvenanceIsLinkedButBreakable(t *testing.T) {
	model := linkageModel("arm")
	model.Parameters["track_arg"] = 1
	model.StringParams["arg_loci"] = "0-99:1"
	loci := core.EnsureARGLoci(model)

	utils.SeedRNG(2718)
	switches, comparisons := 0, 0
	for tr := 0; tr < 300; tr++ {
		pop := &core.Pop{}
		mask, _ := createMask(model, 0)
		pop.ARGCapture(1, loci, mask, mask)
		prev, _ := pop.ARGStrand(1, 0, 0, len(loci))
		for i := 1; i < len(loci); i++ {
			cur, _ := pop.ARGStrand(1, 0, i, len(loci))
			if cur != prev {
				switches++
			}
			prev = cur
			comparisons++
		}
	}
	rate := float64(switches) / float64(comparisons)
	if switches == 0 {
		t.Error("provenance never switched strand along the chromosome: no recombination recorded")
	}
	if rate > 0.5 {
		t.Errorf("adjacent focal loci switch strand %.3f of the time; provenance is not linked", rate)
	}
	t.Logf("strand switches between adjacent bits: %d/%d (%.4f)", switches, comparisons, rate)
}

// Capture on the real birth path: every non-founder gets a provenance record, and it is
// pruned exactly when its pedigree node is.
func TestARGCaptureOnBirthPath(t *testing.T) {
	utils.SeedRNG(808)
	model := pedModel(1)
	model.Parameters["track_DNA"] = 0
	model.Parameters["track_mutations"] = 0
	model.Parameters["track_arg"] = 1
	model.StringParams["arg_loci"] = "0,17,33"
	model.ChromosomeArms = map[int]map[int][]int{1: {0: {0, 32}, 1: {32, 32}}}
	model.FreeParameters["genome_bits"] = 64

	pop := pedPop(model, 5)
	runPedSim(t, model, pop, 12, 4)

	nLoci := len(core.EnsureARGLoci(model))
	born := 0
	for id := range pop.IndData {
		n := pop.Pedigree[id]
		if n.Dad == core.PedFounder {
			continue // founder: no meiosis produced it, so no provenance by design
		}
		born++
		for g := 0; g < 2; g++ {
			for i := 0; i < nLoci; i++ {
				if _, ok := pop.ARGStrand(id, g, i, nLoci); !ok {
					t.Fatalf("individual %d gamete %d locus %d: no provenance recorded", id, g, i)
				}
			}
		}
	}
	if born == 0 {
		t.Fatal("no individuals were born; the test proved nothing")
	}
	for id := range pop.ARGDB {
		if _, ok := pop.Pedigree[id]; !ok {
			t.Errorf("provenance for %d survived its pedigree node", id)
		}
	}
}

// track_arg without track_pedigree records nothing: the records would be unreadable
// (no parent links) and unbounded (nothing prunes them).
func TestARGRequiresPedigree(t *testing.T) {
	utils.SeedRNG(808)
	model := pedModel(0) // track_pedigree off
	model.Parameters["track_arg"] = 1
	model.StringParams["arg_loci"] = "0,17"
	model.ChromosomeArms = map[int]map[int][]int{1: {0: {0, 32}, 1: {32, 32}}}
	pop := pedPop(model, 4)
	runPedSim(t, model, pop, 6, 5)
	if len(pop.ARGDB) != 0 {
		t.Errorf("recorded %d provenance entries with track_pedigree off", len(pop.ARGDB))
	}
}

// The §2 gate: with the flags off, nothing is captured and the RNG stream is untouched.
func TestARGOffIsByteIdentical(t *testing.T) {
	run := func(argOn bool) ([]byte, int) {
		utils.SeedRNG(606)
		model := pedModel(1)
		model.Parameters["track_DNA"] = 1
		model.ChromosomeArms = map[int]map[int][]int{1: {0: {0, 32}, 1: {32, 32}}}
		model.FreeParameters["genome_bits"] = 64
		if argOn {
			model.Parameters["track_arg"] = 1
			model.StringParams["arg_loci"] = "0,17,33"
		}
		pop := pedPop(model, 5)
		runPedSim(t, model, pop, 10, 4)
		return utils.RNGState(), len(pop.IndData)
	}
	offState, offN := run(false)
	onState, onN := run(true)
	if string(offState) != string(onState) {
		t.Error("track_arg perturbed the RNG stream on a run that already draws meiosis masks")
	}
	if offN != onN {
		t.Errorf("different survivors: off=%d on=%d", offN, onN)
	}
}
