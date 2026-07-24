package simulation

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"testing"
)

// linkageModel builds a one-chromosome model (p arm [0,100), q arm [100,200),
// genome_bits=200) under the given linkage_model, used to drive the real meiosis
// mask through createMask.
func linkageModel(linkage string) *core.Model {
	m := &core.Model{
		Parameters:     map[string]float64{},
		StringParams:   map[string]string{"linkage_model": linkage},
		FreeParameters: map[string]int{"genome_bits": 200},
		ChromosomeArms: map[int]map[int][]int{
			1: {0: {0, 100}, 1: {100, 100}},
		},
	}
	return m
}

// runInheritTrials drives one parent's gamete formation `trials` times through the
// real createMask (mask for copy 0), then applies BOTH meiosis (founder bitfield)
// and InheritMutations (de-novo pool) with that same mask. The parent carries a
// founder bit at each founderPos on copy 0 and a strand-0 mutation at each mutPos.
// Returns, per mutation index, the fraction of trials the mutation was inherited
// on the child gamete, and — for indices with a matching founder bit — the
// fraction where mutation and founder bit AGREE (both inherited or both not).
func runInheritTrials(t *testing.T, model *core.Model, mutPos []int, founderPos []int, trials int) (inherited []float64, agreeWithBit []float64) {
	utils.SeedRNG(1)
	words := (model.FreeParameters["genome_bits"] + 63) / 64
	inheritedCnt := make([]int, len(mutPos))
	agreeCnt := make([]int, len(mutPos))

	// mutation ids 1..len(mutPos), pool positions from mutPos.
	pool := map[int]core.Mutation{}
	strand0 := make([]int, len(mutPos))
	for i, p := range mutPos {
		id := i + 1
		pool[id] = core.Mutation{Id: id, Position: p}
		strand0[i] = id
	}

	for tr := 0; tr < trials; tr++ {
		pop := &core.Pop{
			IndData:      map[int][]int{1: individual.MakeIndData(), 2: individual.MakeIndData()},
			IndMutations: map[int]map[int][]int{1: {0: append([]int{}, strand0...), 1: {}}, 2: {0: {}, 1: {}}},
			MutationPool: pool,
			Chromosomes:  map[int][][]uint64{1: {make([]uint64, words), make([]uint64, words)}, 2: {make([]uint64, words), make([]uint64, words)}},
		}
		for _, fp := range founderPos {
			pop.Chromosomes[1][0][fp/64] |= 1 << (uint(fp) % 64)
		}

		mask, _ := createMask(model, 0)
		utils.InheritMutations(model, pop, mask, 1, 2, 0)
		meiosis(pop, mask, 1, 2, 0)

		has := map[int]bool{}
		for _, id := range pop.IndMutations[2][0] {
			has[id] = true
		}
		for i := range mutPos {
			if has[i+1] {
				inheritedCnt[i]++
			}
			// founder bit at the same position (founderPos parallels mutPos here)
			if i < len(founderPos) {
				fp := founderPos[i]
				bit := (pop.Chromosomes[2][0][fp/64]>>(uint(fp)%64))&1 == 1
				if bit == has[i+1] {
					agreeCnt[i]++
				}
			}
		}
	}
	inherited = make([]float64, len(mutPos))
	agreeWithBit = make([]float64, len(mutPos))
	for i := range mutPos {
		inherited[i] = float64(inheritedCnt[i]) / float64(trials)
		agreeWithBit[i] = float64(agreeCnt[i]) / float64(trials)
	}
	return
}

// TestLinkageBitfieldReconciliation is the requirement-(b) check: under
// linkage_model="arm" a de-novo mutation co-segregates with a founder bitfield
// allele at the SAME position (they always agree), whereas under the legacy
// (inverted-polarity) default they always disagree. This is the §6h
// "de-novo mutations bypass the bitfield" defect and its opt-in fix.
func TestLinkageBitfieldReconciliation(t *testing.T) {
	// One mutation and one founder bit, both at position 10.
	pos := []int{10}

	_, agreeLegacy := runInheritTrials(t, linkageModel("legacy"), pos, pos, 5000)
	if agreeLegacy[0] > 0.001 {
		t.Errorf("legacy: mutation agrees with founder bit %.4f of trials, want ~0 (anti-segregation)", agreeLegacy[0])
	}

	_, agreeArm := runInheritTrials(t, linkageModel("arm"), pos, pos, 5000)
	if agreeArm[0] < 0.999 {
		t.Errorf("arm: mutation agrees with founder bit %.4f of trials, want ~1 (co-segregation)", agreeArm[0])
	}
}

// TestLinkageSameArmDecay is the requirement-(a) check: mutations sharing an arm
// are physically linked and recombination breaks that linkage with distance. Two
// close mutations co-inherit far more often than a distant pair. (Polarity does
// not change co-segregation AMONG mutations, so this holds under either model;
// tested under "arm".)
func TestLinkageSameArmDecay(t *testing.T) {
	// positions: reference 10, close 15, far 90 — all on the p arm [0,100).
	mutPos := []int{10, 15, 90}
	inh, _ := runInheritTrials(t, linkageModel("arm"), mutPos, nil, 40000)

	// Co-inheritance with the reference = fraction of trials where the two share
	// the same inherited/not state. Recover it from a joint count instead: rerun
	// measuring agreement between mutation 0 and each other.
	// (Simpler: derive from a dedicated joint pass.)
	closeCo, farCo := jointCoInherit(t, linkageModel("arm"), mutPos)
	_ = inh
	if closeCo < 0.9 {
		t.Errorf("close same-arm pair co-inherit %.3f, want > 0.9 (tight linkage)", closeCo)
	}
	if farCo > closeCo-0.1 {
		t.Errorf("far pair co-inherit %.3f not detectably less than close %.3f (recombination should break linkage)", farCo, closeCo)
	}
}

// jointCoInherit measures the fraction of trials in which mutation[0] and
// mutation[1] (close) / mutation[2] (far) end up in the same inherited state.
func jointCoInherit(t *testing.T, model *core.Model, mutPos []int) (closeCo, farCo float64) {
	utils.SeedRNG(7)
	words := (model.FreeParameters["genome_bits"] + 63) / 64
	pool := map[int]core.Mutation{}
	strand0 := make([]int, len(mutPos))
	for i, p := range mutPos {
		pool[i+1] = core.Mutation{Id: i + 1, Position: p}
		strand0[i] = i + 1
	}
	const trials = 40000
	var closeAgree, farAgree int
	for tr := 0; tr < trials; tr++ {
		pop := &core.Pop{
			IndData:      map[int][]int{1: individual.MakeIndData(), 2: individual.MakeIndData()},
			IndMutations: map[int]map[int][]int{1: {0: append([]int{}, strand0...), 1: {}}, 2: {0: {}, 1: {}}},
			MutationPool: pool,
			Chromosomes:  map[int][][]uint64{1: {make([]uint64, words), make([]uint64, words)}, 2: {make([]uint64, words), make([]uint64, words)}},
		}
		mask, _ := createMask(model, 0)
		utils.InheritMutations(model, pop, mask, 1, 2, 0)
		has := map[int]bool{}
		for _, id := range pop.IndMutations[2][0] {
			has[id] = true
		}
		if has[1] == has[2] {
			closeAgree++
		}
		if has[1] == has[3] {
			farAgree++
		}
	}
	return float64(closeAgree) / trials, float64(farAgree) / trials
}
