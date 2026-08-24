package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"testing"
)

// Regression test for a defect that sterilised every founding-couple scenario.
//
// Fitness is stored SCALED by mu_scale_factor: createChild writes
// int(fitness * mu_scale_factor) and birth.go divides it back out before using it as a
// birth probability. The eden, flood and creation setups passed a literal 1 instead, so
// their founders carried an effective fitness of 1e-6 and essentially never reproduced.
//
// It stayed hidden because selectionMode() collapses to "none" when track_mutations is
// off, and Fitness is then never read at all — every founding-couple smoke model had
// mutations off. The bug only appears in exactly the configuration the TMR4A study needs.

func fitnessModel(scale float64) *core.Model {
	return &core.Model{
		Parameters: map[string]float64{
			"mu_scale_factor": scale,
			"start_pop_size":  2,
			"lifespan":        650,
		},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"genome_bits": 1024, "seed": -1},
		ChromosomeArms: map[int]map[int][]int{1: {0: {0, 512}, 1: {512, 512}}},
	}
}

func emptyPop() *core.Pop {
	return &core.Pop{
		IndData:     map[int][]int{},
		Chromosomes: map[int][][]uint64{},
		Centromeres: map[int][2][]uint64{},
		Tracking:    map[string]int{},
	}
}

func TestFounderFitnessIsScaled(t *testing.T) {
	if got := FounderFitness(fitnessModel(1000000)); got != 1000000 {
		t.Errorf("FounderFitness = %d, want mu_scale_factor 1000000 — an unscaled 1 gives a "+
			"founder a birth probability of 1e-6", got)
	}
}

// A model with no mu_scale_factor must not yield 0, which would be worse than the bug:
// fitness 0 means no founder can ever reproduce at all.
func TestFounderFitnessFallsBackWhenScaleUnset(t *testing.T) {
	for _, scale := range []float64{0, -1} {
		m := fitnessModel(scale)
		if got := FounderFitness(m); got != 1 {
			t.Errorf("mu_scale_factor=%g: FounderFitness = %d, want the safe fallback 1", scale, got)
		}
	}
}

// The integration check: the eden couple must come out of setup carrying the scaled
// value, not the literal 1 the setups used to pass.
func TestEdenFoundersCarryScaledFitness(t *testing.T) {
	utils.SeedRNG(4242)
	model := fitnessModel(1000000)
	pop := emptyPop()
	if err := SetupPopEden(model, pop); err != nil {
		t.Fatalf("SetupPopEden: %v", err)
	}
	if len(pop.IndData) != 2 {
		t.Fatalf("eden built %d founders, want the single couple", len(pop.IndData))
	}
	for id, data := range pop.IndData {
		if got := data[individual.Fitness]; got != 1000000 {
			t.Errorf("founder %d has Fitness %d, want 1000000; birth.go divides by "+
				"mu_scale_factor, so a literal 1 here means the couple never reproduces "+
				"once track_mutations turns fitness selection on", id, got)
		}
	}
}

func TestCreationFoundersCarryScaledFitness(t *testing.T) {
	utils.SeedRNG(4242)
	model := fitnessModel(1000000)
	pop := emptyPop()
	if err := SetupPopCreation(model, pop); err != nil {
		t.Fatalf("SetupPopCreation: %v", err)
	}
	for id, data := range pop.IndData {
		if got := data[individual.Fitness]; got != 1000000 {
			t.Errorf("creation founder %d has Fitness %d, want 1000000", id, got)
		}
	}
}
