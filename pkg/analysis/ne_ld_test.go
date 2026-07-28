package analysis

import (
	"math"
	"testing"

	"drift/pkg/core"
)

// makeNe2ChromModel: two chromosomes (chrom 0 = bits [0,64), chrom 1 = bits
// [64,128)), so sites 0 and 64 lie on different chromosomes (unlinked).
func makeNe2ChromModel() *core.Model {
	return &core.Model{
		ChromosomeArms: map[int]map[int][]int{
			0: {0: {0, 32}, 1: {32, 32}},  // chrom 0 fills bits 0..63
			1: {0: {64, 32}, 1: {96, 32}}, // chrom 1 fills bits 64..127
		},
		FreeParameters: map[string]int{"genome_bits": 128, "run": 1, "year": 10},
		Parameters:     map[string]float64{"ld_max_pairs": 10000, "ne_cM_per_bit": 1.0},
		ModelName:      "Ne2Chrom",
	}
}

// TestRecombFraction covers the three regimes of the bit-distance -> recombination
// fraction map used by the LD leg.
func TestRecombFraction(t *testing.T) {
	model := makeNe2ChromModel()

	// Different chromosomes => unlinked.
	if c := recombFraction(model, nil, 0, 64, 1.0); math.Abs(c-0.5) > 1e-12 {
		t.Errorf("cross-chromosome c = %v, want 0.5", c)
	}

	// Same chromosome, Haldane: dist 10 bits at 1 cM/bit = 0.1 Morgan.
	wantClose := 0.5 * (1 - math.Exp(-0.2))
	if c := recombFraction(model, nil, 0, 10, 1.0); math.Abs(c-wantClose) > 1e-12 {
		t.Errorf("same-chromosome c = %v, want %v", c, wantClose)
	}

	// Monotone increasing with distance, saturating below 0.5.
	cNear := recombFraction(model, nil, 0, 5, 1.0)
	cFar := recombFraction(model, nil, 0, 40, 1.0)
	if !(cNear < cFar && cFar <= 0.5) {
		t.Errorf("expected cNear < cFar <= 0.5, got %v, %v", cNear, cFar)
	}
}

// TestRecombFraction_UsesGeneticMap checks that when a real genetic map is present
// (roadmap §6a) recombFraction reads calibrated map distance instead of the
// ne_cM_per_bit scalar. The map makes the p arm 0.5 cM/bit, so a 10-bit gap is 5 cM
// (0.05 Morgan) — deliberately different from the 1.0 cM/bit scalar (0.1 Morgan).
func TestRecombFraction_UsesGeneticMap(t *testing.T) {
	arms := map[int]map[int][]int{0: {0: {0, 100}, 1: {100, 100}}}
	armCM := map[int]map[int]float64{0: {0: 50, 1: 150}} // p arm 0.5 cM/bit
	model := &core.Model{
		ChromosomeArms: arms,
		ArmCM:          armCM,
		FreeParameters: map[string]int{"genome_bits": 200},
		Parameters:     map[string]float64{},
	}
	gm := core.BuildGeneticMap(arms, armCM, 1.0)

	// Map path: 10 bits at 0.5 cM/bit = 5 cM = 0.05 Morgan.
	wantMap := 0.5 * (1 - math.Exp(-2*0.05))
	if c := recombFraction(model, gm, 10, 20, 1.0); math.Abs(c-wantMap) > 1e-12 {
		t.Errorf("map-calibrated c = %v, want %v", c, wantMap)
	}
	// Scalar fallback (nil map) at the same cMPerBit gives the DIFFERENT 0.1-Morgan value.
	wantScalar := 0.5 * (1 - math.Exp(-2*0.10))
	if c := recombFraction(model, nil, 10, 20, 1.0); math.Abs(c-wantScalar) > 1e-12 {
		t.Errorf("scalar-fallback c = %v, want %v", c, wantScalar)
	}
	if math.Abs(wantMap-wantScalar) < 1e-6 {
		t.Fatal("test misconfigured: map and scalar should differ")
	}
}

// TestLdRecentNe_Disabled: with ne_cM_per_bit = 0 the LD leg is a no-op (no
// recombination map => no defensible calibration).
func TestLdRecentNe_Disabled(t *testing.T) {
	model := makeNe2ChromModel()
	model.Parameters["ne_cM_per_bit"] = 0

	// Two perfectly-linked unlinked-locus haplotypes across 4 individuals.
	pop := ldTestPop()
	ne, pairs := ldRecentNe(model, pop, []int{1, 2, 3, 4})
	if ne != 0 || pairs != 0 {
		t.Errorf("disabled LD leg returned (%v, %d), want (0, 0)", ne, pairs)
	}
}

// TestLdRecentNe_Configured: with a recombination knob set and strong LD between two
// unlinked sites, the estimator returns a finite positive recent Ne over >= 1 pair.
func TestLdRecentNe_Configured(t *testing.T) {
	model := makeNe2ChromModel() // ne_cM_per_bit = 1.0
	pop := ldTestPop()

	ne, pairs := ldRecentNe(model, pop, []int{1, 2, 3, 4})
	if pairs < 1 {
		t.Fatalf("expected >= 1 pair used, got %d", pairs)
	}
	if !(ne > 0) || math.IsInf(ne, 0) || math.IsNaN(ne) {
		t.Errorf("expected finite positive Ne, got %v", ne)
	}
}

// ldTestPop builds 4 individuals in which bit 0 (chrom 0) and bit 64 (chrom 1) are
// perfectly correlated across strand copies (r^2 = 1), both segregating.
func ldTestPop() *core.Pop {
	// strand with both bit0 and bit64 set ("11").
	both := []uint64{uint64(1), uint64(1)}
	none := []uint64{0, 0}
	pop := &core.Pop{
		IndData:     map[int][]int{},
		Chromosomes: map[int][][]uint64{},
	}
	// ind1,2: strand0 = 11, strand1 = 00 ; ind3: strand0 = 00, strand1 = 11 ; ind4: 00/00.
	pop.Chromosomes[1] = [][]uint64{both, none}
	pop.Chromosomes[2] = [][]uint64{both, none}
	pop.Chromosomes[3] = [][]uint64{none, both}
	pop.Chromosomes[4] = [][]uint64{none, none}
	for id := 1; id <= 4; id++ {
		pop.IndData[id] = []int{}
	}
	return pop
}
