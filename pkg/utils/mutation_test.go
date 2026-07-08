package utils

import (
	"math"
	"testing"

	"drift/pkg/core"
)

// newFitnessPop builds a minimal Pop with the given per-strand mutation ids for
// individual 1 and a mutation pool from (id -> effect, dominance-percent).
func newFitnessPop(strand0, strand1 []int, pool map[int]core.Mutation) *core.Pop {
	return &core.Pop{
		IndMutations: map[int]map[int][]int{
			1: {0: strand0, 1: strand1},
		},
		MutationPool: pool,
	}
}

// modelWithFitness returns a model configured for the given fitness_model and
// (for synergistic) epistasis coefficient.
func modelWithFitness(fitnessModel string, beta float64) *core.Model {
	return &core.Model{
		StringParams: map[string]string{"fitness_model": fitnessModel},
		Parameters:   map[string]float64{"epistasis_coefficient": beta},
	}
}

func approx(a, b float64) bool { return math.Abs(a-b) < 1e-12 }

// TestCountFitnessAndMutations checks the per-locus zygosity/dominance scoring
// under the default additive model (fitness = 1 + load), i.e. that the combined
// fitness equals 1 + the expected additive load.
func TestCountFitnessAndMutations(t *testing.T) {
	// Effect is negative for a deleterious mutation; s is its homozygous effect.
	const s = -0.10

	cases := []struct {
		name       string
		strand0    []int
		strand1    []int
		pool       map[int]core.Mutation
		wantCopies int
		wantLoad   float64
	}{
		{
			name:       "no mutations",
			strand0:    nil,
			strand1:    nil,
			pool:       map[int]core.Mutation{},
			wantCopies: 0,
			wantLoad:   0,
		},
		{
			name:    "heterozygous additive h=0.5 -> half effect",
			strand0: []int{7},
			strand1: nil,
			pool: map[int]core.Mutation{
				7: {Id: 7, Effect: s, Dominance: 50},
			},
			wantCopies: 1,
			wantLoad:   0.5 * s,
		},
		{
			name:    "homozygous -> full effect regardless of h",
			strand0: []int{7},
			strand1: []int{7},
			pool: map[int]core.Mutation{
				7: {Id: 7, Effect: s, Dominance: 50},
			},
			wantCopies: 2,
			wantLoad:   s,
		},
		{
			name:    "heterozygous fully recessive h=0 -> no contribution",
			strand0: []int{7},
			strand1: nil,
			pool: map[int]core.Mutation{
				7: {Id: 7, Effect: s, Dominance: 0},
			},
			wantCopies: 1,
			wantLoad:   0,
		},
		{
			name:    "heterozygous fully dominant h=1 -> full effect",
			strand0: nil,
			strand1: []int{7},
			pool: map[int]core.Mutation{
				7: {Id: 7, Effect: s, Dominance: 100},
			},
			wantCopies: 1,
			wantLoad:   s,
		},
		{
			name:    "multiple loci sum additively (het additive + homozygous)",
			strand0: []int{7, 9},
			strand1: []int{9},
			pool: map[int]core.Mutation{
				7: {Id: 7, Effect: s, Dominance: 50}, // het -> 0.5s
				9: {Id: 9, Effect: s, Dominance: 50}, // hom -> s
			},
			wantCopies: 3,
			wantLoad:   0.5*s + s,
		},
		{
			name:       "mutation missing from pool is skipped in load but counted as a copy",
			strand0:    []int{42},
			strand1:    nil,
			pool:       map[int]core.Mutation{},
			wantCopies: 1,
			wantLoad:   0,
		},
	}

	additive := modelWithFitness(FitAdditive, 0)
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			pop := newFitnessPop(tc.strand0, tc.strand1, tc.pool)
			gotCopies, gotFitness := CountFitnessAndMutations(pop, 1, additive)
			if gotCopies != tc.wantCopies {
				t.Errorf("copies = %d, want %d", gotCopies, tc.wantCopies)
			}
			if !approx(gotFitness, 1+tc.wantLoad) {
				t.Errorf("fitness = %.12g, want %.12g", gotFitness, 1+tc.wantLoad)
			}
		})
	}
}

// TestCombineFitnessModels checks that the three fitness models combine per-locus
// contributions as documented, that additive is unchanged, and — critically for
// the §6h neutral harness — that every model returns exactly 1.0 on zero load.
func TestCombineFitnessModels(t *testing.T) {
	contribs := []float64{-0.10, -0.20} // two deleterious loci
	// additive: 1 + (-0.30) = 0.70
	// multiplicative: (1-0.10)(1-0.20) = 0.72
	// synergistic β=1: L=-0.30 -> 1 + L + β·L·|L| = 0.70 - 0.09 = 0.61
	cases := []struct {
		name  string
		model *core.Model
		want  float64
	}{
		{"additive", modelWithFitness(FitAdditive, 0), 0.70},
		{"multiplicative", modelWithFitness(FitMultiplicative, 0), 0.72},
		{"synergistic beta=1", modelWithFitness(FitSynergistic, 1), 0.61},
		{"synergistic beta=0 reduces to additive", modelWithFitness(FitSynergistic, 0), 0.70},
		{"unknown model falls back to additive", modelWithFitness("bogus", 0), 0.70},
		{"empty string defaults to additive", &core.Model{Parameters: map[string]float64{}}, 0.70},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := CombineFitness(contribs, tc.model)
			if !approx(got, tc.want) {
				t.Errorf("CombineFitness = %.12g, want %.12g", got, tc.want)
			}
		})
	}
}

// TestCombineFitnessNeutralNoOp guarantees the compatibility constraint: with no
// load (empty or all-zero contributions) every model returns exactly 1.0, so a
// strictly neutral run stays byte-identical regardless of fitness_model.
func TestCombineFitnessNeutralNoOp(t *testing.T) {
	for _, fm := range []string{FitAdditive, FitMultiplicative, FitSynergistic} {
		model := modelWithFitness(fm, 1)
		for _, contribs := range [][]float64{nil, {}, {0.0, 0.0, 0.0}} {
			got := CombineFitness(contribs, model)
			if got != 1.0 {
				t.Errorf("model %s, contribs %v: fitness = %v, want exactly 1.0", fm, contribs, got)
			}
		}
	}
}

// TestFitnessDominanceDefaultsToAdditive checks that GenerateNewMutations stamps
// h=0.5 (additive) when the fitness_dominance param is absent, rather than the
// bare-map-lookup zero (which would be fully recessive).
func TestFitnessDominanceDefaultsToAdditive(t *testing.T) {
	SeedRNG(12345)
	model := &core.Model{
		Parameters: map[string]float64{
			"mu":              5, // guarantee at least some draws
			"f_neutral":       0, // every mutation non-neutral so Effect != 0
			"f_beneficial":    0,
			"shape":           1,
			"scale":           0.05,
			"Weibull_adj":     1,
			"mu_scale_factor": 1,
			// fitness_dominance intentionally omitted
		},
		FreeParameters: map[string]int{"mutID": 0, "genome_bits": 1000},
	}
	pop := &core.Pop{
		IndMutations: map[int]map[int][]int{},
		MutationPool: map[int]core.Mutation{},
		MutationHist: map[int]int{},
	}
	GenerateNewMutations(model, pop, 1)
	if len(pop.MutationPool) == 0 {
		t.Fatal("expected at least one mutation generated")
	}
	for id, m := range pop.MutationPool {
		if m.Dominance != 50 {
			t.Errorf("mutation %d Dominance = %d, want 50 (additive default)", id, m.Dominance)
		}
	}
}
