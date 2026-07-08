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

func approx(a, b float64) bool { return math.Abs(a-b) < 1e-12 }

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
				7: {Id: 7, Effect: s, Dominance: 50},  // het -> 0.5s
				9: {Id: 9, Effect: s, Dominance: 50},  // hom -> s
			},
			wantCopies: 3,
			wantLoad:   0.5*s + s,
		},
		{
			name:    "mutation missing from pool is skipped in load but counted as a copy",
			strand0: []int{42},
			strand1: nil,
			pool:    map[int]core.Mutation{},
			wantCopies: 1,
			wantLoad:   0,
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			pop := newFitnessPop(tc.strand0, tc.strand1, tc.pool)
			gotCopies, gotLoad := CountFitnessAndMutations(pop, 1)
			if gotCopies != tc.wantCopies {
				t.Errorf("copies = %d, want %d", gotCopies, tc.wantCopies)
			}
			if !approx(gotLoad, tc.wantLoad) {
				t.Errorf("load = %.12g, want %.12g", gotLoad, tc.wantLoad)
			}
		})
	}
}

// TestFitnessDominanceDefaultsToAdditive checks that GenerateNewMutations stamps
// h=0.5 (additive) when the fitness_dominance param is absent, rather than the
// bare-map-lookup zero (which would be fully recessive).
func TestFitnessDominanceDefaultsToAdditive(t *testing.T) {
	SeedRNG(12345)
	model := &core.Model{
		Parameters: map[string]float64{
			"mu":          5, // guarantee at least some draws
			"f_neutral":   0, // every mutation non-neutral so Effect != 0
			"f_beneficial": 0,
			"shape":       1,
			"scale":       0.05,
			"Weibull_adj": 1,
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
