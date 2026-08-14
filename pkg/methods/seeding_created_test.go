package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"math"
	"testing"
)

// Created founder-allele seeding (TMR4A.md W4). The parameter under test is
// `founder_allele_divergence` — "the single most important parameter in the whole study" —
// so what matters is that the realized pairwise divergence between created alleles is the
// number the user asked for, not merely proportional to it.

func createdSeedModel(alleles int, divergence float64, genomeBits int) *core.Model {
	return &core.Model{
		Parameters: map[string]float64{
			"founder_alleles_per_locus": float64(alleles),
			"founder_allele_divergence": divergence,
		},
		StringParams:   map[string]string{"founder_allele_model": "created"},
		FreeParameters: map[string]int{"genome_bits": genomeBits, "seed": -1},
		ChromosomeArms: map[int]map[int][]int{1: {0: {0, genomeBits / 2}, 1: {genomeBits / 2, genomeBits / 2}}},
	}
}

func createdSeedPop(nFounders int) *core.Pop {
	pop := &core.Pop{
		IndData:     map[int][]int{},
		Chromosomes: map[int][][]uint64{},
		Centromeres: map[int][2][]uint64{},
	}
	for i := 0; i < nFounders; i++ {
		d := individual.MakeIndData()
		d[individual.Sex] = i % 2
		pop.IndData[i] = d
	}
	return pop
}

// hammingFrac is the fraction of the first `bits` positions at which two sequences differ.
func hammingFrac(a, b []uint64, bits int) float64 {
	diff := 0
	for i := 0; i < bits; i++ {
		x := (a[i/64] >> (uint(i) % 64)) & 1
		y := (b[i/64] >> (uint(i) % 64)) & 1
		if x != y {
			diff++
		}
	}
	return float64(diff) / float64(bits)
}

// The realized pairwise divergence must hit the requested d. Solving 2q(1−q)=d rather than
// taking q=d/2 is what makes this pass at larger d: the naive choice is low by 25% at
// d=0.5, which would silently mis-scale the study's headline dial.
func TestCreatedSeedingHitsRequestedDivergence(t *testing.T) {
	const bits = 20000
	for _, d := range []float64{0.01, 0.1, 0.25, 0.5} {
		utils.SeedRNG(1234)
		model := createdSeedModel(5, d, bits)
		pop := createdSeedPop(2)
		SeedingCreated(model, pop)

		if len(pop.CreatedAlleleSeqs) != 5 {
			t.Fatalf("d=%g: built %d alleles, want 5", d, len(pop.CreatedAlleleSeqs))
		}
		sum, pairs := 0.0, 0
		for i := 0; i < 5; i++ {
			for j := i + 1; j < 5; j++ {
				sum += hammingFrac(pop.CreatedAlleleSeqs[i], pop.CreatedAlleleSeqs[j], bits)
				pairs++
			}
		}
		got := sum / float64(pairs)
		// Sampling error over 10 pairs × 20k sites is well under 1%.
		if math.Abs(got-d) > 0.01 {
			t.Errorf("founder_allele_divergence=%g: realized mean pairwise divergence %.4f", d, got)
		}
	}
}

// d = 0 is the null arm: the created alleles are identical, so there is no created
// divergence anywhere and every difference a run accumulates is mutational.
func TestCreatedSeedingZeroDivergence(t *testing.T) {
	const bits = 4096
	utils.SeedRNG(99)
	model := createdSeedModel(4, 0, bits)
	pop := createdSeedPop(2)
	before := utils.RNGState()
	SeedingCreated(model, pop)

	for i := 1; i < 4; i++ {
		if f := hammingFrac(pop.CreatedAlleleSeqs[0], pop.CreatedAlleleSeqs[i], bits); f != 0 {
			t.Errorf("allele %d differs from allele 0 at %.4f of sites under d=0", i, f)
		}
	}
	if string(before) != string(utils.RNGState()) {
		t.Error("d=0 drew RNG; with no divergence to place there is nothing to draw")
	}
}

// The founders' somatic genomes must be the first two alleles of their own pools —
// ordinary diploid individuals whose germline happens to hold more.
func TestCreatedSeedingFounderGenomesComeFromTheirPool(t *testing.T) {
	const bits = 4096
	utils.SeedRNG(5)
	model := createdSeedModel(10, 0.2, bits)
	pop := createdSeedPop(2)
	SeedingCreated(model, pop)

	for _, id := range []int{0, 1} {
		pool := pop.CreatedPool(id)
		if len(pool) != 5 {
			t.Fatalf("founder %d pool %v, want 5 alleles", id, pool)
		}
		chrom, ok := pop.Chromosomes[id]
		if !ok {
			t.Fatalf("founder %d has no genome", id)
		}
		for s := 0; s < 2; s++ {
			if f := hammingFrac(chrom[s], pop.CreatedAlleleSeqs[pool[s]], bits); f != 0 {
				t.Errorf("founder %d strand %d differs from pool allele %d at %.4f of sites",
					id, s, pool[s], f)
			}
		}
		// The two somatic strands must be genuinely different sequences at d>0, or the
		// founder is not the heterozygote the model intends.
		if f := hammingFrac(chrom[0], chrom[1], bits); f < 0.1 {
			t.Errorf("founder %d somatic strands differ at only %.4f of sites", id, f)
		}
		if got := pop.IndData[id][individual.AlleleCount]; got == 0 {
			t.Errorf("founder %d has AlleleCount 0; meiosis would skip it", id)
		}
	}
}

// Divergence above 0.5 is not reachable with independent per-site draws, so it is capped
// rather than producing a NaN q.
func TestCreatedSeedingCapsDivergence(t *testing.T) {
	const bits = 4096
	utils.SeedRNG(7)
	model := createdSeedModel(3, 0.9, bits)
	pop := createdSeedPop(2)
	SeedingCreated(model, pop)

	f := hammingFrac(pop.CreatedAlleleSeqs[0], pop.CreatedAlleleSeqs[1], bits)
	if math.IsNaN(f) || f < 0.4 || f > 0.6 {
		t.Errorf("d=0.9 capped at 0.5: realized divergence %.4f, want ≈0.5", f)
	}
}

// Deterministic under a fixed seed — the created genome is part of the run's reproducible
// state, not a per-invocation surprise.
func TestCreatedSeedingIsDeterministic(t *testing.T) {
	const bits = 2048
	build := func() *core.Pop {
		utils.SeedRNG(31415)
		model := createdSeedModel(6, 0.15, bits)
		pop := createdSeedPop(4)
		SeedingCreated(model, pop)
		return pop
	}
	a, b := build(), build()
	for i := range a.CreatedAlleleSeqs {
		if hammingFrac(a.CreatedAlleleSeqs[i], b.CreatedAlleleSeqs[i], bits) != 0 {
			t.Fatalf("created allele %d differs between identically seeded runs", i)
		}
	}
	for id := 0; id < 4; id++ {
		pa, pb := a.CreatedPool(id), b.CreatedPool(id)
		if len(pa) != len(pb) {
			t.Fatalf("founder %d pool sizes differ", id)
		}
		for i := range pa {
			if pa[i] != pb[i] {
				t.Fatalf("founder %d pool differs: %v vs %v", id, pa, pb)
			}
		}
	}
}
