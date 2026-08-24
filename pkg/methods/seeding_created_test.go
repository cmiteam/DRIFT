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

// The "exclusive" site model (founder_allele_site_model). Its defining property is that no
// two created alleles ever carry the derived state at the same site, so an observed shared
// state downstream is genuine shared inheritance rather than coincidence. Both models must
// still hit the requested pairwise divergence — that is the dial the study turns.

func createdSiteModelPop(alleles int, divergence float64, genomeBits int, siteModel string) *core.Model {
	m := createdSeedModel(alleles, divergence, genomeBits)
	m.StringParams["founder_allele_site_model"] = siteModel
	return m
}

// carriersAt counts how many created alleles hold the derived state at one bit.
func carriersAt(seqs [][]uint64, bit int) int {
	n := 0
	for _, seq := range seqs {
		if (seq[bit/64]>>(uint(bit)%64))&1 == 1 {
			n++
		}
	}
	return n
}

func TestCreatedSeedingExclusiveHasNoRecurrentStates(t *testing.T) {
	const bits, alleles = 20000, 5
	utils.SeedRNG(4242)
	pop := createdSeedPop(2)
	SeedingCreated(createdSiteModelPop(alleles, 0.2, bits, "exclusive"), pop)

	variable := 0
	for bit := 0; bit < bits; bit++ {
		switch n := carriersAt(pop.CreatedAlleleSeqs, bit); {
		case n == 0:
		case n == 1:
			variable++
		default:
			t.Fatalf("bit %d is derived in %d alleles; the exclusive model permits at most 1 "+
				"— that coincidence is exactly what it exists to remove", bit, n)
		}
	}
	// p = A*d/2 = 5*0.2/2 = 0.5 of sites variable.
	if got := float64(variable) / bits; math.Abs(got-0.5) > 0.02 {
		t.Errorf("variable-site fraction %.4f, want ~0.5 (= A*d/2)", got)
	}
}

// The independent model is expected to produce recurrent states — this pins the contrast
// so the two models cannot be silently conflated.
func TestCreatedSeedingIndependentDoesProduceRecurrentStates(t *testing.T) {
	const bits, alleles = 20000, 5
	utils.SeedRNG(4242)
	pop := createdSeedPop(2)
	SeedingCreated(createdSiteModelPop(alleles, 0.2, bits, "independent"), pop)

	for bit := 0; bit < bits; bit++ {
		if carriersAt(pop.CreatedAlleleSeqs, bit) > 1 {
			return // found one, as expected
		}
	}
	t.Error("independent model produced no site derived in 2+ alleles; at q~0.1 over 20k sites "+
		"that is essentially impossible and suggests the model is not drawing independently")
}

func TestCreatedSeedingExclusiveHitsRequestedDivergence(t *testing.T) {
	const bits, alleles = 20000, 5
	for _, d := range []float64{0.01, 0.1, 0.25, 0.4} { // ceiling is 2/A = 0.4
		utils.SeedRNG(1234)
		pop := createdSeedPop(2)
		SeedingCreated(createdSiteModelPop(alleles, d, bits, "exclusive"), pop)

		sum, pairs := 0.0, 0
		for i := 0; i < alleles; i++ {
			for j := i + 1; j < alleles; j++ {
				sum += hammingFrac(pop.CreatedAlleleSeqs[i], pop.CreatedAlleleSeqs[j], bits)
				pairs++
			}
		}
		if got := sum / float64(pairs); math.Abs(got-d) > 0.01 {
			t.Errorf("exclusive d=%g: realized mean pairwise divergence %.4f", d, got)
		}
	}
}

// The exclusive model cannot express d > 2/A, because every variable site is spent on a
// single allele. It must cap rather than silently produce something else.
func TestCreatedSeedingExclusiveCapsDivergenceAtTwoOverA(t *testing.T) {
	const bits, alleles = 20000, 10 // ceiling 2/10 = 0.2
	utils.SeedRNG(7)
	pop := createdSeedPop(2)
	SeedingCreated(createdSiteModelPop(alleles, 0.5, bits, "exclusive"), pop)

	sum, pairs := 0.0, 0
	for i := 0; i < alleles; i++ {
		for j := i + 1; j < alleles; j++ {
			sum += hammingFrac(pop.CreatedAlleleSeqs[i], pop.CreatedAlleleSeqs[j], bits)
			pairs++
		}
	}
	if got := sum / float64(pairs); math.Abs(got-0.2) > 0.01 {
		t.Errorf("requested d=0.5 with A=10: realized divergence %.4f, want it capped to 0.2", got)
	}
}

// An unset site model must reproduce the original behaviour exactly, since every existing
// created model omits the parameter.
func TestCreatedSeedingDefaultsToIndependent(t *testing.T) {
	const bits, alleles = 8192, 4
	utils.SeedRNG(31337)
	popDefault := createdSeedPop(2)
	SeedingCreated(createdSeedModel(alleles, 0.1, bits), popDefault)

	utils.SeedRNG(31337)
	popExplicit := createdSeedPop(2)
	SeedingCreated(createdSiteModelPop(alleles, 0.1, bits, "independent"), popExplicit)

	for a := range popDefault.CreatedAlleleSeqs {
		for w := range popDefault.CreatedAlleleSeqs[a] {
			if popDefault.CreatedAlleleSeqs[a][w] != popExplicit.CreatedAlleleSeqs[a][w] {
				t.Fatalf("allele %d word %d differs: unset site model must equal \"independent\"", a, w)
			}
		}
	}
}
