package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"fmt"
	"math"
	"sort"
)

// Created founder-allele seeding — TMR4A.md work item W4, the seeding half.
// Registered as seed_style "created"; see pkg/core/created.go for the model itself and
// pkg/simulation/created_meiosis.go for the transmission half.
//
// WHAT IT BUILDS. *A* = `founder_alleles_per_locus` created HAPLOTYPES, each a full-genome
// bitfield, any two differing at an expected fraction *d* = `founder_allele_divergence` of
// sites. Then a pool of those alleles per founder, and each founder's somatic genome from
// the first two of its own pool.
//
// WHY DIVERGENCE IS SET HERE AND NOT ACCUMULATED. This is the dial TMR4A.md calls "the
// single most important parameter in the whole study". The created alleles differ by
// construction, at t=0, with no elapsed time and no mutations. A coalescent method reading
// those sequences has no way to know that: it converts *d* into an age through the assumed
// mutation rate and reports ~*d*/(2μ) generations. DRIFT knows the true answer is zero
// generations, which is precisely the comparison W7 is built to make.
//
// THE SITE MODEL. Each created allele draws each bit independently with probability *q*,
// chosen so the expected pairwise difference between any two alleles is exactly *d*:
//
//	2q(1−q) = d   ⇒   q = (1 − √(1 − 2d)) / 2
//
// Solving for q rather than setting q = d/2 matters because the naive choice gives
// d(1 − d/2), which is wrong by 25% at d = 0.5. The model is symmetric — no allele is a
// reference and no pair is closer than another — so the created alleles form a star with
// no internal structure for an inference method to mistake for a genealogy. d is capped at
// 0.5, the maximum pairwise difference independent per-site draws can produce.
//
// RNG ORDER is allele-major then bit-minor over sorted founder ids, documented here so the
// stream is reproducible. This seeder draws RNG — it is a genome seeder, that is its job —
// but it runs only when explicitly selected, so no default run is affected.
func SeedingCreated(model *core.Model, pop *core.Pop) {
	model.FreeParameters["seed"] = 1

	totalBits := model.FreeParameters["genome_bits"]
	if totalBits <= 0 {
		fmt.Println("   Created seeding skipped: genome_bits is 0")
		return
	}
	words := (totalBits + 63) / 64

	alleles := int(model.Parameters["founder_alleles_per_locus"])
	if alleles < 1 {
		alleles = 4
	}
	d := model.Parameters["founder_allele_divergence"]
	if d < 0 {
		d = 0
	}
	if d > 0.5 {
		fmt.Printf("   founder_allele_divergence %.4f capped at 0.5 (the maximum pairwise "+
			"difference independent per-site draws can produce)\n", d)
		d = 0.5
	}
	q := (1 - math.Sqrt(1-2*d)) / 2

	// Build the created haplotypes. Allele-major, bit-minor.
	pop.CreatedAlleleSeqs = make([][]uint64, alleles)
	for a := 0; a < alleles; a++ {
		seq := make([]uint64, words)
		if q > 0 {
			for bit := 0; bit < totalBits; bit++ {
				if utils.RandFloat64() < q {
					seq[bit/64] |= 1 << (uint(bit) % 64)
				}
			}
		}
		pop.CreatedAlleleSeqs[a] = seq
	}

	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	sort.Ints(ids)
	pop.AssignCreatedPools(ids, alleles)

	// Each founder's SOMATIC genome is the first two alleles of its own germline pool —
	// an ordinary diploid individual whose germline happens to hold more. A one-allele
	// pool gives a founder homozygous for that created allele.
	for _, id := range ids {
		pool := pop.CreatedPool(id)
		second := 0
		if len(pool) > 1 {
			second = 1
		}
		c0 := append([]uint64(nil), pop.CreatedAlleleSeqs[pool[0]]...)
		c1 := append([]uint64(nil), pop.CreatedAlleleSeqs[pool[second]]...)
		pop.Chromosomes[id] = [][]uint64{c0, c1}

		var blankCentromeres [2][]uint64
		blankCentromeres[0] = make([]uint64, (len(model.ChromosomeArms)+63)/64)
		blankCentromeres[1] = make([]uint64, (len(model.ChromosomeArms)+63)/64)
		pop.Centromeres[id] = blankCentromeres

		pop.IndData[id][individual.YGens] = 0
		pop.IndData[id][individual.MtGens] = 0
		pop.IndData[id][individual.MinGenealoGens] = 0
		pop.IndData[id][individual.MaxGenealoGens] = 0
		pop.IndData[id][individual.NumCentromeres] = 0
		pop.IndData[id][individual.AlleleCount] =
			utils.CountSetBits(c0) + utils.CountSetBits(c1)
	}

	fmt.Printf("   Created %d founder alleles over %d founders (pairwise divergence d=%.4g, "+
		"per-site q=%.4g); germline pools of up to %d alleles per founder\n",
		alleles, len(ids), d, q, (alleles+len(ids)-1)/max2(len(ids), 1))
}

func max2(a, b int) int {
	if a > b {
		return a
	}
	return b
}
