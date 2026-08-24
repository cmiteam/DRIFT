package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"fmt"
	"math"
	"sort"
	"strings"
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
// THE SITE MODEL, selected by `founder_allele_site_model`. Both produce the same expected
// pairwise divergence d; they differ in whether two alleles may carry the derived state at
// the SAME site, which is the identity-by-state / identity-by-descent question W3 exists
// to keep straight.
//
// "independent" (default, the original model). Each created allele draws each bit
// independently with probability q, chosen so the expected pairwise difference between any
// two alleles is exactly d:
//
//	2q(1−q) = d   ⇒   q = (1 − √(1 − 2d)) / 2
//
// Solving for q rather than setting q = d/2 matters because the naive choice gives
// d(1 − d/2), which is wrong by 25% at d = 0.5. Here two alleles DO coincide at a derived
// site by chance, at rate q² — recurrent states that look like shared ancestry but are not.
// d is capped at 0.5, the maximum independent per-site draws can produce.
//
// "exclusive". Each variable site is carried by EXACTLY ONE allele; every other allele is
// ancestral there. So no two created alleles ever share a derived state, and any sharing
// observed downstream is genuine shared inheritance rather than coincidence. A site is
// variable with probability p and its variant is assigned to a uniformly-chosen allele:
//
//	p = A·d/2   ⇒   each allele holds d·L/2 private sites   ⇒   pairwise difference = d
//
// Because p is a probability, this model caps d at 2/A — tighter than 0.5 whenever A > 4.
// The A alleles form a PERFECT star: only private variants, no shared derived states, and
// therefore no internal structure whatsoever for an inference method to read as a
// genealogy.
//
// RNG ORDER differs between the two and is documented here because reproducibility depends
// on it: "independent" is allele-major then bit-minor; "exclusive" is bit-major, drawing
// one uniform per site plus one allele index per variable site.
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
	siteModel := createdSiteModel(model)

	d := model.Parameters["founder_allele_divergence"]
	if d < 0 {
		d = 0
	}
	// Each site model has its own ceiling on the pairwise divergence it can express.
	maxD := 0.5
	if siteModel == siteExclusive {
		maxD = 2 / float64(alleles)
	}
	if d > maxD {
		fmt.Printf("   founder_allele_divergence %.4f capped at %.4g (the maximum pairwise "+
			"difference the %q site model can produce with A=%d)\n", d, maxD, siteModel, alleles)
		d = maxD
	}

	pop.CreatedAlleleSeqs = make([][]uint64, alleles)
	for a := 0; a < alleles; a++ {
		pop.CreatedAlleleSeqs[a] = make([]uint64, words)
	}

	q := 0.0
	if siteModel == siteExclusive {
		// Bit-major: one draw decides whether the site varies at all, a second picks the
		// single allele that carries it. No site is ever derived in two alleles, so the
		// created alleles form a perfect star.
		p := float64(alleles) * d / 2
		for bit := 0; bit < totalBits && p > 0; bit++ {
			if utils.RandFloat64() < p {
				a := utils.RandIntn(alleles)
				pop.CreatedAlleleSeqs[a][bit/64] |= 1 << (uint(bit) % 64)
			}
		}
	} else {
		// Allele-major, bit-minor \u2014 the original model, byte-identical to pre-existing runs.
		q = (1 - math.Sqrt(1-2*d)) / 2
		for a := 0; a < alleles; a++ {
			seq := pop.CreatedAlleleSeqs[a]
			if q > 0 {
				for bit := 0; bit < totalBits; bit++ {
					if utils.RandFloat64() < q {
						seq[bit/64] |= 1 << (uint(bit) % 64)
					}
				}
			}
		}
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

	pool := (alleles + len(ids) - 1) / max2(len(ids), 1)
	if siteModel == siteExclusive {
		fmt.Printf("   Created %d founder alleles over %d founders (pairwise divergence d=%.4g, "+
			"site model %q: %.4g of sites variable, each carried by exactly ONE allele, no "+
			"recurrent states); germline pools of up to %d alleles per founder\n",
			alleles, len(ids), d, siteModel, float64(alleles)*d/2, pool)
	} else {
		fmt.Printf("   Created %d founder alleles over %d founders (pairwise divergence d=%.4g, "+
			"site model %q, per-site q=%.4g); germline pools of up to %d alleles per founder\n",
			alleles, len(ids), d, siteModel, q, pool)
	}
}

// Site-model names for founder_allele_site_model.
const (
	siteIndependent = "independent"
	siteExclusive   = "exclusive"
)

// createdSiteModel resolves founder_allele_site_model, defaulting to the original
// "independent" draw so an existing model keeps its behaviour untouched.
func createdSiteModel(model *core.Model) string {
	switch v := strings.ToLower(strings.TrimSpace(
		model.StringParam("founder_allele_site_model", siteIndependent))); v {
	case "", siteIndependent:
		return siteIndependent
	case siteExclusive:
		return siteExclusive
	default:
		fmt.Printf("WARNING: founder_allele_site_model %q is not recognised, using %q "+
			"(want %q or %q)\n", v, siteIndependent, siteIndependent, siteExclusive)
		return siteIndependent
	}
}

func max2(a, b int) int {
	if a > b {
		return a
	}
	return b
}
