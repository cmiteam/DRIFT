package analysis

// LD-based recent-Ne companion (roadmap §6d, "GONE-style" leg). Under
// drift–recombination balance the expected squared correlation between two loci at
// recombination fraction c is
//     E[r^2] ≈ 1/(4*Ne*c) + 1/S      (Sved 1971; Waples 2006 sampling term)
// so a sample's r^2 at a chosen c reports Ne about 1/(2c) generations ago — loosely
// linked pairs (c ≈ 0.5) see the RECENT Ne, which is exactly the regime creationist
// bottleneck questions care about and the regime SFS methods are blind to. Solving
// for Ne over a band of loosely-linked pairs:
//     Ne ≈ 1 / ( 4 * cBar * ( rBar^2 - 1/S ) ).
//
// WHY THIS IS DELIBERATELY APPROXIMATE. Turning a bit-distance into a recombination
// fraction needs a genetic map, which DRIFT does not have yet (real recombination
// maps are §6a; the current createMask model is a placeholder). So the conversion is
// a single explicit knob, ne_cM_per_bit (cM per genome bit), and the LD leg is OFF
// by default (ne_cM_per_bit = 0). Loci on different chromosomes are treated as
// unlinked (c = 0.5); on the same chromosome, c follows Haldane's map function from
// the bit distance. The rigorous, calibration-free estimate is the temporal one in
// ne_temporal.go; this leg exists so the GONE-style number is available and directly
// comparable once §6a lands, and it is labelled approximate everywhere it surfaces.
//
// No RNG: segregating sites and the pair walk are enumerated deterministically (a
// fixed stride, not random sampling), and CalculatePairwiseLD only counts.

import (
	"drift/pkg/core"
	"drift/pkg/utils"
	"math"
)

// recentCBand is the lower bound on recombination fraction for a pair to count
// toward the "recent" Ne estimate. c >= 0.2 corresponds to roughly <= 2.5
// generations ago (time ≈ 1/(2c)), i.e. the recent history GONE targets.
const recentCBand = 0.2

// ldRecentNe estimates the recent effective size from LD (r^2) among loosely-linked
// site pairs. Returns (0, 0) when disabled (ne_cM_per_bit <= 0), when there are too
// few segregating sites, or when the estimate is unresolved (sampling noise exceeds
// the LD signal). pooled is the pooled member set (used only to require a non-empty
// sample); r^2 itself is computed over the full tracked population, matching
// CalculatePairwiseLD.
func ldRecentNe(model *core.Model, pop *core.Pop, pooled []int) (float64, int) {
	cMPerBit := model.Parameters["ne_cM_per_bit"]
	if cMPerBit <= 0 || len(pooled) < 2 {
		return 0, 0 // LD leg disabled or nothing to work with
	}

	sites := getSegregatingSites(model, pop)
	if len(sites) < 2 {
		return 0, 0
	}

	// Deterministic pair subsampling: enumerate all unordered pairs in a fixed
	// order and take every stride-th one so the work is bounded by ~ld_max_pairs
	// without drawing any RNG.
	maxPairs := int(model.Parameters["ld_max_pairs"])
	if maxPairs <= 0 {
		maxPairs = 10000
	}
	total := len(sites) * (len(sites) - 1) / 2
	stride := 1
	if total > maxPairs {
		stride = (total + maxPairs - 1) / maxPairs
	}

	// Haploid sample size for the 1/S sampling correction, over the population
	// CalculatePairwiseLD actually counts (all individuals with chromosomes).
	nChrom := 0
	for id := range pop.Chromosomes {
		c := pop.Chromosomes[id]
		if len(c) >= 2 && c[0] != nil && c[1] != nil {
			nChrom++
		}
	}
	if nChrom < 2 {
		return 0, 0
	}
	invS := 1.0 / float64(2*nChrom)

	var sumR2, sumC float64
	used := 0
	idx := 0
	for i := 0; i < len(sites); i++ {
		for j := i + 1; j < len(sites); j++ {
			take := idx%stride == 0
			idx++
			if !take {
				continue
			}
			c := recombFraction(model, sites[i], sites[j], cMPerBit)
			if c < recentCBand {
				continue // too tightly linked to inform recent Ne
			}
			ld := CalculatePairwiseLD(model, pop, sites[i], sites[j])
			if ld == nil {
				continue
			}
			sumR2 += ld.RSquared
			sumC += c
			used++
		}
	}
	if used == 0 {
		return 0, 0
	}
	rBar := sumR2 / float64(used)
	cBar := sumC / float64(used)
	adj := rBar - invS
	if adj <= 0 || cBar <= 0 {
		return 0, used // LD below sampling noise: unresolved
	}
	return 1.0 / (4.0 * cBar * adj), used
}

// recombFraction maps a pair of genome-bit positions to a recombination fraction.
// Different chromosomes are unlinked (c = 0.5); on the same chromosome the bit
// distance is converted to centiMorgans via cMPerBit and then to a recombination
// fraction with Haldane's map function c = 0.5*(1 - e^{-2M}). Positions off the arm
// map are treated as unlinked.
func recombFraction(model *core.Model, s1, s2 int, cMPerBit float64) float64 {
	chrom1, _, ok1 := utils.PositionArm(model, s1)
	chrom2, _, ok2 := utils.PositionArm(model, s2)
	if !ok1 || !ok2 || chrom1 != chrom2 {
		return 0.5
	}
	dist := s2 - s1
	if dist < 0 {
		dist = -dist
	}
	morgans := float64(dist) * cMPerBit / 100.0
	c := 0.5 * (1 - math.Exp(-2*morgans))
	if c > 0.5 {
		c = 0.5
	}
	return c
}
