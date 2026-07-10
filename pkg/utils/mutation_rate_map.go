package utils

import (
	"sort"
	"strconv"
	"strings"
)

// RateMap makes the per-position mutation rate non-uniform across the genome
// (roadmap §1, "variable mutation rate"). DRIFT historically drew every new
// mutation's position from a flat RandIntn(genome_bits); real genomes have
// hotspots and coldspots (CpG sites, recombination-associated hotspots,
// mutation-cold housekeeping regions, ...). A RateMap partitions [0, genome_bits)
// into segments each carrying a rate multiplier (1 = baseline outside every
// named region) and draws positions in proportion to segment weight = multiplier
// * length, so mutations concentrate where the multiplier is high.
//
// STRICT NO-OP: an empty spec parses to a nil *RateMap, and GenerateNewMutations
// then keeps the literal RandIntn(genome_bits) draw — the identical single RNG
// draw as before this feature existed. Only runs that configure a map draw the
// two-step weighted stream (RandFloat64 to pick a segment, RandIntn within it),
// so the §6h neutral-validation harness stays byte-identical.
type RateMap struct {
	segStart  []int     // inclusive start bit of each segment
	segLen    []int     // length in bits of each segment (>0)
	cumWeight []float64 // cumulative segment weights; cumWeight[i] = Σ weight[0..i]
	total     float64   // total weight (== cumWeight[len-1]); 0 only if degenerate
}

// rateRegion is one parsed START-END:MULT region (END exclusive), pre-clamp.
type rateRegion struct {
	start, end int
	mult       float64
}

// parseRateMap parses a rate-map spec over a genome of genomeBits positions and
// returns the sampler, or nil when the spec yields no usable non-uniform map (so
// the caller falls back to the strict-no-op uniform draw).
//
// Regions are separated by ';' or ',' (so the SAME parser serves both the global
// `mutation_rate_map` param, which uses ';', and a per-class inline `ratemap=`
// override, which must use ',' to avoid the class spec's ';' separator). Each
// region is START-END:MULT with no interior whitespace — START/END are genome-bit
// positions (END exclusive) and MULT is the rate multiplier for that span:
//
//	1000-2000:10 ; 50000-51000:0.1
//
// A 10x hotspot at bits 1000..2000 and a 0.1x coldspot at 50000..51000; every
// other position keeps the baseline multiplier 1. Regions are clamped to
// [0, genomeBits); later regions override earlier ones on overlap. Returns nil if
// the spec is empty, has no parseable region, or leaves the genome uniform
// (multiplier 1 everywhere) — all of which are equivalent to the uniform no-op.
func parseRateMap(spec string, genomeBits int) *RateMap {
	spec = strings.TrimSpace(spec)
	if spec == "" || genomeBits <= 0 {
		return nil
	}

	regions := make([]rateRegion, 0, 4)
	tokens := strings.FieldsFunc(spec, func(r rune) bool { return r == ';' || r == ',' })
	for _, tok := range tokens {
		tok = strings.TrimSpace(tok)
		if tok == "" {
			continue
		}
		colon := strings.LastIndex(tok, ":")
		if colon < 0 {
			continue
		}
		rng := strings.TrimSpace(tok[:colon])
		multStr := strings.TrimSpace(tok[colon+1:])
		dash := strings.Index(rng, "-")
		if dash < 0 {
			continue
		}
		start, err1 := strconv.Atoi(strings.TrimSpace(rng[:dash]))
		end, err2 := strconv.Atoi(strings.TrimSpace(rng[dash+1:]))
		mult, err3 := strconv.ParseFloat(multStr, 64)
		if err1 != nil || err2 != nil || err3 != nil {
			continue
		}
		// Clamp to the genome; skip empty, backwards, or negative-multiplier spans.
		if start < 0 {
			start = 0
		}
		if end > genomeBits {
			end = genomeBits
		}
		if end <= start || mult < 0 {
			continue
		}
		regions = append(regions, rateRegion{start: start, end: end, mult: mult})
	}
	if len(regions) == 0 {
		return nil
	}

	// Build elementary segments over [0, genomeBits) at every region boundary; each
	// segment's multiplier is that of the LAST spec region covering it, or 1 if
	// none (baseline). This handles sparse, adjacent, and overlapping regions.
	bounds := map[int]struct{}{0: {}, genomeBits: {}}
	for _, r := range regions {
		bounds[r.start] = struct{}{}
		bounds[r.end] = struct{}{}
	}
	cuts := make([]int, 0, len(bounds))
	for b := range bounds {
		cuts = append(cuts, b)
	}
	sort.Ints(cuts)

	rm := &RateMap{}
	nonBaseline := false
	for i := 0; i+1 < len(cuts); i++ {
		a, b := cuts[i], cuts[i+1]
		if b <= a {
			continue
		}
		mult := 1.0
		for _, r := range regions { // later regions override earlier on overlap
			if r.start <= a && b <= r.end {
				mult = r.mult
			}
		}
		if mult != 1.0 {
			nonBaseline = true
		}
		w := mult * float64(b-a)
		if w <= 0 {
			continue // zero-multiplier (dead) span: never drawn
		}
		rm.total += w
		rm.segStart = append(rm.segStart, a)
		rm.segLen = append(rm.segLen, b-a)
		rm.cumWeight = append(rm.cumWeight, rm.total)
	}

	// A map that leaves every position at the baseline multiplier (e.g. only mult=1
	// regions) is indistinguishable from uniform: return nil so the no-op draw is
	// used and the RNG stream is unchanged. Also guard the degenerate all-dead case.
	if !nonBaseline || rm.total <= 0 || len(rm.segStart) == 0 {
		return nil
	}
	return rm
}

// Draw returns a genome position sampled in proportion to the per-segment rate
// weight. It consumes RandFloat64 (segment) + RandIntn (position within segment) —
// deliberately NOT the single RandIntn(genome_bits) of the uniform path, so it is
// only ever reached for explicitly configured (non-default) rate maps.
func (rm *RateMap) Draw() int {
	r := RandFloat64() * rm.total
	i := sort.Search(len(rm.cumWeight), func(k int) bool { return rm.cumWeight[k] > r })
	if i >= len(rm.segStart) {
		i = len(rm.segStart) - 1
	}
	return rm.segStart[i] + RandIntn(rm.segLen[i])
}
