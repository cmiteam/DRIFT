package core

import "sort"

// GeneticMap is a piecewise-linear genetic map over the genome bit-coordinate space
// (roadmap §6a, "read real recombination maps and real human chromosome structure").
// It is built from the existing chromosome-arm scaffold: within each arm the
// recombination rate is treated as constant (arm cM / arm bit-length), so cumulative
// cM is piecewise-linear with a knot at each arm boundary. This is the "per-arm cM"
// layer chosen with Rob; a finer bin map (hotspots) can be substituted behind the
// same accessor later without touching the meiosis kernel or the LD/EHH consumers.
//
// The type is deliberately pure data + pure methods (no RNG, no utils import) so it
// can live in core and be held on Model. createMask (pkg/simulation) does the RNG
// draws and calls BitAtCM to place crossovers; recombFraction (pkg/analysis) calls
// CumCM to turn a bit distance into a real map distance. All bit positions are GLOBAL
// genome-bit indices (arm Start values), matching PositionArm and the founder bitfield.
type GeneticMap struct {
	chroms map[int]*gmChrom
}

type gmChrom struct {
	segs     []gmSeg // ordered by startBit (p arm then q arm)
	totalCM  float64
	startBit int // first bit of the chromosome (min arm start)
	endBit   int // one past the last bit (max arm end)
}

type gmSeg struct {
	startBit int
	length   int     // bits in this arm segment
	cM       float64 // genetic length of this segment
	cumCM    float64 // cumulative cM at startBit, from the chromosome's start
}

// BuildGeneticMap constructs a GeneticMap from the chromosome-arm ranges and an
// optional per-arm cM table. armCM[chrom][arm] gives an arm's genetic length in
// centiMorgans; when an arm is absent from armCM its cM defaults to
// length*defaultCMPerBit (a uniform fallback rate), so a 4-column chromosome_data.csv
// with no cM column still yields a usable uniform map. Arms are laid out in arm-index
// order (0 = p, 1 = q) so the map is contiguous across the centromere.
func BuildGeneticMap(arms map[int]map[int][]int, armCM map[int]map[int]float64, defaultCMPerBit float64) *GeneticMap {
	g := &GeneticMap{chroms: make(map[int]*gmChrom)}
	for chrom, chromArms := range arms {
		gc := &gmChrom{startBit: -1}
		cum := 0.0
		for arm := 0; arm < 2; arm++ {
			armData, ok := chromArms[arm]
			if !ok || len(armData) < 2 {
				continue
			}
			start, length := armData[0], armData[1]
			if length <= 0 {
				continue
			}
			cm := float64(length) * defaultCMPerBit
			if armCM != nil {
				if v, ok := armCM[chrom][arm]; ok {
					cm = v
				}
			}
			gc.segs = append(gc.segs, gmSeg{startBit: start, length: length, cM: cm, cumCM: cum})
			cum += cm
			if gc.startBit < 0 || start < gc.startBit {
				gc.startBit = start
			}
			if end := start + length; end > gc.endBit {
				gc.endBit = end
			}
		}
		if len(gc.segs) == 0 {
			continue
		}
		sort.Slice(gc.segs, func(i, j int) bool { return gc.segs[i].startBit < gc.segs[j].startBit })
		gc.totalCM = cum
		g.chroms[chrom] = gc
	}
	return g
}

// Has reports whether the map has a (non-degenerate) entry for the chromosome.
func (g *GeneticMap) Has(chrom int) bool {
	_, ok := g.chroms[chrom]
	return ok
}

// TotalCM returns the chromosome's total genetic length in centiMorgans.
func (g *GeneticMap) TotalCM(chrom int) float64 {
	if gc, ok := g.chroms[chrom]; ok {
		return gc.totalCM
	}
	return 0
}

// Span returns the chromosome's global bit range [start, end).
func (g *GeneticMap) Span(chrom int) (int, int) {
	if gc, ok := g.chroms[chrom]; ok {
		return gc.startBit, gc.endBit
	}
	return 0, 0
}

// CumCM returns the cumulative genetic distance in centiMorgans from the
// chromosome's start to global bit position pos, interpolating linearly within the
// containing arm segment. Positions before the first segment clamp to 0; positions
// past the last clamp to totalCM.
func (g *GeneticMap) CumCM(chrom, pos int) float64 {
	gc, ok := g.chroms[chrom]
	if !ok {
		return 0
	}
	for _, s := range gc.segs {
		if pos < s.startBit {
			return s.cumCM // in a gap before this segment
		}
		if pos < s.startBit+s.length {
			if s.length == 0 {
				return s.cumCM
			}
			return s.cumCM + float64(pos-s.startBit)*s.cM/float64(s.length)
		}
	}
	return gc.totalCM
}

// BitAtCM inverts CumCM: given a cM offset in [0, totalCM] from the chromosome's
// start, it returns the global bit position at that genetic distance. Used to place a
// crossover drawn uniformly in genetic (map) space, so crossovers concentrate where
// cM/bit is high (the hotspot behavior a real map encodes). Zero-cM segments are
// skipped. The result is clamped to the chromosome's bit span.
func (g *GeneticMap) BitAtCM(chrom int, cm float64) int {
	gc, ok := g.chroms[chrom]
	if !ok {
		return 0
	}
	if cm < 0 {
		cm = 0
	}
	for _, s := range gc.segs {
		if s.cM <= 0 {
			continue
		}
		if cm <= s.cumCM+s.cM {
			off := int((cm - s.cumCM) / s.cM * float64(s.length))
			if off < 0 {
				off = 0
			}
			if off >= s.length {
				off = s.length - 1
			}
			return s.startBit + off
		}
	}
	// cm beyond the last positive segment: last bit of the last segment.
	last := gc.segs[len(gc.segs)-1]
	return last.startBit + last.length - 1
}
