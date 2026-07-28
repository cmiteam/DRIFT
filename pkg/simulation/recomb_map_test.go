package simulation

import (
	"bytes"
	"drift/pkg/core"
	"drift/pkg/utils"
	"sort"
	"testing"
)

// recombModel builds a model over the given arms with the requested recombination
// model and an optional per-arm cM table.
func recombModel(recomb string, arms map[int]map[int][]int, armCM map[int]map[int]float64) *core.Model {
	genomeBits := 0
	for _, ca := range arms {
		for _, a := range ca {
			if end := a[0] + a[1]; end > genomeBits {
				genomeBits = end
			}
		}
	}
	return &core.Model{
		Parameters:     map[string]float64{},
		StringParams:   map[string]string{"recombination_model": recomb},
		FreeParameters: map[string]int{"genome_bits": genomeBits},
		ChromosomeArms: arms,
		ArmCM:          armCM,
	}
}

// legacyCreateMaskOracle is a verbatim copy of the ORIGINAL single-interior-segment
// createMask body (pre-§6a), used to prove the legacy path is byte-identical: same
// mask AND same RNG consumption. If the production legacy path ever drifts from this,
// TestRecombLegacyByteIdentical fails.
func legacyCreateMaskOracle(model *core.Model, sex int) ([]uint64, []uint64) {
	genomeBits := model.FreeParameters["genome_bits"]
	if genomeBits <= 0 {
		return []uint64{}, []uint64{0}
	}
	genomeArrSize := (model.FreeParameters["genome_bits"] + 63) / 64
	genomemask := make([]uint64, genomeArrSize)
	centromask := []uint64{0}

	chroms := make([]int, 0, len(model.ChromosomeArms))
	for chrom := range model.ChromosomeArms {
		chroms = append(chroms, chrom)
	}
	sort.Ints(chroms)

	for _, chrom := range chroms {
		arm0, ok0 := model.ChromosomeArms[chrom][0]
		arm1, ok1 := model.ChromosomeArms[chrom][1]
		if !ok0 || !ok1 {
			continue
		}
		pstart := arm0[0]
		plen := arm0[1]
		qstart := arm1[0]
		qlen := arm1[1]
		if plen <= 0 || qlen <= 0 {
			continue
		}
		ploc := utils.RandIntn(plen)
		qloc := utils.RandIntn(qlen)
		whichCopy := utils.RandIntn(2)
		if whichCopy == 1 {
			startBit := pstart + ploc
			endBit := qstart + qloc
			for i := startBit; i < endBit; i++ {
				genomemask[i/64] |= (1 << (i % 64))
			}
		}
	}
	return genomemask, centromask
}

// TestRecombLegacyByteIdentical is the compatibility guard: with recombination_model
// unset/"legacy" the production createMask must produce the identical mask and consume
// the identical RNG stream as the pre-§6a code. This is what keeps drift -validate PASS.
func TestRecombLegacyByteIdentical(t *testing.T) {
	// Realistic multi-chromosome layout (Default-like tiling).
	arms := map[int]map[int][]int{
		1: {0: {0, 125}, 1: {125, 126}},
		2: {0: {251, 93}, 1: {344, 153}},
		3: {0: {497, 91}, 1: {588, 112}},
	}
	for _, recomb := range []string{"legacy", ""} { // "" => StringParam default "legacy"
		model := recombModel(recomb, arms, nil)

		utils.SeedRNG(99)
		var prodMasks [][]uint64
		for i := 0; i < 200; i++ {
			m, _ := createMask(model, i%2)
			prodMasks = append(prodMasks, m)
		}
		prodState := utils.RNGState()

		// Fresh model (avoid any cached map) + same seed for the oracle.
		oracleModel := recombModel(recomb, arms, nil)
		utils.SeedRNG(99)
		var oracleMasks [][]uint64
		for i := 0; i < 200; i++ {
			m, _ := legacyCreateMaskOracle(oracleModel, i%2)
			oracleMasks = append(oracleMasks, m)
		}
		oracleState := utils.RNGState()

		if !bytes.Equal(prodState, oracleState) {
			t.Fatalf("recomb=%q: RNG state diverged from legacy oracle (stream not byte-identical)", recomb)
		}
		for i := range prodMasks {
			if len(prodMasks[i]) != len(oracleMasks[i]) {
				t.Fatalf("recomb=%q trial %d: mask length %d != %d", recomb, i, len(prodMasks[i]), len(oracleMasks[i]))
			}
			for w := range prodMasks[i] {
				if prodMasks[i][w] != oracleMasks[i][w] {
					t.Fatalf("recomb=%q trial %d word %d: mask %#x != %#x", recomb, i, w, prodMasks[i][w], oracleMasks[i][w])
				}
			}
		}
	}
}

// maskBit reads genome bit p from a mask word array (LSB = position 0), matching the
// meiosis / createMask convention.
func maskBit(mask []uint64, p int) int { return int((mask[p/64] >> (uint(p) % 64)) & 1) }

// countCrossovers counts transitions (recombination breakpoints) within [lo,hi).
func countCrossovers(mask []uint64, lo, hi int) int {
	x := 0
	prev := maskBit(mask, lo)
	for p := lo + 1; p < hi; p++ {
		b := maskBit(mask, p)
		if b != prev {
			x++
		}
		prev = b
	}
	return x
}

// TestRecombMapAssortmentAndCount exercises the two structural fixes §6a delivers over
// the legacy kernel: (F3) independent assortment — a chromosome starts on either
// homolog 50/50 — and (F1/F2) a variable crossover count that is NOT locked to {0,2}
// bracketing the centromere.
func TestRecombMapAssortmentAndCount(t *testing.T) {
	// One chromosome, high genetic length so extra crossovers are common.
	arms := map[int]map[int][]int{1: {0: {0, 100}, 1: {100, 100}}}
	armCM := map[int]map[int]float64{1: {0: 150, 1: 150}} // 300 cM total
	model := recombModel("map", arms, armCM)

	utils.SeedRNG(2024)
	const trials = 20000
	startCopy0 := 0
	xoverHist := map[int]int{}
	sawOdd, sawGT2 := false, false

	for i := 0; i < trials; i++ {
		m, _ := createMask(model, 0)
		if maskBit(m, 0) == 1 { // chromosome start on copy 0
			startCopy0++
		}
		x := countCrossovers(m, 0, 200)
		xoverHist[x]++
		if x%2 == 1 {
			sawOdd = true
		}
		if x > 2 {
			sawGT2 = true
		}
	}

	frac := float64(startCopy0) / trials
	if frac < 0.45 || frac > 0.55 {
		t.Errorf("start-on-copy0 fraction = %.3f, want ~0.5 (independent assortment)", frac)
	}
	// Legacy would produce only x in {0,1(deg),2}; map mode must show odd counts and >2.
	if !sawOdd {
		t.Error("expected some odd crossover counts (impossible in the legacy kernel)")
	}
	if !sawGT2 {
		t.Error("expected some crossover counts > 2 (impossible in the legacy kernel)")
	}
	// Obligate crossover: with 300 cM, the all-zero (x=0, no recombination) outcome
	// must be rare relative to legacy's ~50%.
	if f0 := float64(xoverHist[0]) / trials; f0 > 0.10 {
		t.Errorf("x=0 fraction = %.3f, want small (obligate crossover + high cM)", f0)
	}
}

// TestRecombMapDistinguishableLinkage is the distinguishable-outcome guard: for two
// sites the SAME bit distance apart, a low cM/bit chromosome keeps them tightly linked
// (they almost always come from the same homolog) while a high cM/bit chromosome
// recombines between them toward independence. The recombination fraction between two
// sites is set by the cM between them, so the identical 10-bit gap yields very
// different linkage under the two maps. Same engine, same seed — only the map differs.
func TestRecombMapDistinguishableLinkage(t *testing.T) {
	// Identical bit geometry; chr1 low cM/bit, chr2 high cM/bit.
	arms := map[int]map[int][]int{
		1: {0: {0, 100}, 1: {100, 100}},
		2: {0: {200, 100}, 1: {300, 100}},
	}
	armCM := map[int]map[int]float64{
		1: {0: 1, 1: 1},     // 0.01 cM/bit -> a 10-bit gap = 0.1 cM (tight)
		2: {0: 400, 1: 400}, // 4 cM/bit    -> a 10-bit gap = 40 cM  (loose)
	}
	model := recombModel("map", arms, armCM)

	utils.SeedRNG(555)
	const trials = 20000
	// Two close sites (10 bits apart) within the p arm of each chromosome.
	agreeLow, agreeHigh := 0, 0
	for i := 0; i < trials; i++ {
		m, _ := createMask(model, 0)
		if maskBit(m, 10) == maskBit(m, 20) { // chr1, 0.1 cM apart
			agreeLow++
		}
		if maskBit(m, 210) == maskBit(m, 220) { // chr2, 40 cM apart
			agreeHigh++
		}
	}
	low := float64(agreeLow) / trials
	high := float64(agreeHigh) / trials
	if low < 0.95 {
		t.Errorf("low cM/bit site agreement = %.3f, want > 0.95 (tight linkage over 0.1 cM)", low)
	}
	if high > 0.80 {
		t.Errorf("high cM/bit site agreement = %.3f, want lower (recombining over 40 cM)", high)
	}
	if !(low > high+0.15) {
		t.Errorf("expected low (%.3f) clearly > high (%.3f) — the map must change linkage", low, high)
	}
}
