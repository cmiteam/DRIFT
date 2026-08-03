package analysis

import (
	"drift/pkg/core"
	"os"
	"testing"
)

// A 16-bit genome, one chromosome, two 8-bit arms: arm0 = [0,8), arm1 = [8,16).
// Bitfields fit in a single uint64 word, so the fixtures are exact and readable.

func rohStrand(positions ...int) []uint64 {
	w := make([]uint64, 1)
	for _, p := range positions {
		w[p/64] |= 1 << uint(p%64)
	}
	return w
}

// rohModel builds a minimal model with the 2-arm genome and the given ROH config.
// armCM (may be nil) supplies a real per-arm genetic map; with nil the module uses
// the uniform recomb_cM_per_bit fallback (set to 1 here ⇒ cM == bits).
func rohModel(minBits int, shortMax, longMin float64, armCM map[int]map[int]float64) *core.Model {
	return &core.Model{
		Parameters: map[string]float64{
			"roh_min_bits":      float64(minBits),
			"roh_short_max_cM":  shortMax,
			"roh_long_min_cM":   longMin,
			"recomb_cM_per_bit": 1,
		},
		FreeParameters: map[string]int{"genome_bits": 16, "run": 1, "year": 100},
		ChromosomeArms: map[int]map[int][]int{0: {0: {0, 8}, 1: {8, 8}}},
		ArmCM:          armCM,
		ModelName:      "ROHFixture",
	}
}

func rohPop() *core.Pop {
	return &core.Pop{Chromosomes: map[int][][]uint64{}}
}

// A fully-homozygous individual: both strands all-zero ⇒ homozygous at all 16 bits.
// Per-arm bounding must split this into TWO 8-bit runs (one per arm), NOT one 16-bit
// run — the run may not cross the centromere/arm boundary.
func TestComputeROH_PerArmBoundingAndFullHomozygosity(t *testing.T) {
	m := rohModel(5, 5, 10, nil)
	pop := rohPop()
	pop.Chromosomes[1] = [][]uint64{rohStrand(), rohStrand()}

	res := ComputeROH(m, pop, []int{1})
	if res.NumSampled != 1 {
		t.Fatalf("NumSampled = %d, want 1", res.NumSampled)
	}
	ind := res.Individuals[0]
	if ind.NumRuns != 2 {
		t.Fatalf("NumRuns = %d, want 2 (per-arm bounding splits the 16-bit homozygous span)", ind.NumRuns)
	}
	if ind.TotalBits != 16 || ind.HetSites != 0 {
		t.Fatalf("TotalBits=%d HetSites=%d, want 16 and 0", ind.TotalBits, ind.HetSites)
	}
	if ind.FROHBits != 1.0 {
		t.Fatalf("FROHBits = %v, want 1.0", ind.FROHBits)
	}
	// Each run is 8 bits ⇒ 8 cM at 1 cM/bit ⇒ medium class (5 ≤ 8 < 10).
	if ind.ClassCount[ROHMedium] != 2 || ind.ClassCount[ROHShort] != 0 || ind.ClassCount[ROHLong] != 0 {
		t.Fatalf("class counts = %v, want [0 2 0]", ind.ClassCount)
	}
	if ind.TotalCM != 16 || ind.FROHCM != 1.0 {
		t.Fatalf("TotalCM=%v FROHCM=%v, want 16 and 1.0", ind.TotalCM, ind.FROHCM)
	}
}

// An outbred individual: heterozygous at every odd bit ⇒ homozygous bits are all
// isolated (length-1 runs), so with roh_min_bits=5 NOTHING is called an ROH.
func TestComputeROH_OutbredNoRuns(t *testing.T) {
	m := rohModel(5, 5, 10, nil)
	pop := rohPop()
	pop.Chromosomes[1] = [][]uint64{rohStrand(), rohStrand(1, 3, 5, 7, 9, 11, 13, 15)}

	res := ComputeROH(m, pop, []int{1})
	ind := res.Individuals[0]
	if ind.HetSites != 8 {
		t.Fatalf("HetSites = %d, want 8", ind.HetSites)
	}
	if ind.NumRuns != 0 || ind.TotalBits != 0 || ind.FROHBits != 0 {
		t.Fatalf("outbred individual should have no called ROH; got NumRuns=%d TotalBits=%d FROHBits=%v",
			ind.NumRuns, ind.TotalBits, ind.FROHBits)
	}
}

// Distinguishable outcome (Rob's requirement): a consanguineous (fully homozygous)
// individual shows more and longer ROH and higher F_ROH than an outbred control,
// through the real ComputeROH, and the population aggregates reflect it.
func TestComputeROH_ConsanguineousVsOutbred(t *testing.T) {
	m := rohModel(5, 5, 10, nil)
	pop := rohPop()
	pop.Chromosomes[1] = [][]uint64{rohStrand(), rohStrand()}                          // consanguineous
	pop.Chromosomes[2] = [][]uint64{rohStrand(), rohStrand(1, 3, 5, 7, 9, 11, 13, 15)} // outbred

	res := ComputeROH(m, pop, []int{1, 2})
	if res.NumSampled != 2 {
		t.Fatalf("NumSampled = %d, want 2", res.NumSampled)
	}
	consang, outbred := res.Individuals[0], res.Individuals[1]
	if !(consang.FROHBits > outbred.FROHBits) {
		t.Fatalf("consanguineous F_ROH (%v) must exceed outbred (%v)", consang.FROHBits, outbred.FROHBits)
	}
	if !(consang.NumRuns > outbred.NumRuns) {
		t.Fatalf("consanguineous N_ROH (%d) must exceed outbred (%d)", consang.NumRuns, outbred.NumRuns)
	}
	// Mean F_ROH is the average of 1.0 and 0.0.
	if res.MeanFROHBits != 0.5 {
		t.Fatalf("MeanFROHBits = %v, want 0.5", res.MeanFROHBits)
	}
	if res.MeanHetFrac != 0.25 { // (0 + 8/16) / 2
		t.Fatalf("MeanHetFrac = %v, want 0.25", res.MeanHetFrac)
	}
}

// The roh_min_bits native floor: a run shorter than the floor is not called.
func TestComputeROH_MinBitsFloor(t *testing.T) {
	pop := rohPop()
	// arm0: homozygous [0,4), het at bit 4, homozygous [5,8) ⇒ runs of 4 and 3 bits.
	// arm1: strand0=0 / strand1=1 across 8..15 ⇒ all het ⇒ no runs.
	pop.Chromosomes[1] = [][]uint64{
		rohStrand(4),                            // strand0: 1 only at bit 4
		rohStrand(8, 9, 10, 11, 12, 13, 14, 15), // strand1: 1s across arm1
	}
	mHigh := rohModel(5, 5, 10, nil) // floor 5 ⇒ neither 4-bit nor 3-bit run counts
	if got := ComputeROH(mHigh, pop, []int{1}).Individuals[0].NumRuns; got != 0 {
		t.Fatalf("min_bits=5: NumRuns = %d, want 0", got)
	}
	mLow := rohModel(3, 5, 10, nil) // floor 3 ⇒ both the 4-bit and 3-bit runs count
	indLow := ComputeROH(mLow, pop, []int{1}).Individuals[0]
	if indLow.NumRuns != 2 || indLow.TotalBits != 7 {
		t.Fatalf("min_bits=3: NumRuns=%d TotalBits=%d, want 2 and 7", indLow.NumRuns, indLow.TotalBits)
	}
	if indLow.HetSites != 9 { // 1 in arm0 + 8 in arm1
		t.Fatalf("HetSites = %d, want 9", indLow.HetSites)
	}
}

// cM via a real per-arm genetic map: run cM is read from the map, and the class
// partition follows cM, not bits. arm0 = 16 cM over 8 bits, arm1 = 4 cM over 8 bits.
func TestComputeROH_GeneticMapCMAndClasses(t *testing.T) {
	armCM := map[int]map[int]float64{0: {0: 16.0, 1: 4.0}}
	m := rohModel(5, 5, 10, armCM)
	pop := rohPop()
	pop.Chromosomes[1] = [][]uint64{rohStrand(), rohStrand()} // fully homozygous

	res := ComputeROH(m, pop, []int{1})
	if !res.HasMap {
		t.Fatalf("HasMap = false, want true (ArmCM supplied)")
	}
	ind := res.Individuals[0]
	// arm0 run = 16 cM ⇒ long (>=10); arm1 run = 4 cM ⇒ short (<5). No medium.
	if ind.ClassCount[ROHLong] != 1 || ind.ClassCount[ROHShort] != 1 || ind.ClassCount[ROHMedium] != 0 {
		t.Fatalf("class counts = %v, want [1 0 1] (short/medium/long)", ind.ClassCount)
	}
	if ind.TotalCM != 20 { // 16 + 4
		t.Fatalf("TotalCM = %v, want 20", ind.TotalCM)
	}
	// Total map cM = 20 ⇒ F_ROH_cM = 1.0 (whole genome homozygous).
	if ind.FROHCM != 1.0 {
		t.Fatalf("FROHCM = %v, want 1.0", ind.FROHCM)
	}
}

// SaveROH writes both CSVs and is a no-op with a clean skip when nothing is eligible.
func TestSaveROH_FilesAndEmptyNoOp(t *testing.T) {
	dir := t.TempDir()
	m := rohModel(5, 5, 10, nil)
	m.ResultsDir = dir

	// Empty sample ⇒ no files, nil error.
	if err := SaveROH(m, ComputeROH(m, rohPop(), nil)); err != nil {
		t.Fatalf("SaveROH empty: %v", err)
	}
	if entries, _ := os.ReadDir(dir); len(entries) != 0 {
		t.Fatalf("empty ROH pass wrote %d files, want 0", len(entries))
	}

	// Populated ⇒ two CSVs.
	pop := rohPop()
	pop.Chromosomes[1] = [][]uint64{rohStrand(), rohStrand()}
	if err := SaveROH(m, ComputeROH(m, pop, []int{1})); err != nil {
		t.Fatalf("SaveROH: %v", err)
	}
	for _, kind := range []string{"per_individual", "distribution"} {
		p := rohPath(m, kind)
		if _, err := os.Stat(p); err != nil {
			t.Fatalf("expected %s to exist: %v", p, err)
		}
	}
}
