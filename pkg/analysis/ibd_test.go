package analysis

import (
	"encoding/csv"
	"os"
	"path/filepath"
	"testing"

	"drift/pkg/core"
)

// makeIBDTestModel builds a minimal model with one chromosome of 128 bits:
// p arm = bits [0,64) (word 0), q arm = bits [64,128) (word 1).
func makeIBDTestModel() *core.Model {
	return &core.Model{
		ChromosomeArms: map[int]map[int][]int{
			0: {
				0: {0, 64},  // p arm: start 0, length 64
				1: {64, 64}, // q arm: start 64, length 64
			},
		},
		FreeParameters: map[string]int{"genome_bits": 128, "run": 1, "year": 10},
		Parameters:     map[string]float64{},
		ModelName:      "IBDTest",
	}
}

// chrom builds a [2][]uint64 chromosome (two strand copies) from raw words.
func chrom(copy0Word0, copy0Word1, copy1Word0, copy1Word1 uint64) [][]uint64 {
	return [][]uint64{
		{copy0Word0, copy0Word1},
		{copy1Word0, copy1Word1},
	}
}

const allOnes = ^uint64(0)
const low32 = uint64(0xFFFFFFFF) // bits 0..31 set

func makeIBDTestPop() *core.Pop {
	return &core.Pop{
		IndData: map[int][]int{1: {}, 2: {}, 3: {}, 4: {}},
		Chromosomes: map[int][][]uint64{
			// Ind 1 & 2: entire p arm (64 bits) descends from founder; q arm empty.
			1: chrom(allOnes, 0, 0, 0),
			2: chrom(allOnes, 0, 0, 0),
			// Ind 3: only the low 32 bits of the p arm.
			3: chrom(low32, 0, 0, 0),
			// Ind 4: no chromosomes at all -> must be skipped, not panic.
		},
	}
}

func TestExtractIBD_LongThresholdKeepsOnlyFullArm(t *testing.T) {
	model := makeIBDTestModel()
	pop := makeIBDTestPop()

	// minLen 40 bits: the 64-bit shared arm between 1 and 2 survives;
	// the 32-bit overlaps (1-3, 2-3) are below threshold and dropped.
	res := ExtractIBD(model, pop, []int{1, 2, 3, 4}, 40)

	if len(res.Segments) != 1 {
		t.Fatalf("expected 1 segment, got %d: %+v", len(res.Segments), res.Segments)
	}
	s := res.Segments[0]
	if s.Ind1 != 1 || s.Ind2 != 2 {
		t.Errorf("expected pair (1,2), got (%d,%d)", s.Ind1, s.Ind2)
	}
	if s.LengthBits != 64 || s.StartBit != 0 || s.EndBit != 64 {
		t.Errorf("expected [0,64) len 64, got [%d,%d) len %d", s.StartBit, s.EndBit, s.LengthBits)
	}
	if s.Chrom != 0 || s.Arm != 0 {
		t.Errorf("expected chrom 0 arm 0, got chrom %d arm %d", s.Chrom, s.Arm)
	}
}

func TestExtractIBD_ShortThresholdKeepsPartialOverlaps(t *testing.T) {
	model := makeIBDTestModel()
	pop := makeIBDTestPop()

	// minLen 10 bits: all three overlaps survive -> 64 (1-2), 32 (1-3), 32 (2-3).
	res := ExtractIBD(model, pop, []int{1, 2, 3, 4}, 10)

	if len(res.Segments) != 3 {
		t.Fatalf("expected 3 segments, got %d: %+v", len(res.Segments), res.Segments)
	}

	// Pair (1,2) should total the full 64-bit arm.
	if got := res.TotalByPair[[2]int{1, 2}]; got != 64 {
		t.Errorf("pair (1,2) total = %d, want 64", got)
	}
	// Pairs with ind 3 should each total 32.
	if got := res.TotalByPair[[2]int{1, 3}]; got != 32 {
		t.Errorf("pair (1,3) total = %d, want 32", got)
	}
	if got := res.TotalByPair[[2]int{2, 3}]; got != 32 {
		t.Errorf("pair (2,3) total = %d, want 32", got)
	}

	// Histogram: one 64-bit segment in bin 6, two 32-bit segments in bin 5.
	if res.LengthHist[6] != 1 {
		t.Errorf("expected 1 segment in log2 bin 6 (>=64), got %d", res.LengthHist[6])
	}
	if res.LengthHist[5] != 2 {
		t.Errorf("expected 2 segments in log2 bin 5 (>=32), got %d", res.LengthHist[5])
	}

	// Ind 4 (no chromosomes) should not appear and should not have caused a panic.
	if res.NumPairs != 3 {
		t.Errorf("expected 3 comparable pairs, got %d", res.NumPairs)
	}
}

func TestExtractIBD_RunsAreBoundedWithinArms(t *testing.T) {
	model := makeIBDTestModel()
	// Both arms fully set in both individuals: must yield TWO segments
	// (one per arm), not one segment spanning the arm boundary.
	pop := &core.Pop{
		IndData: map[int][]int{1: {}, 2: {}},
		Chromosomes: map[int][][]uint64{
			1: chrom(allOnes, allOnes, 0, 0),
			2: chrom(allOnes, allOnes, 0, 0),
		},
	}

	res := ExtractIBD(model, pop, []int{1, 2}, 10)
	if len(res.Segments) != 2 {
		t.Fatalf("expected 2 arm-bounded segments, got %d: %+v", len(res.Segments), res.Segments)
	}
	for _, s := range res.Segments {
		if s.LengthBits != 64 {
			t.Errorf("expected each arm segment to be 64 bits, got %d", s.LengthBits)
		}
	}
}

func TestSaveIBD_WritesReadableCSVs(t *testing.T) {
	dir := t.TempDir()
	model := makeIBDTestModel()
	model.ResultsDir = dir
	pop := makeIBDTestPop()

	res := ExtractIBD(model, pop, []int{1, 2, 3}, 10)
	if err := SaveIBD(model, res); err != nil {
		t.Fatalf("SaveIBD failed: %v", err)
	}

	segPath := filepath.Join(dir, "IBDTest_ibd_segments_run1_year10.csv")
	rows := readCSV(t, segPath)
	if len(rows) != 1+3 { // header + 3 segments
		t.Errorf("segment CSV: expected 4 rows, got %d", len(rows))
	}
	wantHeader := []string{"Ind1", "Ind2", "Chrom", "Arm", "StartBit", "EndBit", "LengthBits"}
	if !equalRow(rows[0], wantHeader) {
		t.Errorf("segment CSV header = %v, want %v", rows[0], wantHeader)
	}

	histPath := filepath.Join(dir, "IBDTest_ibd_lengthdist_run1_year10.csv")
	hrows := readCSV(t, histPath)
	if len(hrows) != 1+2 { // header + 2 occupied bins (32-bit, 64-bit)
		t.Errorf("length-dist CSV: expected 3 rows, got %d", len(hrows))
	}
}

func readCSV(t *testing.T, path string) [][]string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("could not open %s: %v", path, err)
	}
	defer f.Close()
	rows, err := csv.NewReader(f).ReadAll()
	if err != nil {
		t.Fatalf("could not read %s: %v", path, err)
	}
	return rows
}

func equalRow(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
