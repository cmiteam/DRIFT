package analysis

import (
	"math"
	"path/filepath"
	"testing"

	"drift/pkg/core"
	"drift/pkg/individual"
)

// fstInd builds one individual's IndData (with a deme label) and a 128-bit
// chromosome (two 64-bit words per strand copy, matching makeIBDTestModel's arms).
func fstInd(deme int, c0w0, c0w1, c1w0, c1w1 uint64) ([]int, [][]uint64) {
	data := individual.MakeIndData()
	data[individual.Deme] = deme
	chr := [][]uint64{{c0w0, c0w1}, {c1w0, c1w1}}
	return data, chr
}

// buildFstPop assembles a Pop from per-individual (data, chromosome) specs, keyed
// by sequential IDs starting at 1.
func buildFstPop(specs [][2]interface{}) (*core.Pop, []int) {
	pop := &core.Pop{
		IndData:     map[int][]int{},
		Chromosomes: map[int][][]uint64{},
	}
	ids := make([]int, 0, len(specs))
	for i, s := range specs {
		id := i + 1
		pop.IndData[id] = s[0].([]int)
		pop.Chromosomes[id] = s[1].([][]uint64)
		ids = append(ids, id)
	}
	return pop, ids
}

const bit0 = uint64(1) // only genome bit 0 (p arm, word 0) is set

// TestComputeFst_CompleteDifferentiation: deme 0 fixed-derived, deme 1 fixed-
// ancestral at a single site. Both estimators must equal 1.
func TestComputeFst_CompleteDifferentiation(t *testing.T) {
	model := makeIBDTestModel()

	d0a, c0a := fstInd(0, bit0, 0, bit0, 0) // homozygous derived
	d0b, c0b := fstInd(0, bit0, 0, bit0, 0)
	d1a, c1a := fstInd(1, 0, 0, 0, 0) // homozygous ancestral
	d1b, c1b := fstInd(1, 0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b},
	})

	res := ComputeFst(model, pop, ids)
	if res.NumSites != 1 {
		t.Fatalf("NumSites = %d, want 1", res.NumSites)
	}
	if len(res.Pairs) != 1 {
		t.Fatalf("expected 1 deme pair, got %d", len(res.Pairs))
	}
	p := res.Pairs[0]
	assertClose(t, "Hudson", p.Hudson, 1.0)
	assertClose(t, "WC", p.WC, 1.0)
	assertClose(t, "GlobalWC", res.GlobalWC, 1.0)
}

// TestComputeFst_NoDifferentiation: both demes all-heterozygous (p = 0.5
// everywhere). Weir & Cockerham θ is exactly 0 (no between-deme variance and
// observed het exactly balances pq); Hudson's estimator is negative (its
// within-sample correction with no real differentiation), which must not read as
// positive structure.
func TestComputeFst_NoDifferentiation(t *testing.T) {
	model := makeIBDTestModel()

	d0a, c0a := fstInd(0, bit0, 0, 0, 0) // heterozygous
	d0b, c0b := fstInd(0, bit0, 0, 0, 0)
	d1a, c1a := fstInd(1, bit0, 0, 0, 0)
	d1b, c1b := fstInd(1, bit0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b},
	})

	res := ComputeFst(model, pop, ids)
	if res.NumSites != 1 {
		t.Fatalf("NumSites = %d, want 1", res.NumSites)
	}
	p := res.Pairs[0]
	assertClose(t, "WC", p.WC, 0.0)
	assertClose(t, "GlobalWC", res.GlobalWC, 0.0)
	if p.Hudson > 1e-9 {
		t.Errorf("Hudson = %v, want <= 0 (no differentiation)", p.Hudson)
	}
}

// TestComputeFst_Intermediate: deme 0 fixed-derived (p=1), deme 1 all-
// heterozygous (p=0.5). Hand-computed at the single site: Hudson = 1/3, WC = 1/2.
func TestComputeFst_Intermediate(t *testing.T) {
	model := makeIBDTestModel()

	d0a, c0a := fstInd(0, bit0, 0, bit0, 0) // homozygous derived
	d0b, c0b := fstInd(0, bit0, 0, bit0, 0)
	d1a, c1a := fstInd(1, bit0, 0, 0, 0) // heterozygous
	d1b, c1b := fstInd(1, bit0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b},
	})

	res := ComputeFst(model, pop, ids)
	p := res.Pairs[0]
	assertClose(t, "Hudson", p.Hudson, 1.0/3.0)
	assertClose(t, "WC", p.WC, 0.5)
}

// TestComputeFst_SingleDemeGuard: with only one deme present, the result is empty
// (Fst undefined) and nothing panics.
func TestComputeFst_SingleDemeGuard(t *testing.T) {
	model := makeIBDTestModel()

	d0a, c0a := fstInd(0, bit0, 0, bit0, 0)
	d0b, c0b := fstInd(0, bit0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}})

	res := ComputeFst(model, pop, ids)
	if len(res.Pairs) != 0 || res.NumSites != 0 {
		t.Errorf("single-deme result should be empty, got %d pairs / %d sites", len(res.Pairs), res.NumSites)
	}
}

// TestComputeFst_SmallDemeExcluded: a deme with fewer than 2 sampled diploids is
// dropped; with only one eligible deme left the result is empty.
func TestComputeFst_SmallDemeExcluded(t *testing.T) {
	model := makeIBDTestModel()

	d0a, c0a := fstInd(0, bit0, 0, bit0, 0)
	d0b, c0b := fstInd(0, bit0, 0, bit0, 0)
	d1a, c1a := fstInd(1, 0, 0, 0, 0) // lone member of deme 1 -> excluded
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}, {d1a, c1a}})

	res := ComputeFst(model, pop, ids)
	if len(res.Pairs) != 0 {
		t.Errorf("expected no pairs (deme 1 excluded), got %d", len(res.Pairs))
	}
}

// TestSaveFst_WritesReadableCSV mirrors TestSaveIBD_WritesReadableCSVs: header +
// one row per pair + the global-θ row.
func TestSaveFst_WritesReadableCSV(t *testing.T) {
	dir := t.TempDir()
	model := makeIBDTestModel()
	model.ResultsDir = dir

	d0a, c0a := fstInd(0, bit0, 0, bit0, 0)
	d0b, c0b := fstInd(0, bit0, 0, bit0, 0)
	d1a, c1a := fstInd(1, 0, 0, 0, 0)
	d1b, c1b := fstInd(1, 0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b}})

	res := ComputeFst(model, pop, ids)
	if err := SaveFst(model, res); err != nil {
		t.Fatalf("SaveFst failed: %v", err)
	}

	path := filepath.Join(dir, "IBDTest_fst_run1_year10.csv")
	rows := readCSV(t, path)
	// header + 1 pair + global row = 3
	if len(rows) != 3 {
		t.Errorf("Fst CSV: expected 3 rows, got %d: %v", len(rows), rows)
	}
	wantHeader := []string{"DemeA", "DemeB", "Fst_Hudson", "Fst_WC", "NumSites"}
	if !equalRow(rows[0], wantHeader) {
		t.Errorf("Fst CSV header = %v, want %v", rows[0], wantHeader)
	}
	if rows[2][0] != "ALL" {
		t.Errorf("expected global row marked ALL, got %v", rows[2])
	}
}

func assertClose(t *testing.T, name string, got, want float64) {
	t.Helper()
	if math.Abs(got-want) > 1e-9 {
		t.Errorf("%s = %v, want %v", name, got, want)
	}
}
