package analysis

import (
	"path/filepath"
	"testing"
)

// The f-stats tests reuse fstInd / buildFstPop / makeIBDTestModel / assertClose /
// readCSV / equalRow from fst_test.go and ibd_test.go (same package). bit0 (genome
// bit 0) is the single polymorphic site in every fixture.

// homozygous-derived / -ancestral single-site individuals in a given deme.
func homDerived(deme int) ([]int, [][]uint64)   { return fstInd(deme, bit0, 0, bit0, 0) }
func homAncestral(deme int) ([]int, [][]uint64) { return fstInd(deme, 0, 0, 0, 0) }

// --- Per-site kernel unit tests (exact) ---------------------------------------

func TestF2Site_Kernel(t *testing.T) {
	// Complete differentiation, n=4 haplotypes each: (1-0)^2 - 0 - 0 = 1.
	assertClose(t, "f2", f2Site(1.0, 4, 0.0, 4), 1.0)
	// Identical frequencies -> negative (pure sampling correction).
	if f2Site(0.5, 4, 0.5, 4) >= 0 {
		t.Errorf("f2Site identical demes should be < 0, got %v", f2Site(0.5, 4, 0.5, 4))
	}
}

func TestF3Site_Kernel(t *testing.T) {
	// Target p=0.5 midway between B(p=1) and C(p=0), n_T=4:
	// (0.5-1)(0.5-0) - 0.5*0.5/3 = -0.25 - 1/12 = -1/3.
	assertClose(t, "f3", f3Site(0.5, 4, 1.0, 0.0), -1.0/3.0)
}

func TestF4AndABBA_Kernel(t *testing.T) {
	// A=0,B=1,C=1,D=0: f4 = (0-1)(1-0) = -1.
	assertClose(t, "f4", f4Site(0, 1, 1, 0), -1.0)
	abba, baba := abbaBaba(0, 1, 1, 0)
	assertClose(t, "abba", abba, 1.0)
	assertClose(t, "baba", baba, 0.0)
	// Identity ABBA-BABA = (p2-p1)(p3-p4).
	assertClose(t, "identity", abba-baba, (1.0-0.0)*(1.0-0.0))
}

// --- End-to-end ComputeFStats -------------------------------------------------

func TestComputeFStats_F2CompleteDifferentiation(t *testing.T) {
	model := makeIBDTestModel()
	d0a, c0a := homDerived(0)
	d0b, c0b := homDerived(0)
	d1a, c1a := homAncestral(1)
	d1b, c1b := homAncestral(1)
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b}})

	res := ComputeFStats(model, pop, ids)
	if res.NumSites != 1 {
		t.Fatalf("NumSites = %d, want 1", res.NumSites)
	}
	if len(res.F2) != 1 {
		t.Fatalf("expected 1 f2 pair, got %d", len(res.F2))
	}
	assertClose(t, "f2(0,1)", res.F2[0].Value, 1.0)
}

func TestComputeFStats_F3Admixture(t *testing.T) {
	model := makeIBDTestModel()
	model.StringParams = map[string]string{"fstats_demes": "0;1,2"}

	// deme 0 (target) p=0.5, deme 1 p=1, deme 2 p=0.
	d0a, c0a := homDerived(0)   // 2 derived
	d0b, c0b := homAncestral(0) // 0 derived -> deme0 p = 2/4 = 0.5
	d1a, c1a := homDerived(1)
	d1b, c1b := homDerived(1)
	d2a, c2a := homAncestral(2)
	d2b, c2b := homAncestral(2)
	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b}, {d2a, c2a}, {d2b, c2b},
	})

	res := ComputeFStats(model, pop, ids)
	if len(res.F3) != 1 {
		t.Fatalf("explicit spec should give exactly 1 f3, got %d", len(res.F3))
	}
	f3 := res.F3[0]
	if f3.Target != 0 || f3.B != 1 || f3.C != 2 {
		t.Errorf("f3 slots = (%d;%d,%d), want (0;1,2)", f3.Target, f3.B, f3.C)
	}
	assertClose(t, "f3(0;1,2)", f3.Value, -1.0/3.0)
	// Explicit mode computes no f4.
	if len(res.F4) != 0 {
		t.Errorf("explicit f3 spec should yield no f4, got %d", len(res.F4))
	}
}

func TestComputeFStats_F4AndD(t *testing.T) {
	model := makeIBDTestModel()
	model.StringParams = map[string]string{"fstats_demes": "0,1;2,3"}

	// A=0 p=0, B=1 p=1, C=2 p=1, D=3 p=0 -> clean ABBA: f4=-1, D=+1.
	d0a, c0a := homAncestral(0)
	d0b, c0b := homAncestral(0)
	d1a, c1a := homDerived(1)
	d1b, c1b := homDerived(1)
	d2a, c2a := homDerived(2)
	d2b, c2b := homDerived(2)
	d3a, c3a := homAncestral(3)
	d3b, c3b := homAncestral(3)
	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b},
		{d2a, c2a}, {d2b, c2b}, {d3a, c3a}, {d3b, c3b},
	})

	res := ComputeFStats(model, pop, ids)
	if len(res.F4) != 1 {
		t.Fatalf("explicit spec should give exactly 1 f4/D, got %d", len(res.F4))
	}
	q := res.F4[0]
	if q.A != 0 || q.B != 1 || q.C != 2 || q.D != 3 {
		t.Errorf("f4 slots = (%d,%d;%d,%d), want (0,1;2,3)", q.A, q.B, q.C, q.D)
	}
	assertClose(t, "f4", q.F4, -1.0)
	assertClose(t, "D", q.Dstat, 1.0)
}

func TestComputeFStats_AutoEnumerationCounts(t *testing.T) {
	model := makeIBDTestModel() // no StringParams -> auto mode

	// Four demes with differing frequencies so bit0 is polymorphic.
	specs := [][2]interface{}{}
	add := func(fn func(int) ([]int, [][]uint64), d int) {
		data, chr := fn(d)
		specs = append(specs, [2]interface{}{data, chr})
	}
	add(homDerived, 0)
	add(homDerived, 0)
	add(homAncestral, 1)
	add(homAncestral, 1)
	add(homDerived, 2)
	add(homAncestral, 2)
	add(homAncestral, 3)
	add(homDerived, 3)
	pop, ids := buildFstPop(specs)

	res := ComputeFStats(model, pop, ids)
	if res.Mode != "auto" {
		t.Errorf("Mode = %q, want auto", res.Mode)
	}
	// 4 demes: C(4,2)=6 f2, 4*C(3,2)=12 f3, C(4,4)*3=3 f4/D.
	if len(res.F2) != 6 {
		t.Errorf("f2 count = %d, want 6", len(res.F2))
	}
	if len(res.F3) != 12 {
		t.Errorf("f3 count = %d, want 12", len(res.F3))
	}
	if len(res.F4) != 3 {
		t.Errorf("f4 count = %d, want 3", len(res.F4))
	}
}

func TestComputeFStats_SingleDemeGuard(t *testing.T) {
	model := makeIBDTestModel()
	d0a, c0a := homDerived(0)
	d0b, c0b := homAncestral(0)
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}})

	res := ComputeFStats(model, pop, ids)
	if len(res.F2) != 0 || len(res.F3) != 0 || len(res.F4) != 0 {
		t.Errorf("single-deme result should be empty, got %d/%d/%d f2/f3/f4",
			len(res.F2), len(res.F3), len(res.F4))
	}
}

func TestComputeFStats_MalformedSpecFallsBackToF2Only(t *testing.T) {
	model := makeIBDTestModel()
	model.StringParams = map[string]string{"fstats_demes": "garbage"}
	d0a, c0a := homDerived(0)
	d0b, c0b := homDerived(0)
	d1a, c1a := homAncestral(1)
	d1b, c1b := homAncestral(1)
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b}})

	res := ComputeFStats(model, pop, ids)
	// f2 is always computed; malformed f3/f4 spec yields none but does not panic.
	if len(res.F2) != 1 {
		t.Errorf("f2 should still compute despite malformed spec, got %d", len(res.F2))
	}
	if len(res.F3) != 0 || len(res.F4) != 0 {
		t.Errorf("malformed spec should yield no f3/f4, got %d/%d", len(res.F3), len(res.F4))
	}
}

func TestSaveFStats_WritesUniformCSV(t *testing.T) {
	dir := t.TempDir()
	model := makeIBDTestModel()
	model.ResultsDir = dir // no StringParams -> auto mode

	specs := [][2]interface{}{}
	add := func(fn func(int) ([]int, [][]uint64), d int) {
		data, chr := fn(d)
		specs = append(specs, [2]interface{}{data, chr})
	}
	add(homDerived, 0)
	add(homDerived, 0)
	add(homAncestral, 1)
	add(homAncestral, 1)
	add(homDerived, 2)
	add(homAncestral, 2)
	add(homAncestral, 3)
	add(homDerived, 3)
	pop, ids := buildFstPop(specs)

	res := ComputeFStats(model, pop, ids)
	if err := SaveFStats(model, res); err != nil {
		t.Fatalf("SaveFStats failed: %v", err)
	}

	path := filepath.Join(dir, "IBDTest_fstats_run1_year10.csv")
	rows := readCSV(t, path)
	// header + 6 f2 + 12 f3 + 3 f4 rows + 3 D rows = 25.
	if len(rows) != 25 {
		t.Errorf("f-stats CSV: expected 25 rows, got %d", len(rows))
	}
	wantHeader := []string{"Stat", "DemeA", "DemeB", "DemeC", "DemeD", "Value", "NumSites"}
	if !equalRow(rows[0], wantHeader) {
		t.Errorf("f-stats CSV header = %v, want %v", rows[0], wantHeader)
	}
	// Count the statistic kinds present.
	kinds := map[string]int{}
	for _, r := range rows[1:] {
		kinds[r[0]]++
	}
	if kinds["f2"] != 6 || kinds["f3"] != 12 || kinds["f4"] != 3 || kinds["D"] != 3 {
		t.Errorf("row kind counts = %v, want f2:6 f3:12 f4:3 D:3", kinds)
	}
}
