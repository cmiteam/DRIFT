package analysis

import (
	"path/filepath"
	"testing"
)

// Joint-SFS tests reuse fstInd / buildFstPop / makeIBDTestModel / readCSV /
// equalRow. Every fixture has a single polymorphic site (genome bit 0).

// hetInd is one individual heterozygous at bit0 (one strand copy derived) in the
// given deme — contributes exactly one derived allele.
func hetInd(deme int) ([]int, [][]uint64) { return fstInd(deme, bit0, 0, 0, 0) }

func TestComputeJointSFS_2D(t *testing.T) {
	model := makeIBDTestModel() // no joint_sfs_demes -> defaults to the eligible demes

	// deme 0: one homozygous-derived (2) + one het (1) = 3 derived of 4.
	// deme 1: one het (1) + one homozygous-ancestral (0) = 1 derived of 4.
	d0a, c0a := homDerived(0)
	d0b, c0b := hetInd(0)
	d1a, c1a := hetInd(1)
	d1b, c1b := homAncestral(1)
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b}})

	res := ComputeJointSFS(model, pop, ids)
	if res.Dim != 2 {
		t.Fatalf("Dim = %d, want 2", res.Dim)
	}
	if res.NumSites != 1 {
		t.Fatalf("NumSites = %d, want 1", res.NumSites)
	}
	if len(res.Sizes) != 2 || res.Sizes[0] != 4 || res.Sizes[1] != 4 {
		t.Fatalf("Sizes = %v, want [4 4]", res.Sizes)
	}
	// The single site lands in cell (3,1).
	if got := res.Cells[[3]int{3, 1, 0}]; got != 1 {
		t.Errorf("cell (3,1) = %d, want 1 (cells: %v)", got, res.Cells)
	}
}

func TestComputeJointSFS_3D(t *testing.T) {
	model := makeIBDTestModel()

	// Three demes, derived counts 2 / 1 / 0 at bit0.
	d0a, c0a := homDerived(0)
	d0b, c0b := homAncestral(0) // deme0: 2 derived
	d1a, c1a := hetInd(1)
	d1b, c1b := homAncestral(1) // deme1: 1 derived
	d2a, c2a := homAncestral(2)
	d2b, c2b := homAncestral(2) // deme2: 0 derived
	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b}, {d2a, c2a}, {d2b, c2b},
	})

	res := ComputeJointSFS(model, pop, ids)
	if res.Dim != 3 {
		t.Fatalf("Dim = %d, want 3", res.Dim)
	}
	if got := res.Cells[[3]int{2, 1, 0}]; got != 1 {
		t.Errorf("cell (2,1,0) = %d, want 1 (cells: %v)", got, res.Cells)
	}
}

func TestComputeJointSFS_ExplicitDemes(t *testing.T) {
	model := makeIBDTestModel()
	model.StringParams = map[string]string{"joint_sfs_demes": "0,2"} // skip deme 1

	d0a, c0a := homDerived(0)
	d0b, c0b := homDerived(0) // deme0: 4 derived
	d1a, c1a := hetInd(1)
	d1b, c1b := homAncestral(1)
	d2a, c2a := homAncestral(2)
	d2b, c2b := homAncestral(2) // deme2: 0 derived
	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b}, {d2a, c2a}, {d2b, c2b},
	})

	res := ComputeJointSFS(model, pop, ids)
	if res.Dim != 2 {
		t.Fatalf("Dim = %d, want 2 (two named demes)", res.Dim)
	}
	if len(res.Demes) != 2 || res.Demes[0] != 0 || res.Demes[1] != 2 {
		t.Fatalf("axis demes = %v, want [0 2]", res.Demes)
	}
	// deme0 has 4 derived, deme2 has 0 -> cell (4,0).
	if got := res.Cells[[3]int{4, 0, 0}]; got != 1 {
		t.Errorf("cell (4,0) = %d, want 1 (cells: %v)", got, res.Cells)
	}
}

func TestComputeJointSFS_SingleDemeGuard(t *testing.T) {
	model := makeIBDTestModel()
	d0a, c0a := homDerived(0)
	d0b, c0b := homAncestral(0)
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}})

	res := ComputeJointSFS(model, pop, ids)
	if res.Dim != 0 || res.NumSites != 0 {
		t.Errorf("single-deme result should be empty, got Dim %d / %d sites", res.Dim, res.NumSites)
	}
}

func TestSaveJointSFS_WritesLongCSV(t *testing.T) {
	dir := t.TempDir()
	model := makeIBDTestModel()
	model.ResultsDir = dir

	d0a, c0a := homDerived(0)
	d0b, c0b := hetInd(0)
	d1a, c1a := hetInd(1)
	d1b, c1b := homAncestral(1)
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b}})

	res := ComputeJointSFS(model, pop, ids)
	if err := SaveJointSFS(model, res); err != nil {
		t.Fatalf("SaveJointSFS failed: %v", err)
	}

	path := filepath.Join(dir, "IBDTest_joint_sfs_run1_year10.csv")
	rows := readCSV(t, path)
	// header + one non-empty cell.
	if len(rows) != 2 {
		t.Fatalf("joint SFS CSV: expected 2 rows, got %d: %v", len(rows), rows)
	}
	wantHeader := []string{"Derived_deme0", "Derived_deme1", "Count"}
	if !equalRow(rows[0], wantHeader) {
		t.Errorf("joint SFS header = %v, want %v", rows[0], wantHeader)
	}
	if !equalRow(rows[1], []string{"3", "1", "1"}) {
		t.Errorf("joint SFS data row = %v, want [3 1 1]", rows[1])
	}
}
