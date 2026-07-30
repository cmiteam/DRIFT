package core

import "testing"

// migModel builds a Model whose migration_matrix param is spec, so the lazy accessor
// parses that spec.
func migModel(spec string) *Model {
	return &Model{StringParams: map[string]string{"migration_matrix": spec}}
}

// findEdge returns the rate of the from→to edge, or -1 if absent.
func findEdge(mm *MigrationMatrix, from, to int) float64 {
	for _, e := range mm.Out[from] {
		if e.Dest == to {
			return e.Rate
		}
	}
	return -1
}

func TestParseMigrationMatrix_Basic(t *testing.T) {
	// Directional, asymmetric, multiple separators (space / ';' / ',').
	mm := parseMigrationMatrix("0>1:0.01 1>0:0.005;1>2:0.02,2>1:0.03")
	if !mm.Active {
		t.Fatal("expected active matrix")
	}
	cases := []struct {
		from, to int
		want     float64
	}{
		{0, 1, 0.01},
		{1, 0, 0.005},
		{1, 2, 0.02},
		{2, 1, 0.03},
	}
	for _, c := range cases {
		if got := findEdge(mm, c.from, c.to); got != c.want {
			t.Errorf("edge %d>%d = %v, want %v", c.from, c.to, got, c.want)
		}
	}
	// Asymmetry: 0>1 exists but 1>0 has a different rate; 2>0 does not exist.
	if findEdge(mm, 0, 1) == findEdge(mm, 1, 0) {
		t.Error("expected asymmetric 0>1 vs 1>0 rates")
	}
	if findEdge(mm, 2, 0) != -1 {
		t.Error("edge 2>0 should be absent")
	}
}

func TestParseMigrationMatrix_EdgesSortedByDest(t *testing.T) {
	// Edges out of deme 0 authored out of Dest order must come back Dest-sorted so the
	// RNG walk is deterministic.
	mm := parseMigrationMatrix("0>3:0.1 0>1:0.1 0>2:0.1")
	edges := mm.Out[0]
	if len(edges) != 3 {
		t.Fatalf("expected 3 outgoing edges, got %d", len(edges))
	}
	for i := 1; i < len(edges); i++ {
		if edges[i-1].Dest >= edges[i].Dest {
			t.Errorf("edges not Dest-sorted: %v", edges)
		}
	}
}

func TestParseMigrationMatrix_MalformedSkippedSelfAndZeroDropped(t *testing.T) {
	// "junk" (no arrow/colon), "5>x:0.1" (bad dest), "6>7:abc" (bad rate), "3>3:0.1"
	// (self-edge), "4>5:0" (zero rate) are all dropped; only 0>1:0.2 survives.
	mm := parseMigrationMatrix("junk 5>x:0.1 6>7:abc 3>3:0.1 4>5:0 0>1:0.2")
	if !mm.Active {
		t.Fatal("expected active matrix (0>1 survives)")
	}
	if findEdge(mm, 0, 1) != 0.2 {
		t.Errorf("0>1 = %v, want 0.2", findEdge(mm, 0, 1))
	}
	for _, bad := range [][2]int{{3, 3}, {4, 5}, {5, 0}, {6, 7}} {
		if findEdge(mm, bad[0], bad[1]) != -1 {
			t.Errorf("edge %d>%d should have been dropped", bad[0], bad[1])
		}
	}
}

func TestParseMigrationMatrix_NegativeClampedDropped(t *testing.T) {
	// A negative rate is clamped to 0 and then dropped as a no-op edge.
	mm := parseMigrationMatrix("0>1:-0.5")
	if mm.Active {
		t.Error("a lone negative (→0) edge leaves the matrix inactive")
	}
}

func TestParseMigrationMatrix_DuplicateLastWins(t *testing.T) {
	mm := parseMigrationMatrix("0>1:0.1 0>1:0.4")
	if findEdge(mm, 0, 1) != 0.4 {
		t.Errorf("0>1 = %v, want 0.4 (last wins)", findEdge(mm, 0, 1))
	}
}

func TestParseMigrationMatrix_EmptyInactive(t *testing.T) {
	for _, spec := range []string{"", "   ", ";,; ", "garbage tokens only"} {
		if mm := parseMigrationMatrix(spec); mm.Active {
			t.Errorf("spec %q should be inactive", spec)
		}
	}
}

func TestMigrationMatrixAccessor_LazyCached(t *testing.T) {
	m := migModel("0>1:0.05")
	if m.Migration != nil {
		t.Fatal("Migration should be nil before first access")
	}
	mm := m.MigrationMatrix()
	if !mm.Active || findEdge(mm, 0, 1) != 0.05 {
		t.Fatalf("accessor returned wrong matrix: %+v", mm)
	}
	if m.Migration == nil {
		t.Error("Migration should be cached after first access")
	}
	// Unset param => inactive.
	if (&Model{StringParams: map[string]string{}}).MigrationMatrix().Active {
		t.Error("unset migration_matrix should be inactive")
	}
}
