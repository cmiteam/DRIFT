package core

import "testing"

// barrierMap builds a 5x5 terrain grid with an OpenWater(3) barrier wall down
// column 2, broken by a Land(1) corridor at row 2 — a strait with a land bridge:
//
//	row0: 1 1 3 1 1
//	row1: 1 1 3 1 1
//	row2: 1 1 1 1 1   <- corridor (gap in the wall)
//	row3: 1 1 3 1 1
//	row4: 1 1 3 1 1
func barrierModel(spec string) *Model {
	return &Model{
		StringParams: map[string]string{"movement_barriers": spec},
		Map: [][]int{
			{1, 1, 3, 1, 1},
			{1, 1, 3, 1, 1},
			{1, 1, 1, 1, 1},
			{1, 1, 3, 1, 1},
			{1, 1, 3, 1, 1},
		},
	}
}

func TestParseBarrierTable(t *testing.T) {
	// Bare code, code:block, and a passable override; separators space/';'/','; a
	// malformed token skipped. Codes 3 and 4 block; 1 is explicitly passable (so it
	// is NOT recorded); "garbage" is skipped.
	tab := parseBarrierTable("3 4:block;1:pass,garbage")
	if !tab.Active {
		t.Fatal("expected active table")
	}
	if !tab.Block[3] || !tab.Block[4] {
		t.Errorf("codes 3 and 4 should block: %v", tab.Block)
	}
	if tab.Block[1] {
		t.Error("code 1 marked :pass should not block")
	}
	if len(tab.Block) != 2 {
		t.Errorf("Block has %d entries, want 2: %v", len(tab.Block), tab.Block)
	}
}

func TestParseBarrierTableEmptyInactive(t *testing.T) {
	// Empty, whitespace, only-separators, no valid codes, and a table whose only
	// entries are explicitly passable all collapse to inactive.
	for _, spec := range []string{"", "   ", ";;", "nope junk", "1:pass 3:open"} {
		tab := parseBarrierTable(spec)
		if tab.Active {
			t.Errorf("parseBarrierTable(%q).Active = true, want inactive", spec)
		}
	}
}

func TestIsBarrierInactive(t *testing.T) {
	// No barrier param: nothing is ever a barrier, including the water cells.
	m := barrierModel("")
	if m.MovementBarriersActive() {
		t.Error("MovementBarriersActive() = true for empty param")
	}
	for _, c := range [][2]int{{0, 2}, {1, 2}, {0, 0}, {9, 9}} {
		if m.IsBarrier(c[0], c[1]) {
			t.Errorf("IsBarrier(%d,%d) = true with no table", c[0], c[1])
		}
	}
}

func TestIsBarrierActive(t *testing.T) {
	m := barrierModel("3")
	if !m.MovementBarriersActive() {
		t.Fatal("MovementBarriersActive() = false, want true")
	}
	cases := []struct {
		lat, lon int
		want     bool
	}{
		{0, 2, true},   // OpenWater in the wall
		{2, 2, false},  // corridor Land
		{0, 0, false},  // Land
		{-1, 0, false}, // out of bounds
		{0, 9, false},  // out of bounds
	}
	for _, c := range cases {
		if got := m.IsBarrier(c.lat, c.lon); got != c.want {
			t.Errorf("IsBarrier(%d,%d) = %v, want %v", c.lat, c.lon, got, c.want)
		}
	}
}

func TestCanTraverseInactive(t *testing.T) {
	// With no barrier table, every move is allowed regardless of what it crosses —
	// this is what keeps Wander's accept predicate (and the whole run) byte-identical.
	m := barrierModel("")
	moves := [][4]int{
		{0, 0, 0, 4}, // straight across the (would-be) wall
		{0, 2, 4, 2}, // straight down the wall
		{0, 0, 4, 4}, // diagonal
	}
	for _, mv := range moves {
		if !m.CanTraverse(mv[0], mv[1], mv[2], mv[3]) {
			t.Errorf("CanTraverse%v = false with no table, want true", mv)
		}
	}
}

func TestCanTraverseLOS(t *testing.T) {
	m := barrierModel("3") // OpenWater blocks crossing
	cases := []struct {
		fromLat, fromLon, toLat, toLon int
		want                           bool
		note                           string
	}{
		{0, 0, 0, 4, false, "horizontal at row0 crosses wall cell (0,2)"},
		{0, 0, 0, 1, true, "adjacent, no interior cell"},
		{2, 0, 2, 4, true, "horizontal along the corridor row (all Land)"},
		{0, 0, 4, 0, true, "vertical down column 0 (all Land)"},
		{0, 2, 4, 2, false, "along the wall: interior (1,2) is a barrier"},
		{0, 0, 2, 2, true, "diagonal to the corridor cell; interior (1,1) is Land"},
		{0, 1, 2, 3, false, "diagonal crossing wall cell (1,2)"},
		{3, 3, 3, 3, true, "no movement"},
	}
	for _, c := range cases {
		got := m.CanTraverse(c.fromLat, c.fromLon, c.toLat, c.toLon)
		if got != c.want {
			t.Errorf("CanTraverse(%d,%d->%d,%d) = %v, want %v (%s)",
				c.fromLat, c.fromLon, c.toLat, c.toLon, got, c.want, c.note)
		}
	}
}

func TestCanTraverseSymmetric(t *testing.T) {
	// The crossing test is a property of the line, so it must be direction-agnostic.
	m := barrierModel("3")
	pairs := [][4]int{{0, 0, 0, 4}, {0, 1, 2, 3}, {2, 0, 2, 4}}
	for _, p := range pairs {
		fwd := m.CanTraverse(p[0], p[1], p[2], p[3])
		rev := m.CanTraverse(p[2], p[3], p[0], p[1])
		if fwd != rev {
			t.Errorf("CanTraverse asymmetric for %v: fwd=%v rev=%v", p, fwd, rev)
		}
	}
}
