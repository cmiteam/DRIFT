package core

import "testing"

// modelWithMap builds a minimal Model carrying a terrain grid and a habitat param,
// for exercising the cell accessors. The map uses codes Land=1, OpenWater=3,
// Desert=5 (matching the utils.Terrain enum).
func modelWithMap(spec string) *Model {
	return &Model{
		StringParams: map[string]string{"habitat_suitability": spec},
		Map: [][]int{
			{1, 3, 5}, // row 0: Land, OpenWater, Desert
			{1, 1, 3}, // row 1
		},
	}
}

func TestParseHabitatTable(t *testing.T) {
	// Separators: space, ';' and ',' all accepted; malformed tokens skipped;
	// negative multipliers clamped to 0.
	tab := parseHabitatTable("1:1.0 5:0.3;6:0.05,4:-2 garbage 3")
	if !tab.Active {
		t.Fatal("expected active table")
	}
	want := map[int]float64{1: 1.0, 5: 0.3, 6: 0.05, 4: 0.0}
	if len(tab.Mult) != len(want) {
		t.Fatalf("Mult has %d entries, want %d: %v", len(tab.Mult), len(want), tab.Mult)
	}
	for code, w := range want {
		if got := tab.Mult[code]; got != w {
			t.Errorf("code %d = %v, want %v", code, got, w)
		}
	}
}

func TestParseHabitatTableEmptyInactive(t *testing.T) {
	for _, spec := range []string{"", "   ", ";;", "junk nototokens"} {
		tab := parseHabitatTable(spec)
		if tab.Active {
			t.Errorf("parseHabitatTable(%q).Active = true, want inactive", spec)
		}
	}
}

func TestCellSuitabilityInactiveIsBinary(t *testing.T) {
	// No habitat param: the accessor must reproduce IsLand exactly — 1.0 on Land,
	// 0.0 everywhere else, 0.0 out of bounds.
	m := modelWithMap("")
	cases := []struct {
		lat, lon int
		want     float64
	}{
		{0, 0, 1.0}, // Land
		{0, 1, 0.0}, // OpenWater
		{0, 2, 0.0}, // Desert (unlisted -> uninhabitable when inactive)
		{1, 0, 1.0}, // Land
		{5, 5, 0.0}, // out of bounds
		{-1, 0, 0.0},
	}
	for _, c := range cases {
		if got := m.CellSuitability(c.lat, c.lon); got != c.want {
			t.Errorf("CellSuitability(%d,%d) = %v, want %v", c.lat, c.lon, got, c.want)
		}
		if got := m.IsHabitable(c.lat, c.lon); got != (c.want > 0) {
			t.Errorf("IsHabitable(%d,%d) = %v, want %v", c.lat, c.lon, got, c.want > 0)
		}
	}
	if m.HabitatActive() {
		t.Error("HabitatActive() = true for empty param")
	}
}

func TestCellSuitabilityActive(t *testing.T) {
	// Desert given 0.3; Land unlisted -> stays 1.0; OpenWater unlisted -> 0.0.
	m := modelWithMap("5:0.3")
	if !m.HabitatActive() {
		t.Fatal("HabitatActive() = false, want true")
	}
	checks := []struct {
		lat, lon int
		want     float64
	}{
		{0, 2, 0.3}, // Desert -> configured
		{0, 0, 1.0}, // Land -> unlisted default 1.0
		{0, 1, 0.0}, // OpenWater -> unlisted default 0.0
	}
	for _, c := range checks {
		if got := m.CellSuitability(c.lat, c.lon); got != c.want {
			t.Errorf("CellSuitability(%d,%d) = %v, want %v", c.lat, c.lon, got, c.want)
		}
	}
	// A Desert cell is now habitable (suitability > 0) where it was not before.
	if !m.IsHabitable(0, 2) {
		t.Error("Desert with suitability 0.3 should be habitable")
	}
}

func TestCellSuitabilityActiveCanDowngradeLand(t *testing.T) {
	// An explicit Land entry overrides the 1.0 default (a harsh homeland).
	m := modelWithMap("1:0.5")
	if got := m.CellSuitability(0, 0); got != 0.5 {
		t.Errorf("downgraded Land suitability = %v, want 0.5", got)
	}
}
