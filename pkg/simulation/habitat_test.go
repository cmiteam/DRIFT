package simulation

import (
	"testing"

	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
)

func TestHabitatMortalityFactor(t *testing.T) {
	cases := []struct {
		suitability, want float64
	}{
		{1.0, 1.0},  // fully suitable -> no change
		{0.5, 2.0},  // half suitability -> double hazard
		{0.25, 4.0}, // quarter -> quadruple
		{2.0, 0.5},  // refugium -> lower hazard
		{0.0, 1.0},  // guard: non-positive -> no change
		{-1.0, 1.0}, // guard
	}
	for _, c := range cases {
		if got := habitatMortalityFactor(c.suitability); got != c.want {
			t.Errorf("habitatMortalityFactor(%v) = %v, want %v", c.suitability, got, c.want)
		}
	}
}

// habitatModel builds a Model with a 3-cell row map (Land=1, Land=1, Desert=5) and
// a suitability table that makes the Desert cell harsh.
func habitatModel(spec string) *core.Model {
	return &core.Model{
		StringParams: map[string]string{"habitat_suitability": spec},
		Map:          [][]int{{1, 1, 5}},
	}
}

// placeAt adds n individuals at the given cell starting at id `start`, returning the
// ids used. Each carries a full IndData slice with Lat/Lon set.
func placeAt(pop *core.Pop, start, n, lat, lon int) []int {
	ids := make([]int, n)
	for i := 0; i < n; i++ {
		id := start + i
		data := individual.MakeIndData()
		data[individual.Lat] = lat
		data[individual.Lon] = lon
		pop.IndData[id] = data
		ids[i] = id
	}
	return ids
}

func TestBuildCellCounts(t *testing.T) {
	pop := &core.Pop{IndData: map[int][]int{}}
	placeAt(pop, 0, 3, 0, 0)
	placeAt(pop, 100, 2, 0, 2)
	counts := buildCellCounts(pop)
	if counts[[2]int{0, 0}] != 3 {
		t.Errorf("cell (0,0) count = %d, want 3", counts[[2]int{0, 0}])
	}
	if counts[[2]int{0, 2}] != 2 {
		t.Errorf("cell (0,2) count = %d, want 2", counts[[2]int{0, 2}])
	}
}

// TestPickVictimsHabitatFavorsPoorCells is the distinguishable-outcome test: with a
// good cell (Land, suitability 1) and an equally-crowded poor cell (Desert,
// suitability 0.2), the weighted cull must remove far more from the poor cell.
func TestPickVictimsHabitatFavorsPoorCells(t *testing.T) {
	m := habitatModel("5:0.2")
	if !m.HabitatActive() {
		t.Fatal("expected active habitat table")
	}

	buildPop := func() (*core.Pop, map[int]bool) {
		pop := &core.Pop{IndData: map[int][]int{}}
		good := placeAt(pop, 0, 20, 0, 0) // Land cell (0,0), suitability 1.0
		placeAt(pop, 100, 20, 0, 2)       // Desert cell (0,2), suitability 0.2
		goodSet := make(map[int]bool, len(good))
		for _, id := range good {
			goodSet[id] = true
		}
		return pop, goodSet
	}

	candidates := make([]int, 0, 40)
	for i := 0; i < 20; i++ {
		candidates = append(candidates, i)     // good
		candidates = append(candidates, 100+i) // poor
	}
	// generateKeyList sorts; mimic that so the draw order matches production.
	sorted := append([]int(nil), candidates...)
	// simple insertion into sorted order
	for i := 1; i < len(sorted); i++ {
		for j := i; j > 0 && sorted[j-1] > sorted[j]; j-- {
			sorted[j-1], sorted[j] = sorted[j], sorted[j-1]
		}
	}

	utils.SeedRNG(4242)
	pop, goodSet := buildPop()
	victims := pickVictimsHabitat(m, pop, 20, sorted)
	if len(victims) != 20 {
		t.Fatalf("got %d victims, want 20", len(victims))
	}
	goodCulled, poorCulled := 0, 0
	for _, id := range victims {
		if goodSet[id] {
			goodCulled++
		} else {
			poorCulled++
		}
	}
	if poorCulled <= goodCulled {
		t.Errorf("poor-cell culls %d not greater than good-cell culls %d (weighting failed)", poorCulled, goodCulled)
	}
	// With a 5:1 weight ratio the poor cell should dominate the 20 victims.
	if poorCulled < 13 {
		t.Errorf("poor-cell culls %d, expected the majority (>=13) with a 5:1 weight", poorCulled)
	}
}

func TestPickVictimsHabitatDeterministic(t *testing.T) {
	m := habitatModel("5:0.2")
	pop := &core.Pop{IndData: map[int][]int{}}
	placeAt(pop, 0, 10, 0, 0)
	placeAt(pop, 100, 10, 0, 2)
	cands := []int{0, 1, 2, 3, 4, 100, 101, 102, 103, 104}

	utils.SeedRNG(777)
	a := pickVictimsHabitat(m, pop, 5, cands)
	utils.SeedRNG(777)
	b := pickVictimsHabitat(m, pop, 5, cands)
	if len(a) != len(b) {
		t.Fatalf("lengths differ: %d vs %d", len(a), len(b))
	}
	for i := range a {
		if a[i] != b[i] {
			t.Fatalf("victim %d differs: %d vs %d (non-deterministic)", i, a[i], b[i])
		}
	}
}

func TestPickVictimsHabitatCullAll(t *testing.T) {
	m := habitatModel("5:0.2")
	pop := &core.Pop{IndData: map[int][]int{}}
	placeAt(pop, 0, 3, 0, 0)
	cands := []int{0, 1, 2}
	utils.SeedRNG(1)
	got := pickVictimsHabitat(m, pop, 5, cands) // count >= len -> all
	if len(got) != 3 {
		t.Fatalf("cull-all returned %d, want 3", len(got))
	}
}
