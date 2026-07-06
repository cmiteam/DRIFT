package demography

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"math"
	"path/filepath"
	"testing"
)

// buildScenario returns a finalized 4-epoch scenario in known units for the
// unit-level checks (Q=10, generation_time=25, mirroring the shipped OoA file).
func buildScenario(t *testing.T) *Scenario {
	t.Helper()
	s := &Scenario{
		GenerationTime: 25,
		Rescale:        10,
		Epochs: []Epoch{
			{Name: "ancestral", DurationGen: 200, Demes: []sizeSpec{{ID: 0, Size: 7300}}},
			{Name: "african", DurationGen: 3200, Demes: []sizeSpec{{ID: 0, Size: 12300}}},
			{Name: "ooa", DurationGen: 4752,
				Split: &splitSpec{Parent: 0, Child: 1, ChildFounders: 2100},
				Demes: []sizeSpec{{ID: 0, Size: 12300}, {ID: 1, Size: 2100}},
				Migration: []migSpec{{A: 0, B: 1, Rate: 0.00025}}},
			{Name: "eurasian", DurationGen: 848,
				Split: &splitSpec{Parent: 1, Child: 2, ChildFounders: 510},
				Demes: []sizeSpec{{ID: 0, Size: 12300}, {ID: 1, Size: 1000, Growth: 0.004}, {ID: 2, Size: 510, Growth: 0.0055}}},
		},
	}
	if err := s.finalize(); err != nil {
		t.Fatalf("finalize: %v", err)
	}
	return s
}

func TestFinalizeTimeline(t *testing.T) {
	s := buildScenario(t)
	// span(gen) = round(dur/Q * gt): 200->500, 3200->8000, 4752->11880, 848->2120.
	wantStart := []int{0, 500, 8500, 20380}
	wantEnd := []int{500, 8500, 20380, 22500}
	for i := range s.Epochs {
		if s.Epochs[i].startYear != wantStart[i] || s.Epochs[i].endYear != wantEnd[i] {
			t.Errorf("epoch %d: got [%d,%d), want [%d,%d)", i,
				s.Epochs[i].startYear, s.Epochs[i].endYear, wantStart[i], wantEnd[i])
		}
	}
	if s.TotalYears() != 22500 {
		t.Errorf("TotalYears = %d, want 22500", s.TotalYears())
	}
}

func TestDemeRanks(t *testing.T) {
	// OoA-shaped chain (origin 0 -> split to 1 -> split to 2): ranks 0,1,2.
	s := buildScenario(t)
	ranks := s.DemeRanks()
	want := map[int]int{0: 0, 1: 1, 2: 2}
	for d, r := range want {
		if ranks[d] != r {
			t.Errorf("DemeRanks[%d] = %d, want %d", d, ranks[d], r)
		}
	}
	if len(ranks) != len(want) {
		t.Errorf("DemeRanks has %d demes, want %d: %v", len(ranks), len(want), ranks)
	}

	// A serial single-origin (Babel-shaped) chain 0->1->2->3 gives strictly
	// increasing ranks — the near-linear distance-from-origin axis.
	serial := &Scenario{
		GenerationTime: 25, Rescale: 10,
		Epochs: []Epoch{
			{Name: "origin", DurationGen: 100, Demes: []sizeSpec{{ID: 0, Size: 1000}}},
			{Name: "d1", DurationGen: 100, Split: &splitSpec{Parent: 0, Child: 1, ChildFounders: 500}, Demes: []sizeSpec{{ID: 0, Size: 1000}, {ID: 1, Size: 500}}},
			{Name: "d2", DurationGen: 100, Split: &splitSpec{Parent: 1, Child: 2, ChildFounders: 250}, Demes: []sizeSpec{{ID: 0, Size: 1000}, {ID: 1, Size: 500}, {ID: 2, Size: 250}}},
			{Name: "d3", DurationGen: 100, Split: &splitSpec{Parent: 2, Child: 3, ChildFounders: 125}, Demes: []sizeSpec{{ID: 0, Size: 1000}, {ID: 1, Size: 500}, {ID: 2, Size: 250}, {ID: 3, Size: 125}}},
		},
	}
	if err := serial.finalize(); err != nil {
		t.Fatalf("finalize: %v", err)
	}
	sr := serial.DemeRanks()
	for d := 0; d <= 3; d++ {
		if sr[d] != d {
			t.Errorf("serial DemeRanks[%d] = %d, want %d", d, sr[d], d)
		}
	}
}

func TestEpochAt(t *testing.T) {
	s := buildScenario(t)
	cases := []struct {
		year int
		name string
	}{
		{0, "ancestral"},
		{499, "ancestral"},
		{500, "african"},
		{8500, "ooa"},
		{20380, "eurasian"},
		{99999, "eurasian"}, // past the end clamps to the final epoch
	}
	for _, c := range cases {
		if got := s.epochAt(c.year).Name; got != c.name {
			t.Errorf("epochAt(%d) = %q, want %q", c.year, got, c.name)
		}
	}
}

func TestCapForRescale(t *testing.T) {
	s := buildScenario(t)
	// Ancestral deme 0: 7300/10 = 730, no growth.
	if got := s.capFor(sizeSpec{ID: 0, Size: 7300}, &s.Epochs[0], 0); got != 730 {
		t.Errorf("ancestral cap = %d, want 730", got)
	}
}

func TestCapForGrowth(t *testing.T) {
	s := buildScenario(t)
	e := &s.Epochs[3] // eurasian, starts year 20380
	spec := sizeSpec{ID: 1, Size: 1000, Growth: 0.004}
	// At epoch start: base only.
	if got := s.capFor(spec, e, e.startYear); got != 100 {
		t.Errorf("EUR founding cap = %d, want 100", got)
	}
	// At epoch end: grown by exp(r * duration_gen) = exp(0.004*848) then /Q.
	want := int(math.Round(1000.0 / 10 * math.Exp(0.004*848)))
	if got := s.capFor(spec, e, e.endYear); got != want {
		t.Errorf("EUR present cap = %d, want %d", got, want)
	}
	// Growth clamps at the scenario end (no runaway past present).
	if s.capFor(spec, e, e.endYear) != s.capFor(spec, e, e.endYear+5000) {
		t.Errorf("growth not clamped past scenario end")
	}
}

func TestPerYearMigration(t *testing.T) {
	s := buildScenario(t)
	// m * Q / gt = 0.00025 * 10 / 25 = 1e-4.
	if got := s.perYearMigration(0.00025); math.Abs(got-1e-4) > 1e-12 {
		t.Errorf("perYearMigration = %g, want 1e-4", got)
	}
	if s.perYearMigration(1e9) != 1 || s.perYearMigration(-1) != 0 {
		t.Errorf("perYearMigration not clamped to [0,1]")
	}
}

// makePop builds a population of n individuals all in deme 0.
func makePop(n int) *core.Pop {
	pop := &core.Pop{IndData: make(map[int][]int)}
	for i := 0; i < n; i++ {
		d := individual.MakeIndData()
		d[individual.Deme] = 0
		pop.IndData[i] = d
	}
	return pop
}

func demeCount(pop *core.Pop, deme int) int {
	n := 0
	for _, d := range pop.IndData {
		if d[individual.Deme] == deme {
			n++
		}
	}
	return n
}

func TestDoSplitRelabels(t *testing.T) {
	utils.SeedRNG(42)
	s := buildScenario(t)
	pop := makePop(100)
	// child_founders 2100/Q = 210, but capped at half the parent (50).
	s.doSplit(pop, &splitSpec{Parent: 0, Child: 1, ChildFounders: 2100})
	if got := demeCount(pop, 1); got != 50 {
		t.Errorf("child deme size = %d, want 50 (half the parent)", got)
	}
	if got := demeCount(pop, 0); got != 50 {
		t.Errorf("parent deme size = %d, want 50", got)
	}
}

func TestDoSplitDeterministic(t *testing.T) {
	s := buildScenario(t)
	pick := func() map[int]bool {
		utils.SeedRNG(7)
		pop := makePop(200)
		s.doSplit(pop, &splitSpec{Parent: 0, Child: 1, ChildFounders: 300}) // 30 founders
		moved := map[int]bool{}
		for id, d := range pop.IndData {
			if d[individual.Deme] == 1 {
				moved[id] = true
			}
		}
		return moved
	}
	a, b := pick(), pick()
	if len(a) != 30 {
		t.Fatalf("moved %d, want 30", len(a))
	}
	for id := range a {
		if !b[id] {
			t.Fatalf("split not deterministic under a fixed seed")
		}
	}
}

func TestMigrateMovesBetweenDemes(t *testing.T) {
	utils.SeedRNG(1)
	s := buildScenario(t)
	// Two demes of 500 each, very high migration so movement is near-certain.
	pop := makePop(1000)
	for i := 500; i < 1000; i++ {
		pop.IndData[i][individual.Deme] = 1
	}
	e := &Epoch{Migration: []migSpec{{A: 0, B: 1, Rate: 1e9}}} // clamps to p=1
	before0 := demeCount(pop, 0)
	s.migrate(pop, e)
	after0 := demeCount(pop, 0)
	// With p=1 everyone swaps to the other deme, so counts trade places exactly.
	if after0 != 1000-before0 {
		t.Errorf("migrate with p=1: deme0 %d -> %d, expected full swap", before0, after0)
	}
}

func TestApplyFiresSplitAndSetsCaps(t *testing.T) {
	utils.SeedRNG(3)
	s := buildScenario(t)
	model := &core.Model{Parameters: map[string]float64{}, FreeParameters: map[string]int{}}
	pop := makePop(300)

	// Year 0 (ancestral): single deme, cap 730.
	model.FreeParameters["year"] = 0
	s.Apply(model, pop)
	if model.DemeCaps[0] != 730 || len(model.DemeCaps) != 1 {
		t.Errorf("ancestral caps = %v, want {0:730}", model.DemeCaps)
	}
	if demeCount(pop, 1) != 0 {
		t.Errorf("deme 1 should not exist before the OOA split")
	}

	// OOA split fires exactly at the epoch-3 start year (8500).
	model.FreeParameters["year"] = s.Epochs[2].startYear
	s.Apply(model, pop)
	if demeCount(pop, 1) == 0 {
		t.Errorf("OOA split did not create deme 1")
	}
	if _, ok := model.DemeCaps[1]; !ok {
		t.Errorf("caps missing deme 1 after OOA split: %v", model.DemeCaps)
	}
}

func TestLoadShippedOoA(t *testing.T) {
	path := filepath.Join("..", "..", "static", "basemodels", "OoA", "demography.json")
	s, err := Load(path)
	if err != nil {
		t.Fatalf("load shipped OoA: %v", err)
	}
	if len(s.Epochs) != 4 {
		t.Fatalf("epochs = %d, want 4", len(s.Epochs))
	}
	if s.Epochs[2].Split == nil || s.Epochs[2].Split.Child != 1 {
		t.Errorf("epoch 2 should split child 1")
	}
	if s.Epochs[3].Split == nil || s.Epochs[3].Split.Child != 2 {
		t.Errorf("epoch 3 should split child 2")
	}
	if s.TotalYears() != 22500 {
		t.Errorf("shipped OoA TotalYears = %d, want 22500", s.TotalYears())
	}
}
