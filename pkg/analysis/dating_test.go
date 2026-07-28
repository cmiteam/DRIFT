package analysis

import (
	"os"
	"path/filepath"
	"reflect"
	"testing"

	"drift/pkg/core"
	"drift/pkg/individual"
)

// buildDatingPop reuses the buildMutPop mutation pool (theta_pi = 1.25 exactly) but
// gives each individual a valid IndData slice (so FindYAdam/FindMtEve can read Sex)
// and empty genealogy DBs (so the true anchors resolve to "not found" / n/a).
func buildDatingPop() (*core.Pop, []int) {
	im := map[int]map[int][]int{
		1: {0: {101, 104}, 1: {102, 104}},
		2: {0: {102, 104}, 1: {104}},
		3: {0: {103, 104}, 1: {103, 104}},
		4: {0: {103, 104}, 1: {103, 104}},
	}
	mk := func(sex int) []int { d := individual.MakeIndData(); d[individual.Sex] = sex; return d }
	pop := &core.Pop{
		IndData:      map[int][]int{1: mk(0), 2: mk(0), 3: mk(1), 4: mk(1)},
		IndMutations: im,
		MutationPool: map[int]core.Mutation{
			101: {Id: 101}, 102: {Id: 102}, 103: {Id: 103}, 104: {Id: 104},
		},
		MaleDB:   map[int]core.Ancestor{},
		FemaleDB: map[int]core.Ancestor{},
	}
	return pop, []int{1, 2, 3, 4}
}

func datingModel(t *testing.T, mu, g float64, strs map[string]string) *core.Model {
	if strs == nil {
		strs = map[string]string{}
	}
	return &core.Model{
		Parameters:     map[string]float64{"mu": mu, "generation_time": g, "dating_sample_size": 0},
		StringParams:   strs,
		FreeParameters: map[string]int{"run": 1, "year": 100},
		ResultsDir:     t.TempDir(),
		ModelName:      "DateTest",
	}
}

// TestComputeDating_PointEstimate: with theta_pi = 1.25 and true mu = 0.5, the
// self-consistent pairwise coalescence time is 1.25/(2*0.5) = 1.25 generations, and
// at generation_time 25 the date is 31.25 years.
func TestComputeDating_PointEstimate(t *testing.T) {
	pop, ids := buildDatingPop()
	model := datingModel(t, 0.5, 25, nil)

	res := ComputeDating(model, pop, ids)
	if res == nil {
		t.Fatal("ComputeDating returned nil")
	}
	approx(t, "ThetaPi", res.ThetaPi, 1.25, 1e-9)
	approx(t, "PointTMRCAGen", res.PointTMRCAGen, 1.25, 1e-9)
	approx(t, "PointDateYears", res.PointDateYears, 31.25, 1e-9)

	// No genealogy DBs => true anchors are "not found".
	if res.YAdamYearsBack != -1 || res.MtEveYearsBack != -1 {
		t.Errorf("anchors = Y %d / Mt %d, want -1 / -1 (no genealogy)", res.YAdamYearsBack, res.MtEveYearsBack)
	}
}

// TestComputeDating_SensitivityGrid: the default grid (mu factors {0.5,1,2} x gen
// times {20,25,29,35}) gives 12 cells and a swing factor of exactly 7.0 — the oldest
// date (mf=0.5, g=35) over the youngest (mf=2, g=20) = (1/0.5*35)/(1/2*20) = 7 —
// independent of theta_pi and mu (the point of the mutation-rate-crisis deliverable).
func TestComputeDating_SensitivityGrid(t *testing.T) {
	pop, ids := buildDatingPop()
	model := datingModel(t, 0.5, 25, nil)

	res := ComputeDating(model, pop, ids)
	if len(res.Cells) != 12 {
		t.Fatalf("cells = %d, want 12 (3 mu x 4 g)", len(res.Cells))
	}
	approx(t, "SwingFactor", res.SwingFactor, 7.0, 1e-9)

	// Halving the assumed mu must exactly double the date at a fixed generation time.
	var dMf1G25, dMf2G25 float64
	for _, c := range res.Cells {
		if c.GenTime == 25 && c.MuFactor == 1 {
			dMf1G25 = c.DateYears
		}
		if c.GenTime == 25 && c.MuFactor == 2 {
			dMf2G25 = c.DateYears
		}
	}
	approx(t, "half-mu doubles date", dMf1G25, 2*dMf2G25, 1e-9)
}

// TestComputeDating_CustomGrid: an explicit single-point grid collapses the swing to
// 1.0 and yields one cell at the true mu / g.
func TestComputeDating_CustomGrid(t *testing.T) {
	pop, ids := buildDatingPop()
	model := datingModel(t, 0.5, 25, map[string]string{
		"dating_mu_factors": "1",
		"dating_gen_times":  "25",
	})

	res := ComputeDating(model, pop, ids)
	if len(res.Cells) != 1 {
		t.Fatalf("cells = %d, want 1", len(res.Cells))
	}
	approx(t, "SwingFactor", res.SwingFactor, 1.0, 1e-9)
	approx(t, "single cell date", res.Cells[0].DateYears, 31.25, 1e-9)
}

// TestComputeDating_NilWhenNoDivergence: an empty mutation pool has nothing to date.
func TestComputeDating_NilWhenNoDivergence(t *testing.T) {
	pop := &core.Pop{
		IndData:      map[int][]int{1: individual.MakeIndData(), 2: individual.MakeIndData()},
		IndMutations: map[int]map[int][]int{},
		MutationPool: map[int]core.Mutation{},
		MaleDB:       map[int]core.Ancestor{},
		FemaleDB:     map[int]core.Ancestor{},
	}
	model := datingModel(t, 0.5, 25, nil)
	if res := ComputeDating(model, pop, []int{1, 2}); res != nil {
		t.Errorf("expected nil result for no divergence, got %+v", res)
	}
}

// TestSaveDating_WritesFiles: the summary and sensitivity CSVs are both written; a
// nil result is a no-op.
func TestSaveDating_WritesFiles(t *testing.T) {
	pop, ids := buildDatingPop()
	model := datingModel(t, 0.5, 25, nil)
	res := ComputeDating(model, pop, ids)

	if err := SaveDating(model, res); err != nil {
		t.Fatalf("SaveDating: %v", err)
	}
	for _, suffix := range []string{"_dating_run1_year100.csv", "_dating_sensitivity_run1_year100.csv"} {
		p := filepath.Join(model.ResultsDir, model.ModelName+suffix)
		if _, err := os.Stat(p); err != nil {
			t.Errorf("expected output %s: %v", p, err)
		}
	}
	// Nil result is a clean no-op.
	if err := SaveDating(model, nil); err != nil {
		t.Errorf("SaveDating(nil) = %v, want nil", err)
	}
}

// TestParseFloatList covers separators, malformed/non-positive skipping, and the
// empty/all-invalid fallback to the default grid.
func TestParseFloatList(t *testing.T) {
	def := []float64{9, 9}
	cases := []struct {
		in   string
		want []float64
	}{
		{"0.5;1;2", []float64{0.5, 1, 2}},
		{"1, 2 3", []float64{1, 2, 3}}, // mixed comma/space separators
		{"2;bad;4", []float64{2, 4}},   // malformed token skipped
		{"", def},                      // empty => default
		{"abc;-1;0", def},              // all invalid/non-positive => default
	}
	for _, c := range cases {
		got := parseFloatList(c.in, def)
		if !reflect.DeepEqual(got, c.want) {
			t.Errorf("parseFloatList(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}
