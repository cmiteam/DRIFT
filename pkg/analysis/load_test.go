package analysis

import (
	"encoding/csv"
	"math"
	"os"
	"path/filepath"
	"testing"

	"drift/pkg/core"
	"drift/pkg/individual"
)

// loadFixture builds a 3-individual population with a known mutation pool so the
// load statistics are hand-verifiable (additive fitness):
//
//	mut 10: Effect -0.1  homozygous in ALL 3  -> FIXED deleterious (copies 6=2N)
//	mut 20: Effect  0.0  homozygous in ALL 3  -> FIXED neutral
//	mut 30: Effect -0.05 het in ind1 only     -> segregating deleterious
//	mut 40: Effect +0.02 homozygous in ind2   -> segregating beneficial
//
// Fitness (additive, h=0.5): ind1 = 1 - 0.1 - 0.5*0.05 = 0.875;
// ind2 = 1 - 0.1 + 0.02 = 0.92; ind3 = 1 - 0.1 = 0.90. Mean = 0.898333.
func loadFixture() (*core.Model, *core.Pop) {
	model := &core.Model{
		Parameters:     map[string]float64{"track_mutations": 1, "mu": 1, "f_neutral": 0.9},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"year": 500, "run": 1},
		ModelName:      "LoadFix",
		ResultsDir:     ".",
	}
	ind := func() []int { return individual.MakeIndData() } // Deme defaults to 0
	pop := &core.Pop{
		IndData: map[int][]int{1: ind(), 2: ind(), 3: ind()},
		IndMutations: map[int]map[int][]int{
			1: {0: {10, 20, 30}, 1: {10, 20}},
			2: {0: {10, 20, 40}, 1: {10, 20, 40}},
			3: {0: {10, 20}, 1: {10, 20}},
		},
		MutationPool: map[int]core.Mutation{
			10: {Id: 10, Effect: -0.1, Dominance: 50, Class: 0},
			20: {Id: 20, Effect: 0.0, Dominance: 50, Class: 0},
			30: {Id: 30, Effect: -0.05, Dominance: 50, Class: 0},
			40: {Id: 40, Effect: 0.02, Dominance: 50, Class: 0},
		},
	}
	return model, pop
}

func loadApprox(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func TestComputeLoadSnapshot(t *testing.T) {
	model, pop := loadFixture()
	s := computeLoadSnapshot(model, pop, 500)

	if s.CensusN != 3 {
		t.Fatalf("CensusN=%d want 3", s.CensusN)
	}
	if !loadApprox(s.MeanFitness, 0.898333, 1e-5) {
		t.Errorf("MeanFitness=%.6f want 0.898333", s.MeanFitness)
	}
	if !loadApprox(s.TotalLoad, 0.101667, 1e-5) {
		t.Errorf("TotalLoad=%.6f want 0.101667", s.TotalLoad)
	}
	if !loadApprox(s.FixedLoad, 0.1, 1e-9) {
		t.Errorf("FixedLoad=%.6f want 0.1", s.FixedLoad)
	}
	if !loadApprox(s.SegLoad, 0.001667, 1e-5) {
		t.Errorf("SegLoad=%.6f want 0.001667", s.SegLoad)
	}
	if s.NumFixed != 2 || s.NumFixedDel != 1 || s.NumFixedBen != 0 {
		t.Errorf("fixed: num=%d del=%d ben=%d want 2/1/0", s.NumFixed, s.NumFixedDel, s.NumFixedBen)
	}
	if s.NumSeg != 2 {
		t.Errorf("NumSeg=%d want 2", s.NumSeg)
	}
	if !loadApprox(s.MeanMutPerInd, 5.0, 1e-9) {
		t.Errorf("MeanMutPerInd=%.4f want 5.0", s.MeanMutPerInd)
	}
	if !loadApprox(s.MeanDelPerInd, 7.0/3.0, 1e-9) {
		t.Errorf("MeanDelPerInd=%.4f want 2.3333", s.MeanDelPerInd)
	}
}

// TestComputeLoadSnapshotSynergistic: switching fitness_model changes the realized
// fitness/load (dominance+epistasis flow through CombineFitness), confirming the
// tracker reflects the active model rather than assuming additive.
func TestComputeLoadSnapshotSynergistic(t *testing.T) {
	model, pop := loadFixture()
	add := computeLoadSnapshot(model, pop, 500)
	model.StringParams["fitness_model"] = "synergistic"
	model.Parameters["epistasis_coefficient"] = 1
	syn := computeLoadSnapshot(model, pop, 500)
	if syn.TotalLoad <= add.TotalLoad {
		t.Errorf("synergistic load %.6f should exceed additive %.6f", syn.TotalLoad, add.TotalLoad)
	}
}

func TestComputeDFE(t *testing.T) {
	_, pop := loadFixture()
	seg, fixed, keys := computeDFE(pop)
	if len(keys) != 4 {
		t.Fatalf("got %d DFE buckets, want 4: %+v", len(keys), keys)
	}
	// mut10 -0.1 => del exp -1, fixed. mut20 0 => neutral, fixed.
	// mut30 -0.05 => del exp -2, seg. mut40 +0.02 => ben exp -2, seg.
	kDelFixed := dfeKey{0, -1, -1}
	kNeuFixed := dfeKey{0, 0, 0}
	kDelSeg := dfeKey{0, -1, -2}
	kBenSeg := dfeKey{0, 1, -2}
	// Selection's fingerprint: deleterious present in the segregating column, absent
	// from the fixed column.
	if fixed[kDelFixed] != 1 || seg[kDelFixed] != 0 {
		t.Errorf("del-fixed bucket: fixed=%d seg=%d want 1/0", fixed[kDelFixed], seg[kDelFixed])
	}
	if fixed[kNeuFixed] != 1 {
		t.Errorf("neutral-fixed count=%d want 1", fixed[kNeuFixed])
	}
	if seg[kDelSeg] != 1 || seg[kBenSeg] != 1 {
		t.Errorf("segregating del=%d ben=%d want 1/1", seg[kDelSeg], seg[kBenSeg])
	}
}

func TestMeltdownVerdict(t *testing.T) {
	rising := map[int]*core.LoadSnapshot{
		100: {TotalLoad: 0.01}, 200: {TotalLoad: 0.02}, 300: {TotalLoad: 0.03},
		400: {TotalLoad: 0.05}, 500: {TotalLoad: 0.08}, 600: {TotalLoad: 0.10},
	}
	if v := meltdownVerdict(rising); v != "rising (load accumulating)" {
		t.Errorf("rising history verdict=%q", v)
	}
	stationary := map[int]*core.LoadSnapshot{
		100: {TotalLoad: 0.05}, 200: {TotalLoad: 0.051}, 300: {TotalLoad: 0.049},
		400: {TotalLoad: 0.050}, 500: {TotalLoad: 0.0505}, 600: {TotalLoad: 0.0495},
	}
	if v := meltdownVerdict(stationary); v != "stationary (quasi-balance)" {
		t.Errorf("stationary history verdict=%q", v)
	}
	if v := meltdownVerdict(map[int]*core.LoadSnapshot{1: {}}); v != "insufficient-windows" {
		t.Errorf("short history verdict=%q", v)
	}
}

func TestSaveLoadTimeSeriesWritesFiles(t *testing.T) {
	model, pop := loadFixture()
	dir := t.TempDir()
	model.ResultsDir = dir
	pop.LoadHistory = map[int]*core.LoadSnapshot{
		100: computeLoadSnapshot(model, pop, 100),
		300: computeLoadSnapshot(model, pop, 300),
		500: computeLoadSnapshot(model, pop, 500),
	}
	if err := SaveLoadTimeSeries(model, pop); err != nil {
		t.Fatalf("SaveLoadTimeSeries: %v", err)
	}
	for _, name := range []string{"LoadFix_load_timeseries_run1.csv", "LoadFix_dfe_run1.csv", "LoadFix_load_summary_run1.csv"} {
		p := filepath.Join(dir, name)
		info, err := os.Stat(p)
		if err != nil {
			t.Errorf("expected %s: %v", name, err)
			continue
		}
		if info.Size() == 0 {
			t.Errorf("%s is empty", name)
		}
	}
	// The time series must have a header + 3 data rows (years 100/300/500).
	f, _ := os.Open(filepath.Join(dir, "LoadFix_load_timeseries_run1.csv"))
	defer f.Close()
	rows, _ := csv.NewReader(f).ReadAll()
	if len(rows) != 4 {
		t.Errorf("time series rows=%d want 4 (header+3)", len(rows))
	}
}

// TestSaveLoadTimeSeriesNoMutations: track_mutations off => graceful no-op (no
// files, no error).
func TestSaveLoadTimeSeriesNoMutations(t *testing.T) {
	model, pop := loadFixture()
	model.Parameters["track_mutations"] = 0
	dir := t.TempDir()
	model.ResultsDir = dir
	if err := SaveLoadTimeSeries(model, pop); err != nil {
		t.Fatalf("SaveLoadTimeSeries no-op: %v", err)
	}
	entries, _ := os.ReadDir(dir)
	if len(entries) != 0 {
		t.Errorf("expected no files when track_mutations off, got %d", len(entries))
	}
}
