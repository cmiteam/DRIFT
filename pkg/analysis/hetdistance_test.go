package analysis

import (
	"path/filepath"
	"testing"

	"drift/pkg/core"
)

// mockScheduler is a core.DemographyScheduler that also exposes DemeRanks, so the
// analysis picks up real serial-founder ranks (via the optional demeRanker
// interface) exactly as pkg/demography.Scenario does.
type mockScheduler struct{ ranks map[int]int }

func (m mockScheduler) Apply(*core.Model, *core.Pop) {}
func (m mockScheduler) DemeRanks() map[int]int       { return m.ranks }

// TestHetDistance_DecliningGradient: three demes at ranks 0,1,2 with sample
// derived-allele frequencies 0.5, 0.25, 0 at one shared polymorphic site give
// unbiased He of 2/3, 1/2, 0. The He~rank OLS slope is exactly -1/3 (a clean
// serial-founder decline) and He falls with rank.
func TestHetDistance_DecliningGradient(t *testing.T) {
	model := makeIBDTestModel()
	model.DemographyScheduler = mockScheduler{ranks: map[int]int{0: 0, 1: 1, 2: 2}}

	// Deme 0 (rank 0): both individuals heterozygous -> derived 2/4 = 0.5.
	d0a, c0a := fstInd(0, bit0, 0, 0, 0)
	d0b, c0b := fstInd(0, bit0, 0, 0, 0)
	// Deme 1 (rank 1): one heterozygous, one ancestral -> derived 1/4 = 0.25.
	d1a, c1a := fstInd(1, bit0, 0, 0, 0)
	d1b, c1b := fstInd(1, 0, 0, 0, 0)
	// Deme 2 (rank 2): both ancestral -> derived 0 (monomorphic within deme).
	d2a, c2a := fstInd(2, 0, 0, 0, 0)
	d2b, c2b := fstInd(2, 0, 0, 0, 0)

	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b}, {d2a, c2a}, {d2b, c2b},
	})

	res := ComputeHetDistance(model, pop, ids)
	if res.NumSites != 1 {
		t.Fatalf("NumSites = %d, want 1", res.NumSites)
	}
	if len(res.Demes) != 3 {
		t.Fatalf("expected 3 demes, got %d", len(res.Demes))
	}
	if !res.HasRanks {
		t.Errorf("HasRanks = false, want true (mock supplies ranks)")
	}
	// Demes are sorted by rank; verify He per deme.
	assertClose(t, "He[rank0]", res.Demes[0].ExpHet, 2.0/3.0)
	assertClose(t, "He[rank1]", res.Demes[1].ExpHet, 0.5)
	assertClose(t, "He[rank2]", res.Demes[2].ExpHet, 0.0)
	// Observed het: deme0 both het -> 1.0; deme1 one het of two -> 0.5; deme2 -> 0.
	assertClose(t, "Ho[rank0]", res.Demes[0].ObsHet, 1.0)
	assertClose(t, "Ho[rank1]", res.Demes[1].ObsHet, 0.5)
	assertClose(t, "Ho[rank2]", res.Demes[2].ObsHet, 0.0)
	// Gradient: slope exactly -1/3, declining.
	assertClose(t, "slope", res.Slope, -1.0/3.0)
	if res.Slope >= 0 {
		t.Errorf("slope = %v, want negative (heterozygosity declines with rank)", res.Slope)
	}
}

// TestHetDistance_NoSchedulerFallback: with no scheduler the deme label stands in
// as the rank and HasRanks is false, but the per-deme He is still computed.
func TestHetDistance_NoSchedulerFallback(t *testing.T) {
	model := makeIBDTestModel() // DemographyScheduler stays nil

	d0a, c0a := fstInd(0, bit0, 0, 0, 0)
	d0b, c0b := fstInd(0, bit0, 0, 0, 0)
	d1a, c1a := fstInd(1, 0, 0, 0, 0)
	d1b, c1b := fstInd(1, 0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b},
	})

	res := ComputeHetDistance(model, pop, ids)
	if res.HasRanks {
		t.Errorf("HasRanks = true, want false (no scheduler)")
	}
	if len(res.Demes) != 2 {
		t.Fatalf("expected 2 demes, got %d", len(res.Demes))
	}
	// Ranks default to the deme labels 0 and 1.
	if res.Demes[0].Rank != 0 || res.Demes[1].Rank != 1 {
		t.Errorf("fallback ranks = %d,%d, want 0,1", res.Demes[0].Rank, res.Demes[1].Rank)
	}
}

// TestHetDistance_SingleDemeGuard: with only one eligible deme the result is empty
// (a gradient is undefined).
func TestHetDistance_SingleDemeGuard(t *testing.T) {
	model := makeIBDTestModel()
	d0a, c0a := fstInd(0, bit0, 0, 0, 0)
	d0b, c0b := fstInd(0, 0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{{d0a, c0a}, {d0b, c0b}})

	res := ComputeHetDistance(model, pop, ids)
	if len(res.Demes) != 0 || res.NumSites != 0 {
		t.Errorf("single-deme result should be empty, got %d demes / %d sites", len(res.Demes), res.NumSites)
	}
}

// TestSaveHetDistance_WritesReadableCSVs: the per-deme table and the one-row
// gradient summary are both written without error.
func TestSaveHetDistance_WritesReadableCSVs(t *testing.T) {
	model := makeIBDTestModel()
	model.ResultsDir = t.TempDir()
	model.DemographyScheduler = mockScheduler{ranks: map[int]int{0: 0, 1: 1}}

	d0a, c0a := fstInd(0, bit0, 0, 0, 0)
	d0b, c0b := fstInd(0, bit0, 0, 0, 0)
	d1a, c1a := fstInd(1, 0, 0, 0, 0)
	d1b, c1b := fstInd(1, 0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{
		{d0a, c0a}, {d0b, c0b}, {d1a, c1a}, {d1b, c1b},
	})

	res := ComputeHetDistance(model, pop, ids)
	if err := SaveHetDistance(model, res); err != nil {
		t.Fatalf("SaveHetDistance failed: %v", err)
	}

	table := readCSV(t, filepath.Join(model.ResultsDir, "IBDTest_hetdistance_run1_year10.csv"))
	// header + one row per deme (2) = 3
	if len(table) != 3 {
		t.Errorf("hetdistance CSV: expected 3 rows, got %d: %v", len(table), table)
	}
	wantHeader := []string{"Deme", "Rank", "SampledDiploids", "ExpHet", "ObsHet", "NumSites"}
	if !equalRow(table[0], wantHeader) {
		t.Errorf("hetdistance header = %v, want %v", table[0], wantHeader)
	}

	grad := readCSV(t, filepath.Join(model.ResultsDir, "IBDTest_hetgradient_run1_year10.csv"))
	// header + one summary row = 2
	if len(grad) != 2 {
		t.Errorf("hetgradient CSV: expected 2 rows, got %d: %v", len(grad), grad)
	}
}
