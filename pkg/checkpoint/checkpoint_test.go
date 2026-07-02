package checkpoint

import (
	"drift/pkg/core"
	"drift/pkg/utils"
	"path/filepath"
	"testing"
)

func sampleModelPop() (*core.Model, *core.Pop) {
	model := &core.Model{
		ModelName:      "CkptTest",
		Scenario:       "default",
		FreeParameters: map[string]int{"year": 42, "run": 1, "indID": 99, "mutID": 7, "seed": -1},
	}
	pop := &core.Pop{
		IndData:      map[int][]int{1: {0, -1, -1, 5}, 2: {1, -1, -1, 6}},
		Chromosomes:  map[int][][]uint64{1: {{0xFF, 0x00}, {0x0F, 0x00}}},
		MutationPool: map[int]core.Mutation{7: {Id: 7, Position: 100, Effect: -0.01, Count: 3}},
		Tracking:     map[string]int{"births": 12, "marriages": 4},
		MaleDB:       map[int]core.Ancestor{1: {ID: 0, BirthYear: -20}},
	}
	return model, pop
}

func TestSaveLoadRoundTrip(t *testing.T) {
	utils.SeedRNG(2024)
	model, pop := sampleModelPop()
	path := filepath.Join(t.TempDir(), "run.ckpt")

	if err := Save(path, model, pop); err != nil {
		t.Fatalf("Save: %v", err)
	}
	cp, err := Load(path)
	if err != nil {
		t.Fatalf("Load: %v", err)
	}

	if cp.Year != 42 || cp.Run != 1 {
		t.Errorf("Year/Run = %d/%d, want 42/1", cp.Year, cp.Run)
	}
	if cp.ModelName != "CkptTest" || cp.Scenario != "default" {
		t.Errorf("model metadata not preserved: %q/%q", cp.ModelName, cp.Scenario)
	}
	if cp.FreeParameters["indID"] != 99 || cp.FreeParameters["mutID"] != 7 {
		t.Errorf("FreeParameters not preserved: %v", cp.FreeParameters)
	}
	if len(cp.Pop.IndData) != 2 || cp.Pop.IndData[2][3] != 6 {
		t.Errorf("Pop.IndData not preserved: %v", cp.Pop.IndData)
	}
	if cp.Pop.Chromosomes[1][0][0] != 0xFF {
		t.Errorf("Pop.Chromosomes not preserved: %v", cp.Pop.Chromosomes)
	}
	if m := cp.Pop.MutationPool[7]; m.Count != 3 || m.Position != 100 {
		t.Errorf("Pop.MutationPool not preserved: %+v", m)
	}
	if cp.Pop.Tracking["births"] != 12 {
		t.Errorf("Pop.Tracking not preserved: %v", cp.Pop.Tracking)
	}
	if len(cp.RNGState) != 32 {
		t.Errorf("RNGState length = %d, want 32", len(cp.RNGState))
	}
}

// Restoring a checkpoint must reproduce the identical downstream RNG stream.
func TestRestoreReproducesStream(t *testing.T) {
	utils.SeedRNG(777)
	for i := 0; i < 5; i++ {
		utils.RandFloat64()
	}
	model, pop := sampleModelPop()
	path := filepath.Join(t.TempDir(), "run.ckpt")
	if err := Save(path, model, pop); err != nil {
		t.Fatalf("Save: %v", err)
	}
	want := []float64{utils.RandFloat64(), utils.RandFloat64()}

	cp, err := Load(path)
	if err != nil {
		t.Fatalf("Load: %v", err)
	}
	fresh := &core.Model{FreeParameters: map[string]int{}}
	if err := Restore(fresh, cp, 0); err != nil { // 0 = exact resume
		t.Fatalf("Restore: %v", err)
	}
	got := []float64{utils.RandFloat64(), utils.RandFloat64()}
	if got[0] != want[0] || got[1] != want[1] {
		t.Errorf("resume stream = %v, want %v", got, want)
	}
	if fresh.FreeParameters["indID"] != 99 {
		t.Errorf("Restore did not overlay FreeParameters: %v", fresh.FreeParameters)
	}
}

// A fork (non-zero forkSeed) must diverge from the exact-resume stream.
func TestForkDiverges(t *testing.T) {
	utils.SeedRNG(555)
	model, pop := sampleModelPop()
	path := filepath.Join(t.TempDir(), "run.ckpt")
	if err := Save(path, model, pop); err != nil {
		t.Fatalf("Save: %v", err)
	}

	cp, _ := Load(path)
	fresh := &core.Model{FreeParameters: map[string]int{}}
	Restore(fresh, cp, 0)
	exact := utils.RandFloat64()

	cp2, _ := Load(path)
	Restore(fresh, cp2, 12345) // fork
	forked := utils.RandFloat64()

	if exact == forked {
		t.Errorf("fork did not diverge from exact resume (both %v)", exact)
	}
}
