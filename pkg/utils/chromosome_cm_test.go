package utils

import (
	"drift/pkg/core"
	"os"
	"path/filepath"
	"testing"
)

func writeTemp(t *testing.T, body string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "chromosome_data.csv")
	if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
		t.Fatalf("write temp: %v", err)
	}
	return p
}

func newModel() *core.Model {
	return &core.Model{
		Parameters:     map[string]float64{},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{},
	}
}

// TestChromosomeCMColumn: the optional 5th column (roadmap §6a) populates ArmCM and,
// via EnsureGeneticMap, yields a genetic map whose arm cM matches the file.
func TestChromosomeCMColumn(t *testing.T) {
	p := writeTemp(t, "Chromosome,Arm,Start,Length,cM\n1,0,0,100,40\n1,1,100,100,60\n")
	model := newModel()
	if err := LoadChromosomesFromPath(model, p); err != nil {
		t.Fatalf("load: %v", err)
	}
	if got := model.FreeParameters["genome_bits"]; got != 200 {
		t.Errorf("genome_bits = %d, want 200", got)
	}
	if model.ArmCM[1][0] != 40 || model.ArmCM[1][1] != 60 {
		t.Errorf("ArmCM = %v, want {0:40, 1:60}", model.ArmCM[1])
	}
	gm := EnsureGeneticMap(model)
	if gm == nil {
		t.Fatal("EnsureGeneticMap returned nil despite a cM column")
	}
	if gm.TotalCM(1) != 100 {
		t.Errorf("TotalCM = %v, want 100 (40+60)", gm.TotalCM(1))
	}
}

// TestChromosomeNoCMColumn: a legacy 4-column file leaves ArmCM empty and, with the
// default legacy recombination model, EnsureGeneticMap stays nil (byte-identical path).
func TestChromosomeNoCMColumn(t *testing.T) {
	p := writeTemp(t, "Chromosome,Arm,Start,Length\n1,0,0,100\n1,1,100,100\n")
	model := newModel()
	if err := LoadChromosomesFromPath(model, p); err != nil {
		t.Fatalf("load: %v", err)
	}
	if len(model.ArmCM) != 0 {
		t.Errorf("ArmCM = %v, want empty for a 4-column file", model.ArmCM)
	}
	if gm := EnsureGeneticMap(model); gm != nil {
		t.Error("EnsureGeneticMap should be nil for a legacy 4-column file in legacy mode")
	}
	// But requesting map mode builds a derived-rate map even without a cM column.
	model.StringParams["recombination_model"] = "map"
	model.GeneticMap = nil // clear the cached nil-decision from the prior call
	if gm := EnsureGeneticMap(model); gm == nil {
		t.Error("EnsureGeneticMap should build a derived map in map mode")
	} else if gm.TotalCM(1) != 200 { // 200 bits * default 1.0 cM/bit
		t.Errorf("derived TotalCM = %v, want 200", gm.TotalCM(1))
	}
}

// TestChromosomeMalformedCM: a non-numeric cM cell is a hard load error (fail-fast).
func TestChromosomeMalformedCM(t *testing.T) {
	p := writeTemp(t, "Chromosome,Arm,Start,Length,cM\n1,0,0,100,notanumber\n")
	model := newModel()
	if err := LoadChromosomesFromPath(model, p); err == nil {
		t.Error("expected an error for a malformed cM value, got nil")
	}
}
