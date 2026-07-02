package checkpoint

import (
	"drift/pkg/core"
	"drift/pkg/utils"
	"path/filepath"
	"strings"
	"testing"
)

func namedTestModel(t *testing.T) *core.Model {
	// Results dir under a temp dir; saves live alongside it.
	return &core.Model{
		ModelName:      "NamedTest",
		Scenario:       "default",
		ResultsDir:     filepath.Join(t.TempDir(), "results"),
		FreeParameters: map[string]int{"year": 100, "run": 1},
	}
}

func TestNamedSaveLoadList(t *testing.T) {
	utils.SeedRNG(1)
	model := namedTestModel(t)
	pop := &core.Pop{IndData: map[int][]int{1: {0}, 2: {1}, 3: {0}}}

	path, err := SaveNamed(model, pop, "baseline")
	if err != nil {
		t.Fatalf("SaveNamed: %v", err)
	}
	if !strings.HasSuffix(path, filepath.Join("saves", "baseline.ckpt")) {
		t.Errorf("unexpected save path: %s", path)
	}

	cp, err := LoadNamed(model, "baseline")
	if err != nil {
		t.Fatalf("LoadNamed: %v", err)
	}
	if cp.Year != 100 || len(cp.Pop.IndData) != 3 {
		t.Errorf("loaded state wrong: year=%d n=%d", cp.Year, len(cp.Pop.IndData))
	}

	saves, err := List(model)
	if err != nil {
		t.Fatalf("List: %v", err)
	}
	if len(saves) != 1 || saves[0].Name != "baseline" || saves[0].Year != 100 || saves[0].NumLiving != 3 {
		t.Errorf("List returned %+v", saves)
	}
}

func TestListEmptyWhenNoSavesDir(t *testing.T) {
	model := namedTestModel(t) // saves dir doesn't exist yet
	saves, err := List(model)
	if err != nil {
		t.Fatalf("List on missing dir should not error: %v", err)
	}
	if len(saves) != 0 {
		t.Errorf("expected no saves, got %v", saves)
	}
}

func TestSavePathIsContained(t *testing.T) {
	model := namedTestModel(t)
	// A traversal-y name must be reduced to its base component.
	p := SavePath(model, "../../etc/evil")
	if !strings.HasSuffix(p, filepath.Join("saves", "evil.ckpt")) {
		t.Errorf("SavePath did not contain name to saves dir: %s", p)
	}
}
