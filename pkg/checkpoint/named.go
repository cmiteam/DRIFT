package checkpoint

import (
	"drift/pkg/core"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// ext is the file extension for named saved states.
const ext = ".ckpt"

// SavesDir returns the directory where named run states for a model live. Saved
// states belong to the model, so for a user model the location is canonical —
// users/<user>/models/<model>/saves — independent of where results happen to be
// written this run (ResultsDir is rewritten by SetupUserDirectories and can be
// overridden by -output-dir). For a base/legacy run with no user, it falls back
// to a "saves" folder alongside the results directory.
func SavesDir(model *core.Model) string {
	if model.Username != "" && model.ModelName != "" {
		return filepath.Join("users", model.Username, "models", model.ModelName, "saves")
	}
	return filepath.Join(filepath.Dir(model.ResultsDir), "saves")
}

// SavePath returns the on-disk path for a named save. The name is reduced to its
// base component so it cannot escape the saves directory.
func SavePath(model *core.Model, name string) string {
	return filepath.Join(SavesDir(model), filepath.Base(name)+ext)
}

// SaveNamed writes the current run state under a name in the model's saves
// directory (creating it if needed) and returns the path written.
func SaveNamed(model *core.Model, pop *core.Pop, name string) (string, error) {
	if strings.TrimSpace(name) == "" {
		return "", fmt.Errorf("save name is empty")
	}
	dir := SavesDir(model)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return "", fmt.Errorf("create saves dir %q: %w", dir, err)
	}
	path := SavePath(model, name)
	if err := Save(path, model, pop); err != nil {
		return "", err
	}
	return path, nil
}

// LoadNamed loads a named save from the model's saves directory.
func LoadNamed(model *core.Model, name string) (*Checkpoint, error) {
	return Load(SavePath(model, name))
}

// DeleteNamed removes a named save from the model's saves directory.
func DeleteNamed(model *core.Model, name string) error {
	if strings.TrimSpace(name) == "" {
		return fmt.Errorf("save name is empty")
	}
	return os.Remove(SavePath(model, name))
}

// ModelFor returns a minimal model identifying a user model, for saves
// operations (list/delete) that don't need the full model loaded.
func ModelFor(username, modelName string) *core.Model {
	return &core.Model{Username: username, ModelName: modelName}
}

// SaveInfo is lightweight metadata about a named save, for listing.
type SaveInfo struct {
	Name      string
	Year      int
	Run       int
	Scenario  string
	RNGSeed   int64
	SavedAt   string
	NumLiving int
}

// List returns metadata for every named save in the model's saves directory,
// sorted by name. Returns nil (no error) if the directory does not exist yet.
func List(model *core.Model) ([]SaveInfo, error) {
	dir := SavesDir(model)
	entries, err := os.ReadDir(dir)
	if os.IsNotExist(err) {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}

	var saves []SaveInfo
	for _, e := range entries {
		if e.IsDir() || !strings.HasSuffix(e.Name(), ext) {
			continue
		}
		name := strings.TrimSuffix(e.Name(), ext)
		cp, err := Load(filepath.Join(dir, e.Name()))
		if err != nil {
			// Skip unreadable / version-mismatched saves rather than failing the list.
			continue
		}
		living := 0
		if cp.Pop != nil {
			living = len(cp.Pop.IndData)
		}
		saves = append(saves, SaveInfo{
			Name:      name,
			Year:      cp.Year,
			Run:       cp.Run,
			Scenario:  cp.Scenario,
			RNGSeed:   cp.RNGSeed,
			SavedAt:   cp.SavedAt,
			NumLiving: living,
		})
	}
	sort.Slice(saves, func(i, j int) bool { return saves[i].Name < saves[j].Name })
	return saves, nil
}
