package utils

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"drift/pkg/core"
)

// newBareModel returns a model with just the numeric map initialized, mirroring
// the empty models the loaders build before applying parameters.
func newBareModel() *core.Model {
	return &core.Model{Parameters: make(map[string]float64)}
}

// TestSetParameterRoutesStringParams is the core correctness fix: a genuinely-string
// parameter must land in StringParams even when its value is all-numeric, so
// StringParam() can read it. Previously strconv.ParseFloat success routed an
// all-numeric value into the numeric Parameters map, silently disabling the feature.
func TestSetParameterRoutesStringParams(t *testing.T) {
	cases := []struct{ name, value string }{
		{"movement_barriers", "3"},             // bare terrain code — the reported bug
		{"movement_barriers", "3:block"},       // non-numeric form that always worked
		{"migration_matrix", "0>1:0.5"},        // has ':' but leading token is numeric-ish
		{"mutation_classes", "point rate=8"},   // whitespace spec
		{"habitat_suitability", "1:1.0 5:0.3"}, // table
		{"geographic_fst_grid", "3"},           // bare-numeric grid, all-numeric
		{"fstats_demes", "0,1,2"},              // list
		{"linkage_model", "arm"},               // plain word
		{"deme_inheritance", "paternal"},       // plain word
	}
	for _, c := range cases {
		m := newBareModel()
		if err := SetParameter(m, c.name, c.value); err != nil {
			t.Fatalf("SetParameter(%s=%s) error: %v", c.name, c.value, err)
		}
		if got, ok := m.StringParams[c.name]; !ok || got != c.value {
			t.Errorf("%s=%s: StringParams[%q]=%q,ok=%v; want %q", c.name, c.value, c.name, got, ok, c.value)
		}
		if _, inNum := m.Parameters[c.name]; inNum {
			t.Errorf("%s=%s: leaked into numeric Parameters", c.name, c.value)
		}
		// StringParam() must now return it (the downstream read path).
		if got := m.StringParam(c.name, ""); got != c.value {
			t.Errorf("%s=%s: StringParam()=%q, want %q", c.name, c.value, got, c.value)
		}
	}
}

// TestSetParameterRoutesNumeric guards the byte-identical contract: the two legacy
// numeric-coded dropdowns and ordinary numeric params must STILL land in Parameters.
// mating_style=1 in Parameters is what resolves to distance-mating on the default run.
func TestSetParameterRoutesNumeric(t *testing.T) {
	cases := []struct {
		name  string
		value string
		want  float64
	}{
		{"mating_style", "1", 1},          // default distance mating — must stay numeric
		{"selection", "0", 0},             // legacy dead dropdown — stays numeric
		{"carrying_capacity", "250", 250}, // ordinary float param
		{"num_demes", "3", 3},
		{"mu", "10", 10},
	}
	for _, c := range cases {
		m := newBareModel()
		if err := SetParameter(m, c.name, c.value); err != nil {
			t.Fatalf("SetParameter(%s=%s) error: %v", c.name, c.value, err)
		}
		if got, ok := m.Parameters[c.name]; !ok || got != c.want {
			t.Errorf("%s=%s: Parameters[%q]=%v,ok=%v; want %v", c.name, c.value, c.name, got, ok, c.want)
		}
		if _, inStr := m.StringParams[c.name]; inStr {
			t.Errorf("%s=%s: should not be in StringParams", c.name, c.value)
		}
	}
}

// TestMovementBarriersEngagesAfterRouting is the distinguishable-outcome case: the
// exact scenario from the bug report — movement_barriers=3 set through the loader —
// now activates the barrier table end to end, whereas the old ParseFloat routing left
// it inert. Drives the real core accessor.
func TestMovementBarriersEngagesAfterRouting(t *testing.T) {
	m := newBareModel()
	if err := SetParameter(m, "movement_barriers", "3"); err != nil {
		t.Fatalf("SetParameter: %v", err)
	}
	if !m.MovementBarriersActive() {
		t.Fatal("movement_barriers=3 did not activate the barrier table (still inert)")
	}
	// Contrast: the old ParseFloat routing put "3" in Parameters, leaving StringParam
	// empty so the table parsed to inactive. Confirm it is now the string path.
	if m.StringParam("movement_barriers", "") != "3" {
		t.Errorf("movement_barriers not readable via StringParam after routing")
	}
}

// TestStringParamAllowlist keeps the hand-maintained stringParams set in sync with
// the authoritative schema: every format=string row in parameter_defaults.csv must be
// routed as a string, EXCEPT the two legacy numeric-coded dropdowns (mating_style,
// selection) and the specially-handled identifier rows. A new string param that is
// forgotten fails here loudly rather than shipping silently broken.
func TestStringParamAllowlist(t *testing.T) {
	path := filepath.Join("..", "..", "static", "parameter_defaults.csv")
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}

	// Rows handled specially in SetParameter/LoadParameters before type routing.
	specialCased := map[string]bool{
		"model_name": true, "username": true, "scenario": true, "map_name": true,
	}
	// Legacy dropdowns declared string but storing an integer code in Parameters.
	legacyNumericDropdowns := map[string]bool{"mating_style": true, "selection": true}

	lines := strings.Split(strings.ReplaceAll(string(data), "\r\n", "\n"), "\n")
	for i, line := range lines {
		if i == 0 || strings.TrimSpace(line) == "" {
			continue // header / blank
		}
		cols := strings.Split(line, ",")
		if len(cols) < 5 {
			continue
		}
		name := strings.TrimSpace(cols[0])
		format := strings.TrimSpace(cols[3])
		if format != "string" {
			// A numeric param must not be in the allowlist.
			if stringParams[name] {
				t.Errorf("%s is format=%s but is in stringParams", name, format)
			}
			continue
		}
		if specialCased[name] || legacyNumericDropdowns[name] {
			if stringParams[name] {
				t.Errorf("%s is a special/legacy row and must NOT be in stringParams", name)
			}
			continue
		}
		if !stringParams[name] {
			t.Errorf("%s is format=string in parameter_defaults.csv but missing from stringParams "+
				"(add it, or it will be silently mis-routed when set to an all-numeric value)", name)
		}
	}
}

// TestLoadParametersRoutesByType guards the second loader (csv.go LoadParameters, used
// by the base-model and legacy paths): after loading the real defaults, mating_style
// is numeric in Parameters and the empty-default string params live in StringParams.
func TestLoadParametersRoutesByType(t *testing.T) {
	m := newBareModel()
	if err := LoadParameters(m, filepath.Join("..", "..", "static")); err != nil {
		t.Fatalf("LoadParameters: %v", err)
	}
	if got, ok := m.Parameters["mating_style"]; !ok || got != 1 {
		t.Errorf("mating_style: Parameters=%v,ok=%v; want 1", got, ok)
	}
	if _, ok := m.Parameters["movement_barriers"]; ok {
		t.Errorf("movement_barriers leaked into numeric Parameters")
	}
	if _, ok := m.StringParams["movement_barriers"]; !ok {
		t.Errorf("movement_barriers (format=string) not in StringParams after LoadParameters")
	}
}
