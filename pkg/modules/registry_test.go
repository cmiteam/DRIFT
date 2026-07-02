package modules

import (
	"drift/pkg/core"
	"testing"
)

// The real phase modules register from their owning packages (simulation,
// events, config). This test package registers its own probes so it does not
// depend on those imports; separate tests below cover resolution logic.

func TestRegisterAndAvailable(t *testing.T) {
	RegisterBirth("test_std", func(*core.Model, *core.Pop) {})
	defer delete(births, "test_std")

	found := false
	for _, n := range AvailableBirth() {
		if n == "test_std" {
			found = true
		}
	}
	if !found {
		t.Errorf("AvailableBirth() = %v, missing registered module", AvailableBirth())
	}
}

func TestDispatchSelectsModule(t *testing.T) {
	called := ""
	RegisterBirth("probe", func(*core.Model, *core.Pop) { called = "probe" })
	defer delete(births, "probe")

	model := &core.Model{StringParams: map[string]string{"birth_style": "probe"}}
	DispatchBirth(model, &core.Pop{})
	if called != "probe" {
		t.Errorf("DispatchBirth did not call the selected module (called=%q)", called)
	}
}

func TestLegacyIntegerResolution(t *testing.T) {
	// A model using the old integer mating_style must resolve to a name.
	model := &core.Model{Parameters: map[string]float64{"mating_style": 1}}
	if got := MatingStyle(model); got != "distance" {
		t.Errorf("MatingStyle(mating_style=1) = %q, want distance", got)
	}
	model2 := &core.Model{Parameters: map[string]float64{"seed_style": 2}}
	if got := SeedStyle(model2); got != "max_het" {
		t.Errorf("SeedStyle(seed_style=2) = %q, want max_het", got)
	}
	// String name takes precedence over any legacy int.
	model3 := &core.Model{
		Parameters:   map[string]float64{"mating_style": 1},
		StringParams: map[string]string{"mating_style": "random"},
	}
	if got := MatingStyle(model3); got != "random" {
		t.Errorf("string style should win over legacy int, got %q", got)
	}
}

func TestSetupStyleScenarioDriven(t *testing.T) {
	RegisterSetup("flood", func(*core.Model, *core.Pop) error { return nil })
	defer delete(setups, "flood")

	if got := SetupStyle(&core.Model{Scenario: "flood"}); got != "flood" {
		t.Errorf("SetupStyle(scenario=flood) = %q, want flood", got)
	}
	// Unknown scenario falls back to default.
	if got := SetupStyle(&core.Model{Scenario: "nonesuch"}); got != "default" {
		t.Errorf("SetupStyle(unknown scenario) = %q, want default", got)
	}
	// Explicit setup_style overrides scenario.
	m := &core.Model{Scenario: "flood", StringParams: map[string]string{"setup_style": "creation"}}
	if got := SetupStyle(m); got != "creation" {
		t.Errorf("setup_style override = %q, want creation", got)
	}
}

func TestValidateStyles(t *testing.T) {
	RegisterBirth(DefaultStyle, func(*core.Model, *core.Pop) {})
	RegisterDeath(DefaultStyle, func(*core.Model, *core.Pop) int { return 0 })
	RegisterMating("random", func(*core.Model, *core.Pop, []int, []int) {})
	RegisterSeed("single", func(*core.Model, *core.Pop) {})
	RegisterSetup("default", func(*core.Model, *core.Pop) error { return nil })
	defer func() {
		delete(births, DefaultStyle)
		delete(deaths, DefaultStyle)
		delete(matings, "random")
		delete(seeds, "single")
		delete(setups, "default")
	}()

	// Unset styles resolve to registered defaults → valid.
	if err := ValidateStyles(&core.Model{}); err != nil {
		t.Errorf("ValidateStyles(defaults) returned error: %v", err)
	}
	// Explicit unknown string style → error.
	bad := &core.Model{StringParams: map[string]string{"birth_style": "nope"}}
	if err := ValidateStyles(bad); err == nil {
		t.Errorf("ValidateStyles(unknown birth_style) returned nil, want error")
	}
}

func TestModelStringParam(t *testing.T) {
	m := &core.Model{StringParams: map[string]string{"birth_style": "bobs"}}
	if got := m.StringParam("birth_style", "standard"); got != "bobs" {
		t.Errorf("StringParam set = %q, want bobs", got)
	}
	if got := m.StringParam("death_style", "standard"); got != "standard" {
		t.Errorf("StringParam unset = %q, want standard", got)
	}
	var empty core.Model // nil map must not panic
	if got := empty.StringParam("x", "def"); got != "def" {
		t.Errorf("StringParam on nil map = %q, want def", got)
	}
}
