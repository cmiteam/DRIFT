// Package modules is the central registry of interchangeable simulation phase
// algorithms (birth, death, mating, seeding, setup).
//
// A user selects an algorithm per phase by name via a "<phase>_style" parameter
// (e.g. birth_style=standard, mating_style=distance); the Dispatch* functions
// route to the chosen implementation. Legacy models that set mating_style /
// seed_style as INTEGERS keep working — resolveStyle maps the old numeric codes
// to names.
//
// Discovery: Go cannot enumerate functions by name at runtime, so a naming
// convention alone can't populate a dropdown. Each module SELF-REGISTERS (from
// an init() in the package that owns its implementation — simulation, events,
// config) under a string key; the registry keys are the authoritative list of
// available modules (Available* / AllStyles feed the GUI dropdowns). To add an
// algorithm: write the function, register it in an init(), recompile — it then
// appears automatically.
//
// This package imports only core, so any package can register into it without an
// import cycle.
package modules

import (
	"drift/pkg/core"
	"fmt"
	"sort"
)

const DefaultStyle = "standard"

// Phase function signatures.
type (
	BirthFunc  func(model *core.Model, pop *core.Pop)
	DeathFunc  func(model *core.Model, pop *core.Pop) int
	MatingFunc func(model *core.Model, pop *core.Pop, availableMen, availableWomen []int)
	SeedFunc   func(model *core.Model, pop *core.Pop)
	SetupFunc  func(model *core.Model, pop *core.Pop) error
)

var (
	births  = map[string]BirthFunc{}
	deaths  = map[string]DeathFunc{}
	matings = map[string]MatingFunc{}
	seeds   = map[string]SeedFunc{}
	setups  = map[string]SetupFunc{}
)

// Registration — call from an init() in the package that owns the implementation.
func RegisterBirth(name string, f BirthFunc)   { births[name] = f }
func RegisterDeath(name string, f DeathFunc)   { deaths[name] = f }
func RegisterMating(name string, f MatingFunc) { matings[name] = f }
func RegisterSeed(name string, f SeedFunc)     { seeds[name] = f }
func RegisterSetup(name string, f SetupFunc)   { setups[name] = f }

// Available* return the registered module names (sorted) for each phase.
func AvailableBirth() []string  { return sortedKeys(births) }
func AvailableDeath() []string  { return sortedKeys(deaths) }
func AvailableMating() []string { return sortedKeys(matings) }
func AvailableSeed() []string   { return sortedKeys(seeds) }
func AvailableSetup() []string  { return sortedKeys(setups) }

// AllStyles returns the available module names per phase parameter — the data
// source for the GUI dropdowns.
func AllStyles() map[string][]string {
	return map[string][]string{
		"birth_style":  AvailableBirth(),
		"death_style":  AvailableDeath(),
		"mating_style": AvailableMating(),
		"seed_style":   AvailableSeed(),
		"setup_style":  AvailableSetup(),
	}
}

// Legacy integer codes (pre-registry) mapped to module names, for backward
// compatibility with existing models. See pkg/core/constants.go.
var (
	legacyMating = map[int]string{0: "random", 1: "distance", 2: "age_distance"}
	legacySeed   = map[int]string{0: "single", 1: "population", 2: "max_het"}
)

// resolveStyle picks the module name for a phase: the string param if set, else a
// legacy integer param mapped through legacyInts, else def.
func resolveStyle(model *core.Model, key string, legacyInts map[int]string, def string) string {
	if s := model.StringParam(key, ""); s != "" {
		return s
	}
	if legacyInts != nil {
		if v, ok := model.Parameters[key]; ok {
			if name, ok := legacyInts[int(v)]; ok {
				return name
			}
		}
	}
	return def
}

// Selected style name per phase.
func BirthStyle(model *core.Model) string  { return resolveStyle(model, "birth_style", nil, DefaultStyle) }
func DeathStyle(model *core.Model) string  { return resolveStyle(model, "death_style", nil, DefaultStyle) }
func MatingStyle(model *core.Model) string { return resolveStyle(model, "mating_style", legacyMating, "random") }
func SeedStyle(model *core.Model) string   { return resolveStyle(model, "seed_style", legacySeed, "single") }

// SetupStyle is scenario-driven by default: setup_style if set, else the model's
// scenario when a matching setup is registered, else "default".
func SetupStyle(model *core.Model) string {
	if s := model.StringParam("setup_style", ""); s != "" {
		return s
	}
	if _, ok := setups[model.Scenario]; ok {
		return model.Scenario
	}
	return "default"
}

// Dispatch* run the selected module for a phase, warning and falling back to a
// safe default if the chosen name is not registered.
func DispatchBirth(model *core.Model, pop *core.Pop) {
	if f, ok := births[BirthStyle(model)]; ok {
		f(model, pop)
		return
	}
	warnFallback("birth_style", BirthStyle(model), AvailableBirth(), DefaultStyle)
	births[DefaultStyle](model, pop)
}

func DispatchDeath(model *core.Model, pop *core.Pop) int {
	if f, ok := deaths[DeathStyle(model)]; ok {
		return f(model, pop)
	}
	warnFallback("death_style", DeathStyle(model), AvailableDeath(), DefaultStyle)
	return deaths[DefaultStyle](model, pop)
}

func DispatchMating(model *core.Model, pop *core.Pop, men, women []int) {
	if f, ok := matings[MatingStyle(model)]; ok {
		f(model, pop, men, women)
		return
	}
	warnFallback("mating_style", MatingStyle(model), AvailableMating(), "random")
	matings["random"](model, pop, men, women)
}

func DispatchSeed(model *core.Model, pop *core.Pop) {
	if f, ok := seeds[SeedStyle(model)]; ok {
		f(model, pop)
		return
	}
	warnFallback("seed_style", SeedStyle(model), AvailableSeed(), "single")
	seeds["single"](model, pop)
}

func DispatchSetup(model *core.Model, pop *core.Pop) error {
	if f, ok := setups[SetupStyle(model)]; ok {
		return f(model, pop)
	}
	warnFallback("setup_style", SetupStyle(model), AvailableSetup(), "default")
	return setups["default"](model, pop)
}

// ValidateStyles checks that every explicitly-selected phase module is
// registered, returning a clear error (listing what IS available) if not. Legacy
// integer codes and scenario-driven setup resolve safely and never error here;
// only an explicit unknown *_style name fails, so a GUI/param typo fails fast.
func ValidateStyles(model *core.Model) error {
	if s := BirthStyle(model); births[s] == nil {
		return unknown("birth_style", s, AvailableBirth())
	}
	if s := DeathStyle(model); deaths[s] == nil {
		return unknown("death_style", s, AvailableDeath())
	}
	if s := MatingStyle(model); matings[s] == nil {
		return unknown("mating_style", s, AvailableMating())
	}
	if s := SeedStyle(model); seeds[s] == nil {
		return unknown("seed_style", s, AvailableSeed())
	}
	if s := SetupStyle(model); setups[s] == nil {
		return unknown("setup_style", s, AvailableSetup())
	}
	return nil
}

func unknown(key, name string, avail []string) error {
	return fmt.Errorf("unknown %s %q; available: %v", key, name, avail)
}

func warnFallback(key, name string, avail []string, def string) {
	fmt.Printf("WARNING: unknown %s %q; available %v — using %q\n", key, name, avail, def)
}

func sortedKeys[V any](m map[string]V) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}
