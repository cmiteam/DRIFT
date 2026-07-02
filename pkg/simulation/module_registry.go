package simulation

import (
	"drift/pkg/core"
	"drift/pkg/methods"
	"drift/pkg/modules"
)

// Register this package's phase modules into the central registry (pkg/modules).
// Birth/Death implementations live in this package (birthStandard/deathStandard);
// the mating algorithms live in pkg/methods.
func init() {
	modules.RegisterBirth(modules.DefaultStyle, birthStandard)
	modules.RegisterDeath(modules.DefaultStyle, deathStandard)

	modules.RegisterMating("random", methods.MatingRandom)
	modules.RegisterMating("distance", methods.MatingDistance)
	modules.RegisterMating("age_distance", methods.MatingAgeDistance)
}

// Birth dispatches to the birth module selected by the birth_style parameter.
func Birth(model *core.Model, pop *core.Pop) { modules.DispatchBirth(model, pop) }

// Death dispatches to the death module selected by the death_style parameter.
func Death(model *core.Model, pop *core.Pop) int { return modules.DispatchDeath(model, pop) }
