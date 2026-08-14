// File: pkg/methods/setup_pop_eden.go
package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
)

// SetupPopEden creates a SINGLE long-lived founding couple — the created-origin scenario
// TMR4A.md W4 is built around, and the one the TMR4A claim is actually about.
//
// It is the flood scenario with one couple instead of three, which matters more than it
// sounds: with one couple the founder generation is exactly four somatic strands, so
// "four alleles per locus" is a real ceiling under ordinary diploidy and the only way past
// it is created germline heterogeneity (founder_allele_model="created", A > 4). That makes
// this scenario the clean test bed — with three couples the ceiling is twelve and the
// question TMR4A poses does not arise.
//
// Pair it with seed_style="created" so the couple's germline carries its created-allele
// pool; with any other seeder the couple is an ordinary diploid pair and every locus is
// capped at four founder haplotypes by construction.
//
// Consumes RNG only for the couple's placement, exactly as SetupPopFlood does.
func SetupPopEden(model *core.Model, pop *core.Pop) error {
	landCoordinates := utils.FindHabitableCells(model) // habitat-aware (§3); == land cells with no habitat table
	if len(landCoordinates) == 0 {
		landCoordinates = [][2]int{{0, 0}}
	}

	// Long lifespans so one couple can found a population before dying, matching the
	// flood template. Ages/lifespans are deliberately the same numbers as
	// SetupPopFlood's so the two scenarios differ only in the number of couples.
	const (
		founderAge      = 100
		founderLifespan = 650
		founderBirth    = -100
	)

	loc := landCoordinates[utils.RandIntn(len(landCoordinates))]

	const manID, womanID = 0, 1

	CreateFounder(pop, manID, founderAge, founderLifespan, 1, model, landCoordinates, 0)
	pop.IndData[manID][individual.Sex] = 0
	pop.IndData[manID][individual.MarriageState] = womanID
	pop.IndData[manID][individual.Lat] = loc[0]
	pop.IndData[manID][individual.Lon] = loc[1]
	pop.IndData[manID][individual.BirthYear] = founderBirth
	pop.IndData[manID][individual.Lifespan] = founderLifespan

	CreateFounder(pop, womanID, founderAge, founderLifespan, 1, model, landCoordinates, 0)
	pop.IndData[womanID][individual.Sex] = 1
	pop.IndData[womanID][individual.MarriageState] = manID
	pop.IndData[womanID][individual.Lat] = loc[0]
	pop.IndData[womanID][individual.Lon] = loc[1]
	pop.IndData[womanID][individual.BirthYear] = founderBirth
	pop.IndData[womanID][individual.Lifespan] = founderLifespan

	model.FreeParameters["last_pop_size"] = len(pop.IndData)
	model.FreeParameters["indID"] = 2

	return nil
}
