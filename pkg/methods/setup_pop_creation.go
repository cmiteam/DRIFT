// File: pkg/methods/setup_pop_creation.go
package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
)

// SetupPopCreation creates 2 founding individuals (Adam and Eve)
// Both start at age 20 with extended lifespans
func SetupPopCreation(model *core.Model, pop *core.Pop) error {
	landCoordinates := utils.FindLandCells(model)
	print("Found", len(landCoordinates), "land coordinates\n")

	// If no map tracking, use default location [0,0]
	if len(landCoordinates) == 0 {
		landCoordinates = [][2]int{{0, 0}}
	}

	// Pick a random location for the couple
	loc := landCoordinates[utils.RandIntn(len(landCoordinates))]

	// Create Adam (male)
	adamID := 0
	CreateFounder(pop, adamID, 20, 900, 1, model, landCoordinates)
	pop.IndData[adamID][individual.Sex] = 0 // Male
	pop.IndData[adamID][individual.MarriageState] = 1 // Married to Eve
	pop.IndData[adamID][individual.Lat] = loc[0]
	pop.IndData[adamID][individual.Lon] = loc[1]
	pop.IndData[adamID][individual.BirthYear] = -20
	pop.IndData[adamID][individual.Lifespan] = 900

	// Create Eve (female)
	eveID := 1
	CreateFounder(pop, eveID, 20, 900, 1, model, landCoordinates)
	pop.IndData[eveID][individual.Sex] = 1 // Female
	pop.IndData[eveID][individual.MarriageState] = 0 // Married to Adam
	pop.IndData[eveID][individual.Lat] = loc[0]
	pop.IndData[eveID][individual.Lon] = loc[1]
	pop.IndData[eveID][individual.BirthYear] = -20
	pop.IndData[eveID][individual.Lifespan] = 900
	pop.IndData[eveID][individual.LastBirthYear] = -10 // Can have children immediately

	model.FreeParameters["last_pop_size"] = len(pop.IndData)
	model.FreeParameters["indID"] = 2

	return nil
}
