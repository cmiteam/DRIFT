// File: pkg/methods/setup_pop_flood.go
package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
)

// SetupPopFlood creates a small number of long-lived couples
// This scenario simulates starting with just a few founding couples who have
// extended lifespans, allowing them to build up a population over many generations.
// Default is 3 couples, but can be configured via parameters.
func SetupPopFlood(model *core.Model, pop *core.Pop) error {
	landCoordinates := utils.FindHabitableCells(model) // habitat-aware (roadmap §3); == land cells with no habitat table
	print("Found", len(landCoordinates), "land coordinates\n")

	// If no map tracking, use default location [0,0]
	if len(landCoordinates) == 0 {
		landCoordinates = [][2]int{{0, 0}}
	}

	// Get configuration from parameters
	numCouples := 3

	// Create the founding couples
	for couple := 0; couple < numCouples; couple++ {
		manID := couple * 2
		womanID := couple*2 + 1

		// Pick a random location for this couple
		loc := landCoordinates[utils.RandIntn(len(landCoordinates))]

		// Create the man
		CreateFounder(pop, manID, 100, 650, 1, model, landCoordinates, 0)
		pop.IndData[manID][individual.Sex] = 0 // Ensure male
		pop.IndData[manID][individual.MarriageState] = womanID
		pop.IndData[manID][individual.Lat] = loc[0]
		pop.IndData[manID][individual.Lon] = loc[1]
		pop.IndData[manID][individual.BirthYear] = -100
		pop.IndData[manID][individual.Lifespan] = 650

		// Create the woman
		CreateFounder(pop, womanID, 100, 650, 1, model, landCoordinates, 0)
		pop.IndData[womanID][individual.Sex] = 1 // Ensure female
		pop.IndData[womanID][individual.MarriageState] = manID
		pop.IndData[womanID][individual.Lat] = loc[0]
		pop.IndData[womanID][individual.Lon] = loc[1]
		pop.IndData[womanID][individual.BirthYear] = -100
		pop.IndData[womanID][individual.Lifespan] = 650
	}

	model.FreeParameters["last_pop_size"] = len(pop.IndData)
	model.FreeParameters["indID"] = numCouples * 2

	return nil
}
