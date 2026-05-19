// File: pkg/methods/setup_pop_default.go
package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"math/rand"
	"time"
)

// SetupPopDefault creates a population with age distribution
// This creates the initial population at the start of the simulation with
// individuals distributed across different ages according to the cumulative
// probability distribution defined in the model.
func SetupPopDefault(model *core.Model, pop *core.Pop) error {
	// Initialize random number generator
	rand.Seed(time.Now().UnixNano())

	// Get land coordinates (empty if track_map == 0)
	landCoordinates := utils.FindLandCells(model)

	// If no map tracking, use default location [0,0]
	if len(landCoordinates) == 0 {
		landCoordinates = [][2]int{{0, 0}}
	}

	popSize := int(model.Parameters["start_pop_size"])
	fitness := int(model.Parameters["mu_scale_factor"])
	lifespan := int(model.Parameters["lifespan"])

	for i := 0; i < popSize; i++ {
		// Assign age based on cumulative probability distribution
		age := 0
		r := rand.Float64()
		for a, cumProb := range model.CumulativeProb {
			if r <= cumProb {
				age = a
				break
			}
		}

		// Create individual with age-appropriate birth year
		CreateFounder(pop, i, age, lifespan, fitness, model, landCoordinates)
	}

	model.FreeParameters["last_pop_size"] = len(pop.IndData)
	model.FreeParameters["indID"] = popSize

	return nil
}

// CreateFounder creates a founding individual for the initial population
func CreateFounder(pop *core.Pop, id int, age int, lifespan int, fitness int, model *core.Model, landCoordinates [][2]int) {
	pop.IndData[id] = individual.MakeIndData()
	pop.IndData[id][individual.BirthYear] = -age // Negative means born before simulation start
	pop.IndData[id][individual.Lifespan] = lifespan
	pop.IndData[id][individual.Sex] = rand.Intn(2) // 0 = male, 1 = female
	pop.IndData[id][individual.MarriageState] = -1 // Unmarried initially
	pop.IndData[id][individual.NumBirths] = 0
	pop.IndData[id][individual.LastBirthYear] = 0
	pop.IndData[id][individual.Fitness] = fitness
	pop.IndData[id][individual.AlleleCount] = 0
	pop.IndData[id][individual.NumCentromeres] = 0
	pop.IndData[id][individual.YGens] = -1
	pop.IndData[id][individual.MtGens] = -1
	pop.IndData[id][individual.MinGenealoGens] = -1
	pop.IndData[id][individual.MaxGenealoGens] = -1
	pop.IndData[id][individual.Sons] = 0
	pop.IndData[id][individual.Daughters] = 0

	// Assign random land location
	loc := landCoordinates[rand.Intn(len(landCoordinates))]
	pop.IndData[id][individual.Lat] = loc[0]
	pop.IndData[id][individual.Lon] = loc[1]
}
