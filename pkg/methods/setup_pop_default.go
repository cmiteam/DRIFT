// File: pkg/methods/setup_pop_default.go
package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"sort"
)

// SetupPopDefault creates a population with age distribution
// This creates the initial population at the start of the simulation with
// individuals distributed across different ages according to the cumulative
// probability distribution defined in the model.
func SetupPopDefault(model *core.Model, pop *core.Pop) error {
	// Get habitable coordinates (empty if track_map == 0). With no habitat_suitability
	// table this is exactly the land cells; with one it also includes harsh terrains
	// given a positive suitability (roadmap §3).
	landCoordinates := utils.FindHabitableCells(model)

	// If no map tracking, use default location [0,0]
	if len(landCoordinates) == 0 {
		landCoordinates = [][2]int{{0, 0}}
	}

	popSize := int(model.Parameters["start_pop_size"])
	fitness := int(model.Parameters["mu_scale_factor"])
	lifespan := int(model.Parameters["lifespan"])

	// Number of demes (subpopulations). Default 1 = single-population behavior,
	// under which deme assignment is a no-op (every founder is deme 0). >1 enables
	// the island model consumed by Fst / f-stats (roadmap §6b).
	numDemes := int(model.Parameters["num_demes"])
	if numDemes < 1 {
		numDemes = 1
	}
	// When a real map exists, split its land cells into numDemes contiguous spatial
	// blocks so demes are also geographically separated (composes with spatial
	// mating and the §6c serial-founder thread). With no map there is a single
	// [0,0] cell and demes are purely abstract (isolation carried by mating).
	demeCells := partitionLandByDeme(landCoordinates, numDemes)

	// Ages sorted ascending so the inverse-CDF age draw is deterministic (map
	// iteration order is randomized) and always picks the smallest matching age.
	ages := make([]int, 0, len(model.CumulativeProb))
	for a := range model.CumulativeProb {
		ages = append(ages, a)
	}
	sort.Ints(ages)

	for i := 0; i < popSize; i++ {
		// Assign age based on cumulative probability distribution
		age := 0
		r := utils.RandFloat64()
		for _, a := range ages {
			if r <= model.CumulativeProb[a] {
				age = a
				break
			}
		}

		// Round-robin deme assignment keeps demes balanced and deterministic.
		deme := i % numDemes

		// Create individual with age-appropriate birth year, located within its
		// deme's spatial block.
		CreateFounder(pop, i, age, lifespan, fitness, model, demeCells[deme], deme)
	}

	model.FreeParameters["last_pop_size"] = len(pop.IndData)
	model.FreeParameters["indID"] = popSize

	return nil
}

// CreateFounder creates a founding individual for the initial population.
// deme is the founder's subpopulation label; landCoordinates are the land cells
// this founder may be placed in (its deme's spatial block when a map exists).
func CreateFounder(pop *core.Pop, id int, age int, lifespan int, fitness int, model *core.Model, landCoordinates [][2]int, deme int) {
	pop.IndData[id] = individual.MakeIndData()
	pop.IndData[id][individual.Deme] = deme
	pop.IndData[id][individual.BirthYear] = -age // Negative means born before simulation start
	pop.IndData[id][individual.Lifespan] = lifespan
	pop.IndData[id][individual.Sex] = utils.RandIntn(2) // 0 = male, 1 = female
	pop.IndData[id][individual.MarriageState] = -1      // Unmarried initially
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
	loc := landCoordinates[utils.RandIntn(len(landCoordinates))]
	pop.IndData[id][individual.Lat] = loc[0]
	pop.IndData[id][individual.Lon] = loc[1]
}

// partitionLandByDeme splits land cells into numDemes contiguous blocks so each
// deme's founders occupy a distinct geographic region. Cells are sorted (Lat, then
// Lon) for a deterministic, spatially-contiguous split. With a single cell (no map)
// or numDemes == 1, every deme shares the whole set — geographic separation is then
// a no-op and deme isolation is carried by partitioned mating instead.
func partitionLandByDeme(landCoordinates [][2]int, numDemes int) [][][2]int {
	blocks := make([][][2]int, numDemes)
	if numDemes <= 1 || len(landCoordinates) <= 1 {
		for d := 0; d < numDemes; d++ {
			blocks[d] = landCoordinates
		}
		return blocks
	}

	cells := make([][2]int, len(landCoordinates))
	copy(cells, landCoordinates)
	sort.Slice(cells, func(i, j int) bool {
		if cells[i][0] != cells[j][0] {
			return cells[i][0] < cells[j][0]
		}
		return cells[i][1] < cells[j][1]
	})

	// Even split; the last block absorbs any remainder. Every deme gets at least
	// one cell so no founder placement fails.
	per := len(cells) / numDemes
	if per < 1 {
		per = 1
	}
	for d := 0; d < numDemes; d++ {
		start := d * per
		end := start + per
		if d == numDemes-1 || end > len(cells) {
			end = len(cells)
		}
		if start >= len(cells) {
			// More demes than cells: wrap so every deme still has a cell.
			start = start % len(cells)
			end = start + 1
		}
		blocks[d] = cells[start:end]
	}
	return blocks
}
