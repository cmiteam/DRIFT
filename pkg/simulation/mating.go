package simulation

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/modules"
	"sort"
)

func Mating(model *core.Model, pop *core.Pop) {
	var availableMen, availableWomen []int

	// Find eligible individuals (unmarried and of mature age)
	for id, data := range pop.IndData {
		if data[individual.MarriageState] == -1 &&
			model.FreeParameters["year"]-data[individual.BirthYear] >= int(model.Parameters["maturity"]) {
			if data[individual.Sex] == 0 {
				availableMen = append(availableMen, id)
			} else if data[individual.Sex] == 1 {
				availableWomen = append(availableWomen, id)
			}
		}
	}

	// Sort so the lists have a deterministic order before the mating module
	// consumes randomness (Go map iteration order is randomized; without this a
	// fixed rng_seed would still produce different runs).
	sort.Ints(availableMen)
	sort.Ints(availableWomen)

	// Dispatch to the mating module selected by the mating_style parameter.
	modules.DispatchMating(model, pop, availableMen, availableWomen)
}

// createInfluenceGrid creates a grid where each cell contains IDs of women who can mate there
func createInfluenceGrid(model *core.Model, pop *core.Pop, availableWomen []int, maxDistance int) [][][]int {
	// Get map dimensions
	mapHeight := 101
	mapWidth := 101
	if model.Parameters["track_map"] == 1 {
		mapHeight = len(model.Map)
		mapWidth = len(model.Map[0])
	}
	// Create influence grid - each cell contains list of women who can mate there
	influenceGrid := make([][][]int, mapHeight)
	for i := range influenceGrid {
		influenceGrid[i] = make([][]int, mapWidth)
		for j := range influenceGrid[i] {
			influenceGrid[i][j] = make([]int, 0)
		}
	}

	// For each woman, add her ID to all cells within mating range
	for _, womanID := range availableWomen {
		womanLat := pop.IndData[womanID][individual.Lat]
		womanLon := pop.IndData[womanID][individual.Lon]

		// Add woman's ID to all cells within Manhattan distance
		for lat := womanLat - maxDistance; lat <= womanLat+maxDistance; lat++ {
			for lon := womanLon - maxDistance; lon <= womanLon+maxDistance; lon++ {
				// Check bounds
				if lat >= 0 && lat < mapHeight && lon >= 0 && lon < mapWidth {
					influenceGrid[lat][lon] = append(influenceGrid[lat][lon], womanID)
				}
			}
		}
	}

	return influenceGrid
}

func removeWomanFromGrid(influenceGrid [][][]int, womanID int, womanLat, womanLon, maxDistance int) {
	// Remove woman's ID from all cells within Manhattan distance
	for lat := womanLat - maxDistance; lat <= womanLat+maxDistance; lat++ {
		for lon := womanLon - maxDistance; lon <= womanLon+maxDistance; lon++ {
			// Check bounds
			if lat >= 0 && lat < len(influenceGrid) && lon >= 0 && lon < len(influenceGrid[0]) {
				influenceGrid[lat][lon] = removeFromSlice(influenceGrid[lat][lon], womanID)
			}
		}
	}
}

func removeFromSlice(slice []int, value int) []int {
	for i, v := range slice {
		if v == value {
			// Remove by swapping with last element and truncating
			slice[i] = slice[len(slice)-1]
			return slice[:len(slice)-1]
		}
	}
	return slice // Not found, return unchanged
}
