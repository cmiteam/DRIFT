package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
)

func MatingDistance(model *core.Model, pop *core.Pop, availableMen, availableWomen []int) {
	maxDistance := int(model.Parameters["max_mating_distance"])
	matched := make(map[int]bool) // Track who's already matched

	// Create influence grid for fast lookup
	influenceGrid := createInfluenceGrid(model, pop, availableWomen, maxDistance)

	// Shuffle for random processing order
	utils.RandShuffle(len(availableMen), func(i, j int) {
		availableMen[i], availableMen[j] = availableMen[j], availableMen[i]
	})

	// For each man, pick a random eligible woman from his cell
	for _, manID := range availableMen {
		if matched[manID] {
			continue
		}

		manLat := pop.IndData[manID][individual.Lat]
		manLon := pop.IndData[manID][individual.Lon]

		// Get all available women in this man's cell
		availableInCell := influenceGrid[manLat][manLon]

		// Filter out already matched women
		var eligibleWomen []int
		for _, womanID := range availableInCell {
			if !matched[womanID] {
				eligibleWomen = append(eligibleWomen, womanID)
			}
		}

		// Pick random woman if any available
		if len(eligibleWomen) > 0 {
			randomIndex := utils.RandIntn(len(eligibleWomen))
			bestWoman := eligibleWomen[randomIndex]

			// Create marriage
			pop.IndData[manID][individual.MarriageState] = bestWoman
			pop.IndData[bestWoman][individual.MarriageState] = manID

			// Woman moves to man's location (with terrain validation). Habitat-aware
			// (roadmap §3): == IsLand with no habitat table, but also admits harsh
			// habitable terrains when a suitability table is configured.
			if model.IsHabitable(manLat, manLon) {
				pop.IndData[bestWoman][individual.Lat] = manLat
				pop.IndData[bestWoman][individual.Lon] = manLon
			}

			matched[manID] = true
			matched[bestWoman] = true
			pop.Tracking["marriages"]++
		}
	}
}
