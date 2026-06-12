package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"math/rand"
	"time"
)

func MatingRandom(model *core.Model, pop *core.Pop, availableMen, availableWomen []int) {

	// Randomize people
	rand.Seed(time.Now().UnixNano())
	rand.Shuffle(len(availableMen), func(i, j int) {
		availableMen[i], availableMen[j] = availableMen[j], availableMen[i]
	})
	rand.Shuffle(len(availableWomen), func(i, j int) {
		availableWomen[i], availableWomen[j] = availableWomen[j], availableWomen[i]
	})

	// Pair up to the shortest list
	pairCount := utils.Min(len(availableMen), len(availableWomen))

	for i := 0; i < pairCount; i++ {
		man := availableMen[i]
		woman := availableWomen[i]

		// Create marriage
		pop.IndData[man][individual.MarriageState] = woman
		pop.IndData[woman][individual.MarriageState] = man

		// Woman moves to man's location
		manLat := pop.IndData[man][individual.Lat]
		manLon := pop.IndData[man][individual.Lon]

		if utils.IsLand(model, manLat, manLon) {
			pop.IndData[woman][individual.Lat] = manLat
			pop.IndData[woman][individual.Lon] = manLon
		}
		// If man's location is invalid, woman stays where she is

		pop.Tracking["marriages"]++
	}
}
