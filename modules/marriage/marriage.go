package marriage

import (
	"drift/types"
	"math"
	"math/rand"
	"time"
)

func Marriage(model *types.Model, pop *types.Pop) {
	var availableMen, availableWomen []int
	// Find eligible individuals
	for id, data := range pop.IndData {
		if data["marriage_state"] == -1 && model.FreeParameters["year"]-data["birth_year"] >= int(model.Parameters["maturity"]) {
			if data["sex"] == 0 {
				availableMen = append(availableMen, id)
			} else if data["sex"] == 1 {
				availableWomen = append(availableWomen, id)
			}
		}
	}

	// Get mating style
	if model.MatingStyle == "random" {
		randomMating(model, pop, availableMen, availableWomen)
	} else if model.MatingStyle == "distance" {
		distanceMating(model, pop, availableMen, availableWomen)
	}
}

// randomMating performs the original random mating algorithm
func randomMating(model *types.Model, pop *types.Pop, availableMen, availableWomen []int) {

	// Randomize people (note this will create unusual age gaps among married couples)
	rand.Seed(time.Now().UnixNano())
	rand.Shuffle(len(availableMen), func(i, j int) {
		availableMen[i], availableMen[j] = availableMen[j], availableMen[i]
	})
	rand.Shuffle(len(availableWomen), func(i, j int) {
		availableWomen[i], availableWomen[j] = availableWomen[j], availableWomen[i]
	})

	// Trim male and female lists to the shortest of the two
	if len(availableMen) > len(availableWomen) {
		availableMen = availableMen[:len(availableWomen)]
	} else {
		availableWomen = availableWomen[:len(availableMen)]
	}

	// Assign spouses
	for i := 0; i < len(availableMen); i++ {
		pop.IndData[availableMen[i]]["marriage_state"] = availableWomen[i]
		pop.IndData[availableWomen[i]]["marriage_state"] = availableMen[i]
		pop.Tracking["marriages"]++
	}
}

// distanceMating performs matching based on Manhattan distance
func distanceMating(model *types.Model, pop *types.Pop, availableMen, availableWomen []int) {
	maxDistance := int(model.Parameters["max_mating_distance"])

	// Debug information
	//	fmt.Printf("Distance mating: %d eligible men, %d eligible women, max distance: %d\n",
	//		len(availableMen), len(availableWomen), maxDistance)

	// Create matches based on distance constraints
	var matches [][2]int // Stores pairs of [womanID, manID]

	// For each woman, find eligible men within distance threshold
	for _, womanID := range availableWomen {
		womanLat := pop.IndData[womanID]["lat"]
		womanLon := pop.IndData[womanID]["lon"]

		// Debug: Count potential matches for this woman
		potentialMatches := 0

		// For each man, check if within range using Manhattan distance
		for _, manID := range availableMen {
			manLat := pop.IndData[manID]["lat"]
			manLon := pop.IndData[manID]["lon"]

			// Calculate Manhattan distance
			manhattanDist := abs(womanLat-manLat) + abs(womanLon-manLon)

			if manhattanDist <= maxDistance {
				// Within range, add as potential match
				matches = append(matches, [2]int{womanID, manID})
				potentialMatches++
			}
		}

		// If this woman had no matches, print her location
		if potentialMatches == 0 && len(availableMen) > 0 {
			//			fmt.Printf("Woman %d at (%d,%d) had no eligible matches within distance %d\n",
			//				womanID, womanLat, womanLon, maxDistance)

			// Print the closest man's distance
			closestDist := math.MaxInt32
			closestManID := -1

			for _, manID := range availableMen {
				manLat := pop.IndData[manID]["lat"]
				manLon := pop.IndData[manID]["lon"]

				latDiff := abs(womanLat - manLat)
				lonDiff := abs(womanLon - manLon)
				dist := latDiff + lonDiff

				if dist < closestDist {
					closestDist = dist
					closestManID = manID
				}
			}

			if closestManID >= 0 {
				//				fmt.Printf("  Closest man %d at (%d,%d) with distance %d\n",
				//					closestManID, pop.IndData[closestManID]["lat"],
				//					pop.IndData[closestManID]["lon"], closestDist)
			}
		}
	}

	//	fmt.Printf("Found %d potential matches within distance %d\n", len(matches), maxDistance)

	// Process matches (simple greedy algorithm for now)
	matched := make(map[int]bool) // Track who's already matched

	for _, match := range matches {
		womanID, manID := match[0], match[1]

		// If neither person is matched yet, create the marriage
		if !matched[womanID] && !matched[manID] {
			pop.IndData[womanID]["marriage_state"] = manID
			pop.IndData[manID]["marriage_state"] = womanID
			pop.Tracking["marriages"]++

			matched[womanID] = true
			matched[manID] = true
		}
	}

	// fmt.Printf("Created %d marriages based on distance constraints\n",
	//
	//	pop.Tracking["marriages"])
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
