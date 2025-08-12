package marriage

import (
	"drift/modules/individual"
	"drift/modules/maploader"
	"drift/modules/utils"
	"drift/types"
	"math"
	"math/rand"
	"time"
)

// Marriage handles pairing of eligible individuals in the population
func Marriage(model *types.Model, pop *types.Pop) {
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

	// Get mating style from parameters and execute directly
	if styleFloat, exists := model.Parameters["mating_style"]; exists {
		styleInt := int(styleFloat)
		switch styleInt {
		case 0:
			randomMating(model, pop, availableMen, availableWomen)
		case 1:
			distanceMating(model, pop, availableMen, availableWomen)
		case 2:
			ageDistanceMating(model, pop, availableMen, availableWomen)
		default:
			randomMating(model, pop, availableMen, availableWomen) // fallback
			print("Unknown mating style, defaulting to random\n")
		}
	} else {
		// Default if parameter doesn't exist
		ageDistanceMating(model, pop, availableMen, availableWomen)
		print("mating style not defined, using random mating (default)\n")
	}
}

// randomMating pairs individuals randomly (no distance constraints)
func randomMating(model *types.Model, pop *types.Pop, availableMen, availableWomen []int) {
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

		if maploader.IsLand(model, manLat, manLon) {
			pop.IndData[woman][individual.Lat] = manLat
			pop.IndData[woman][individual.Lon] = manLon
		}
		// If man's location is invalid, woman stays where she is

		pop.Tracking["marriages"]++
	}
}

// distanceMating pairs individuals based on proximity using influence grid
func distanceMating(model *types.Model, pop *types.Pop, availableMen, availableWomen []int) {
	maxDistance := int(model.Parameters["max_mating_distance"])
	matched := make(map[int]bool) // Track who's already matched

	// Create influence grid for fast lookup
	influenceGrid := createInfluenceGrid(model, pop, availableWomen, maxDistance)

	// Shuffle for random processing order
	rand.Seed(time.Now().UnixNano())
	rand.Shuffle(len(availableMen), func(i, j int) {
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
			randomIndex := rand.Intn(len(eligibleWomen))
			bestWoman := eligibleWomen[randomIndex]

			// Create marriage
			pop.IndData[manID][individual.MarriageState] = bestWoman
			pop.IndData[bestWoman][individual.MarriageState] = manID

			// Woman moves to man's location (with terrain validation)
			if maploader.IsLand(model, manLat, manLon) {
				pop.IndData[bestWoman][individual.Lat] = manLat
				pop.IndData[bestWoman][individual.Lon] = manLon
			}

			matched[manID] = true
			matched[bestWoman] = true
			pop.Tracking["marriages"]++
		}
	}
}

// ageDistanceMating pairs individuals based on age preference and distance
func ageDistanceMating(model *types.Model, pop *types.Pop, availableMen, availableWomen []int) {
	maxDistance := int(model.Parameters["max_mating_distance"])
	currentYear := model.FreeParameters["year"]

	// Create influence grid for fast lookup
	influenceGrid := createInfluenceGrid(model, pop, availableWomen, maxDistance)

	// Shuffle men for random processing order
	rand.Seed(time.Now().UnixNano())
	rand.Shuffle(len(availableMen), func(i, j int) {
		availableMen[i], availableMen[j] = availableMen[j], availableMen[i]
	})

	for _, manID := range availableMen {

		manLat := pop.IndData[manID][individual.Lat]
		manLon := pop.IndData[manID][individual.Lon]
		manAge := currentYear - pop.IndData[manID][individual.BirthYear]

		// Get all available women in this man's cell
		availableInCell := influenceGrid[manLat][manLon]

		// Find eligible candidates from women in this cell
		type candidate struct {
			id    int
			age   int
			score float64
		}

		var candidates []candidate

		for _, womanID := range availableInCell {

			womanAge := currentYear - pop.IndData[womanID][individual.BirthYear]

			// Calculate age preference score (higher = better)
			ageDiff := manAge - womanAge
			var ageScore float64
			if ageDiff >= 0 && ageDiff <= 5 {
				// Ideal: woman 0-5 years younger
				ageScore = 5.0 - float64(ageDiff)
			} else if ageDiff < 0 {
				// Woman older than man
				ageScore = math.Max(0.1, 1.0+float64(ageDiff)/2.0)
			} else {
				// Woman much younger
				ageScore = math.Max(0.1, 2.0-float64(ageDiff-5)/10.0)
			}

			candidates = append(candidates, candidate{
				id:    womanID,
				age:   womanAge,
				score: ageScore,
			})
		}

		// Select best candidate based on score
		if len(candidates) > 0 {
			// Use weighted random selection
			totalScore := 0.0
			for _, c := range candidates {
				totalScore += c.score
			}

			var bestWoman int = -1
			if totalScore > 0.001 {
				// Weighted random selection
				r := rand.Float64() * totalScore
				cumulative := 0.0
				for _, c := range candidates {
					cumulative += c.score
					if r <= cumulative {
						bestWoman = c.id
						break
					}
				}
			}

			// Fallback to first candidate if selection failed
			if bestWoman == -1 {
				bestWoman = candidates[0].id
			}

			// Create marriage
			pop.IndData[manID][individual.MarriageState] = bestWoman
			pop.IndData[bestWoman][individual.MarriageState] = manID

			// Remove bestWoman from all cells in influence grid
			removeWomanFromGrid(influenceGrid, bestWoman, pop.IndData[bestWoman][individual.Lat], pop.IndData[bestWoman][individual.Lon], maxDistance)

			// Woman moves to man's location
			pop.IndData[bestWoman][individual.Lat] = manLat
			pop.IndData[bestWoman][individual.Lon] = manLon

			pop.Tracking["marriages"]++
		} else {
			// No eligible women found - make man wander to find new territory
			newLat, newLon := utils.Wander(model, manLat, manLon)
			pop.IndData[manID][individual.Lat] = newLat
			pop.IndData[manID][individual.Lon] = newLon
		}
	}
}

// createInfluenceGrid creates a grid where each cell contains IDs of women who can mate there
func createInfluenceGrid(model *types.Model, pop *types.Pop, availableWomen []int, maxDistance int) [][][]int {
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
