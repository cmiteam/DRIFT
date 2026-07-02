package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"math"
)

func MatingAgeDistance(model *core.Model, pop *core.Pop, availableMen, availableWomen []int) {
	maxDistance := int(model.Parameters["max_mating_distance"])
	currentYear := model.FreeParameters["year"]

	// Create influence grid for fast lookup
	influenceGrid := createInfluenceGrid(model, pop, availableWomen, maxDistance)

	// Shuffle men for random processing order
	utils.RandShuffle(len(availableMen), func(i, j int) {
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
				r := utils.RandFloat64() * totalScore
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
