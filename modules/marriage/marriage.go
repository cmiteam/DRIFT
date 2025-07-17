package marriage

import (
	"drift/modules/individual"
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

	// Get mating style from parameters (default to 0 = random)
	matingStyle := 0
	if style, exists := model.Parameters["mating_style"]; exists {
		matingStyle = int(style)
	}
	switch matingStyle {
	case 0:
		randomMating(model, pop, availableMen, availableWomen)
	case 1:
		distanceMating(model, pop, availableMen, availableWomen)
	case 2:
		ageDistanceMating(model, pop, availableMen, availableWomen)
	default:
		randomMating(model, pop, availableMen, availableWomen) // fallback
	}
}

// randomMating pairs individuals randomly
func randomMating(model *types.Model, pop *types.Pop, availableMen, availableWomen []int) {
	// Randomize people
	rand.Seed(time.Now().UnixNano())
	rand.Shuffle(len(availableMen), func(i, j int) {
		availableMen[i], availableMen[j] = availableMen[j], availableMen[i]
	})
	rand.Shuffle(len(availableWomen), func(i, j int) {
		availableWomen[i], availableWomen[j] = availableWomen[j], availableWomen[i]
	})

	// Trim to the shortest of the two lists
	pairCount := min(len(availableMen), len(availableWomen))

	// Assign spouses
	for i := 0; i < pairCount; i++ {
		man := availableMen[i]
		woman := availableWomen[i]
		pop.IndData[man][individual.MarriageState] = woman
		pop.IndData[woman][individual.MarriageState] = man
		pop.IndData[woman][individual.Lat] = pop.IndData[man][individual.Lat]
		pop.IndData[woman][individual.Lon] = pop.IndData[man][individual.Lon]
		pop.Tracking["marriages"]++
	}
}

// distanceMating pairs individuals based on proximity
func distanceMating(model *types.Model, pop *types.Pop, availableMen, availableWomen []int) {
	maxDistance := int(model.Parameters["max_mating_distance"])

	// Create matches based on distance constraints
	var matches [][2]int // Stores pairs of [womanID, manID]

	// For each woman, find eligible men within distance threshold
	for _, womanID := range availableWomen {
		womanLat := pop.IndData[womanID][individual.Lat]
		womanLon := pop.IndData[womanID][individual.Lon]

		// Debug: Count potential matches for this woman
		potentialMatches := 0

		// For each man, check if within range using Manhattan distance
		for _, manID := range availableMen {
			manLat := pop.IndData[manID][individual.Lat]
			manLon := pop.IndData[manID][individual.Lon]

			// Calculate Manhattan distance
			manhattanDist := abs(womanLat-manLat) + abs(womanLon-manLon)

			if manhattanDist <= maxDistance {
				// Within range, add as potential match
				matches = append(matches, [2]int{womanID, manID})
				potentialMatches++
			}
		}

		if potentialMatches == 0 && len(availableMen) > 0 {
			closestDist := math.MaxInt32
			//closestManID := -1

			for _, manID := range availableMen {
				manLat := pop.IndData[manID][individual.Lat]
				manLon := pop.IndData[manID][individual.Lon]

				latDiff := abs(womanLat - manLat)
				lonDiff := abs(womanLon - manLon)
				dist := latDiff + lonDiff

				if dist < closestDist {
					closestDist = dist
					//closestManID = manID
				}
			}
		}
	}

	// Process matches (simple greedy algorithm for now)
	matched := make(map[int]bool) // Track who's already matched

	for _, match := range matches {
		womanID, manID := match[0], match[1]

		// If neither person is matched yet, create the marriage
		if !matched[womanID] && !matched[manID] {
			pop.IndData[womanID][individual.MarriageState] = manID
			pop.IndData[manID][individual.MarriageState] = womanID
			pop.Tracking["marriages"]++

			matched[womanID] = true
			matched[manID] = true
		}
	}

	// fmt.Printf("Created %d marriages based on distance constraints\n",
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// ageDistanceMatingOptimized pairs individuals based on distance and age preference with optimizations
func ageDistanceMating(model *types.Model, pop *types.Pop, availableMen, availableWomen []int) {
	maxDistance := int(model.Parameters["max_mating_distance"])
	currentYear := model.FreeParameters["year"]

	// Pre-calculate all ages to avoid repeated calculations
	ages := make(map[int]int, len(pop.IndData))
	for id, data := range pop.IndData {
		ages[id] = currentYear - data[individual.BirthYear]
	}

	// Build spatial index for women
	spatialIndex := make(map[spatialKey][]int)
	cellSize := maxDistance // Cell size equal to max mating distance

	for _, womanID := range availableWomen {
		lat := pop.IndData[womanID][individual.Lat]
		lon := pop.IndData[womanID][individual.Lon]
		key := spatialKey{lat / cellSize, lon / cellSize}
		spatialIndex[key] = append(spatialIndex[key], womanID)
	}

	// Shuffle men for random access order
	rand.Seed(time.Now().UnixNano())
	rand.Shuffle(len(availableMen), func(i, j int) {
		availableMen[i], availableMen[j] = availableMen[j], availableMen[i]
	})

	// Track matched women
	matchedWomen := make(map[int]bool, len(availableWomen))

	// For each man, find a suitable woman
	for _, manID := range availableMen {
		manLat := pop.IndData[manID][individual.Lat]
		manLon := pop.IndData[manID][individual.Lon]
		manAge := ages[manID]

		// Get women in nearby cells
		nearbyWomen := getNearbyWomenFromSpatialIndex(spatialIndex, manLat, manLon, cellSize, maxDistance)

		// Filter for unmatched women and calculate actual distances
		var eligibleWomen []struct {
			id       int
			distance int
		}

		for _, womanID := range nearbyWomen {
			if !matchedWomen[womanID] {
				womanLat := pop.IndData[womanID][individual.Lat]
				womanLon := pop.IndData[womanID][individual.Lon]

				// Calculate Manhattan distance
				distance := utils.Abs(manLat-womanLat) + utils.Abs(manLon-womanLon)

				if distance <= maxDistance {
					eligibleWomen = append(eligibleWomen, struct {
						id       int
						distance int
					}{womanID, distance})
				}
			}
		}

		// If no eligible women, make man wander (improved wandering logic)
		if len(eligibleWomen) == 0 {
			newLat, newLon := utils.Wander(model, manLat, manLon)
			pop.IndData[manID][individual.Lat] = newLat
			pop.IndData[manID][individual.Lon] = newLon
			continue
		}

		// Select woman based on age preference and distance
		womanID := selectWomanByAgeDistanceScore(model, eligibleWomen, manAge, ages, pop)

		// Create marriage
		if womanID != -1 {
			pop.IndData[manID][individual.MarriageState] = womanID
			pop.IndData[womanID][individual.MarriageState] = manID
			pop.IndData[womanID][individual.Lat] = pop.IndData[manID][individual.Lat]
			pop.IndData[womanID][individual.Lon] = pop.IndData[manID][individual.Lon]

			matchedWomen[womanID] = true
			pop.Tracking["marriages"]++
		}
	}
}

// Spatial indexing helpers
type spatialKey struct {
	cellX, cellY int
}

// createSpatialIndex builds a simple grid-based spatial index
func createSpatialIndex(pop *types.Pop, individuals []int) map[spatialKey][]int {
	index := make(map[spatialKey][]int)
	cellSize := 50 // Adjust based on typical maxDistance

	for _, id := range individuals {
		lat := pop.IndData[id][individual.Lat]
		lon := pop.IndData[id][individual.Lon]
		key := spatialKey{lat / cellSize, lon / cellSize}
		index[key] = append(index[key], id)
	}

	return index
}

// getNearbyWomen returns women IDs within the maximum distance
func getNearbyWomen(spatialIndex map[spatialKey][]int, lat, lon, maxDistance int) []int {
	cellSize := 50 // Must match the value used in createSpatialIndex
	centerCell := spatialKey{lat / cellSize, lon / cellSize}

	// Determine cell range to search
	cellRange := (maxDistance / cellSize) + 1
	var nearbyWomen []int

	// Search in neighboring cells
	for dx := -cellRange; dx <= cellRange; dx++ {
		for dy := -cellRange; dy <= cellRange; dy++ {
			key := spatialKey{centerCell.cellX + dx, centerCell.cellY + dy}
			if women, exists := spatialIndex[key]; exists {
				nearbyWomen = append(nearbyWomen, women...)
			}
		}
	}

	return nearbyWomen
}

// getNearbyWomenFromSpatialIndex returns women within range using spatial index
func getNearbyWomenFromSpatialIndex(spatialIndex map[spatialKey][]int, lat, lon, cellSize, maxDistance int) []int {
	centerCell := spatialKey{lat / cellSize, lon / cellSize}
	cellRange := (maxDistance / cellSize) + 1

	var results []int
	for dx := -cellRange; dx <= cellRange; dx++ {
		for dy := -cellRange; dy <= cellRange; dy++ {
			key := spatialKey{centerCell.cellX + dx, centerCell.cellY + dy}
			if women, exists := spatialIndex[key]; exists {
				results = append(results, women...)
			}
		}
	}

	return results
}

// selectWomanByAgeDistanceScore selects a woman based on a combined score of age preference and distance
func selectWomanByAgeDistanceScore(model *types.Model, eligibleWomen []struct{ id, distance int }, manAge int, ages map[int]int, pop *types.Pop) int {
	if len(eligibleWomen) == 0 {
		return -1
	}

	// Calculate scores
	type scoredWoman struct {
		id    int
		score float64
	}

	var scoredWomen []scoredWoman
	totalScore := 0.0

	for _, woman := range eligibleWomen {
		womanID := woman.id
		distance := woman.distance

		femaleAge := ages[womanID]
		ageDiff := manAge - femaleAge

		// Age preference score (higher = better)
		var ageScore float64
		if ageDiff >= 0 && ageDiff <= 5 {
			// Ideal age range: 0-5 years younger
			ageScore = 5.0 - float64(ageDiff)
		} else if ageDiff < 0 {
			// Woman is older than man
			ageScore = math.Max(0.1, 1.0+float64(ageDiff)/2.0)
		} else {
			// Woman is more than 5 years younger
			ageScore = math.Max(0.1, 2.0-float64(ageDiff-5)/10.0)
		}

		// Distance score (higher = better)
		maxDistance := int(model.Parameters["max_mating_distance"])
		distanceScore := float64(maxDistance-distance) / float64(maxDistance)

		// Combined score (with weights)
		// Adjust weights to prioritize age or distance
		combinedScore := (ageScore * 0.7) + (distanceScore * 0.3)

		scoredWomen = append(scoredWomen, scoredWoman{womanID, combinedScore})
		totalScore += combinedScore
	}

	// If total score is effectively zero, return random woman
	if totalScore < 0.001 {
		randomIndex := rand.Intn(len(eligibleWomen))
		return eligibleWomen[randomIndex].id
	}

	// Normalize scores and select using weighted probability
	r := rand.Float64()
	cumulativeProb := 0.0

	for _, sw := range scoredWomen {
		probability := sw.score / totalScore
		cumulativeProb += probability
		if r <= cumulativeProb {
			return sw.id
		}
	}

	// Fallback (rarely reached due to floating point rounding)
	return eligibleWomen[len(eligibleWomen)-1].id
}

// MapBoundaries holds the extents of the map
type MapBoundaries struct {
	MinLat, MaxLat, MinLon, MaxLon int
	TerrainMap                     map[spatialKey]int // Using spatialKey for coordinates
}
