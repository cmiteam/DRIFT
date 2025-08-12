package utils

import (
	"drift/modules/individual"
	"drift/modules/maploader"
	"drift/types"
	"math/rand"
)

// Wander moves an individual to a random location within a square area
// centered on the current position.
func Wander(model *types.Model, lat, lon int) (int, int) {
	// Get maximum wander distance
	maxDistance := int(model.Parameters["wander"])
	if maxDistance < 1 {
		return lat, lon // No wandering if distance is too small
	}

	// Choose random offsets between -maxDistance and +maxDistance
	offsetLat := rand.Intn(2*maxDistance+1) - maxDistance
	offsetLon := rand.Intn(2*maxDistance+1) - maxDistance

	// Calculate new position
	newLat := lat + offsetLat
	newLon := lon + offsetLon

	if maploader.IsLand(model, newLat, newLon) {
		return newLat, newLon
	}
	return lat, lon
}

// abs returns the absolute value of x
func Abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// min returns the smaller of a and b
func Min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func CountGenealo(indData *map[int][]int) int {
	var genealo int
	for _, ind := range *indData {
		if ind[individual.MaxGenealoGens] > -1 {
			genealo++
		}
	}
	return genealo
}
