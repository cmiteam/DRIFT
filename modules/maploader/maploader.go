package maploader

import (
	"drift/modules/csvutils"
	"drift/types"
	"fmt"
	"math"
	"sort"
	"strconv"
)

// Terrain types
type Terrain int

const (
	InvalidTerrainLow Terrain = iota // 0
	Land                             // 1
	CoastalWater                     // 2
	OpenWater                        // 3
	HighMountain                     // 4
	Desert                           // 5
	Ice                              // 6
	InvalidTerrainHigh
)

// Load the map from a CSV file and populate the model's Map map.
func LoadMap(model *types.Model, mapRoot string) error {
	minLat, minLon := 1000000.0, 1000000.0
	maxLat, maxLon := -1000000.0, -1000000.0

	// Initialize the map if it's nil
	if model.Map == nil {
		model.Map = make(map[int]map[int]int)
	}

	// Derive the filename from the model and load the CSV file
	filename := fmt.Sprintf("%s.csv", model.MapName)
	csvLoader := csvutils.CSVLoader{
		FileName:   filename,
		Dir:        mapRoot,
		MinRecords: 2,
	}
	records, err := csvLoader.LoadCSV()
	if err != nil {
		return err
	}

	fmt.Printf("Loading map: %s (%d records)\n", filename, len(records)-1)

	// Store all coordinates for processing
	allCoordinates := make([]struct {
		lat, lon float64
		terrain  int
	}, 0, len(records)-1)

	// Sample lat/lon values to determine resolution
	latSamples := make([]float64, 0, 100)
	lonSamples := make([]float64, 0, 100)

	for _, record := range records[1:] {
		lat, err := strconv.ParseFloat(record[0], 64)
		if err != nil {
			return err
		}

		lon, err := strconv.ParseFloat(record[1], 64)
		if err != nil {
			return err
		}

		terrain, err := strconv.ParseInt(record[2], 10, 64)
		if err != nil {
			return err
		}

		// Store for later processing
		allCoordinates = append(allCoordinates, struct {
			lat, lon float64
			terrain  int
		}{lat, lon, int(terrain)})

		// Track min values
		if lat < minLat {
			minLat = lat
		}
		if lon < minLon {
			minLon = lon
		}

		// Add to samples (limit to 100 samples for efficiency)
		if len(latSamples) < 100 {
			latSamples = append(latSamples, lat)
		}
		if len(lonSamples) < 100 {
			lonSamples = append(lonSamples, lon)
		}
	}

	// Determine scale factor based on coordinate patterns
	scaleFactor := calculateScaleFactor(latSamples, lonSamples)
	fmt.Printf("Detected coordinate resolution, using scale factor: %f\n", scaleFactor)
	fmt.Printf("Coordinate offsets: lat=%f, lon=%f\n", minLat, minLon)

	// Process and store the data with offsets
	for _, coord := range allCoordinates {
		// Apply offsets and convert to integers
		adjustedLat := int(math.Round((coord.lat - minLat) * scaleFactor))
		adjustedLon := int(math.Round((coord.lon - minLon) * scaleFactor))

		// Track the maximum values after adjustment
		if float64(adjustedLat) > maxLat {
			maxLat = float64(adjustedLat)
		}
		if float64(adjustedLon) > maxLon {
			maxLon = float64(adjustedLon)
		}

		// Store in model.Map with adjusted coordinates
		if model.Map[adjustedLat] == nil {
			model.Map[adjustedLat] = make(map[int]int)
		}
		model.Map[adjustedLat][adjustedLon] = coord.terrain
	}

	// Store the coordinate info
	model.FreeParameters["maxLat"] = int(maxLat)
	model.FreeParameters["maxLon"] = int(maxLon)

	// Calculate tile size for display
	latRange := int(maxLat) + 1
	lonRange := int(maxLon) + 1
	width := model.Parameters["map_width"]
	height := model.Parameters["map_height"]
	model.FreeParameters["tileWidth"] = int(width / float64(lonRange))
	model.FreeParameters["tileHeight"] = int(height / float64(latRange))

	fmt.Printf("Map loaded with dimensions: %d×%d\n", latRange, lonRange)

	return nil
}

// calculateScaleFactor analyzes coordinate samples to find the appropriate scale factor
func calculateScaleFactor(latSamples, lonSamples []float64) float64 {
	// Find the smallest decimal increment in a slice of floats
	findSmallestIncrement := func(samples []float64) float64 {
		if len(samples) < 2 {
			return 1.0 // Default if not enough samples
		}

		// Sort the samples
		sort.Float64s(samples)

		// Find the smallest non-zero difference
		minDiff := 1000.0
		for i := 1; i < len(samples); i++ {
			diff := math.Abs(samples[i] - samples[i-1])
			if diff > 0 && diff < minDiff {
				minDiff = diff
			}
		}

		if minDiff >= 1000.0 {
			return 1.0 // No meaningful difference found
		}

		return minDiff
	}

	// Get the smallest increments
	latIncrement := findSmallestIncrement(latSamples)
	lonIncrement := findSmallestIncrement(lonSamples)

	// Use the smaller of the two increments to ensure we don't lose precision
	increment := math.Min(latIncrement, lonIncrement)

	// Check for common patterns
	if almostEqual(increment, 0.1) {
		return 10.0
	} else if almostEqual(increment, 0.01) {
		return 100.0
	} else if almostEqual(increment, 0.001) {
		return 1000.0
	} else if almostEqual(increment, 0.5) {
		return 2.0
	} else if almostEqual(increment, 0.25) {
		return 4.0
	} else if almostEqual(increment, 0.2) {
		return 5.0
	} else if increment > 0 {
		// For other patterns, use 1/increment as scale factor
		return 1.0 / increment
	}

	// Default fallback
	return 1.0
}

// almostEqual checks if two floats are nearly equal
func almostEqual(a, b float64) bool {
	epsilon := 1e-9
	return math.Abs(a-b) < epsilon
}
