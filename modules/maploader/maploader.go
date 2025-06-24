package maploader

import (
	"drift/modules/csvutils"
	"drift/types"
	"fmt"
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
func LoadMap(model *types.Model, mapRoot string) (error, int, int, int, int) {
	minLat, minLon := 1000000, 1000000
	maxLat, maxLon := 0, 0

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
		return err, minLat, minLon, maxLat, maxLon
	}

	for _, record := range records[1:] {
		lat, err := strconv.Atoi(record[0])
		if err != nil {
			return err, minLat, minLon, maxLat, maxLon
		}

		lon, err := strconv.Atoi(record[1])
		if err != nil {
			return err, minLat, minLon, maxLat, maxLon
		}

		terrain, err := strconv.ParseInt(record[2], 10, 64)
		if err != nil {
			return err, minLat, minLon, maxLat, maxLon
		}

		// Track min/max values
		if lat < minLat {
			minLat = lat
		}
		if lat > maxLat {
			maxLat = lat
		}
		if lon < minLon {
			minLon = lon
		}
		if lon > maxLon {
			maxLon = lon
		}

		// Store coordinates in model.Map
		if model.Map[lat] == nil {
			model.Map[lat] = make(map[int]int)
		}
		model.Map[lat][lon] = int(terrain)
	}

	fmt.Println("=== MAP VISUALIZATION DIAGNOSTICS ===")
	fmt.Printf("Map dimensions: %d x %d\n", maxLon-minLon+1, maxLat-minLat+1)
	fmt.Printf("Map boundaries: minLat=%d, minLon=%d, maxLat=%d, maxLon=%d\n",
		minLat, minLon, maxLat, maxLon)

	return err, minLat, minLon, maxLat, maxLon
}
