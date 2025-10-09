package utils

import (
	"drift/pkg/core"
	"fmt"
	"strconv"
	"strings"
    "math/rand"
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

// LoadMap loads the map from a CSV file in grid format where each row = latitude, each column = longitude
// Expected format:
// 1,1,1,1,1,1,1
// 1,1,1,1,5,1,1
// 1,1,1,1,1,1,1
func LoadMap(model *core.Model, mapRoot string) error {
	// Derive the filename from the model and load the CSV file
	filename := fmt.Sprintf("%s.csv", model.MapName)
    records, err := LoadCSV(filename)
	if err != nil {
		return err
	}
	fmt.Printf("CSV loaded - %d rows, first row: %v\n", len(records), records[0][:min(5, len(records[0]))])
	// Get dimensions from the CSV data
	mapHeight := len(records)
	if mapHeight == 0 {
		return fmt.Errorf("empty map file")
	}

	mapWidth := len(records[0])
	if mapWidth == 0 {
		return fmt.Errorf("empty map row")
	}

	// Validate that all rows have the same width
	for i, record := range records {
		if len(record) != mapWidth {
			return fmt.Errorf("inconsistent row width at row %d: expected %d, got %d", i, mapWidth, len(record))
		}
	}

	// Initialize the map as a 2D slice
	model.Map = make([][]int, mapHeight)
	for lat := 0; lat < mapHeight; lat++ {
		model.Map[lat] = make([]int, mapWidth)
	}

	// Parse the terrain data
	for csvRow, record := range records {
		flippedLat := mapHeight - 1 - csvRow // Calculate flipped latitude
		for lon, cellStr := range record {
			cellStr = strings.TrimSpace(cellStr)
			if csvRow == 0 && lon == 0 { // Changed from 'lat' to 'csvRow'
				cellStr = strings.TrimPrefix(cellStr, "\ufeff") // Remove UTF-8 BOM
			}
			terrain, err := strconv.Atoi(cellStr)
			if err != nil {
				return fmt.Errorf("invalid terrain value '%s' at CSV position (%d,%d): %v", cellStr, csvRow, lon, err)
			}

			// Validate terrain type
			if terrain < int(InvalidTerrainLow) || terrain >= int(InvalidTerrainHigh) {
				return fmt.Errorf("terrain value %d out of valid range at CSV position (%d,%d)", terrain, csvRow, lon)
			}

			model.Map[flippedLat][lon] = terrain // Store at flipped position
		}
	}

	// Count terrain types for statistics
	terrainCounts := make(map[int]int)
	landCells := 0
	for lat := 0; lat < mapHeight; lat++ {
		for lon := 0; lon < mapWidth; lon++ {
			terrain := model.Map[lat][lon]
			terrainCounts[terrain]++
			if terrain == int(Land) {
				landCells++
			}
		}
	}

	model.Parameters["map_land_cells"] = float64(landCells)
	model.Parameters["map_land_fraction"] = float64(landCells) / float64(mapWidth*mapHeight)

	fmt.Println("=== MAP LOADED (Grid Format) ===")
	fmt.Printf("Map dimensions: %d x %d cells\n", mapWidth, mapHeight)
	fmt.Printf("Total cells: %d\n", mapWidth*mapHeight)
	fmt.Printf("Land cells: %d (%.1f%%)\n", landCells, model.Parameters["map_land_fraction"]*100)
	fmt.Printf("Terrain distribution:\n")
	for terrain, count := range terrainCounts {
		percentage := float64(count) / float64(mapWidth*mapHeight) * 100
		fmt.Printf("  Type %d: %d cells (%.1f%%)\n", terrain, count, percentage)
	}
	return nil
}

// IsValidPosition checks if the given coordinates are within map bounds
func IsValidPosition(model *core.Model, lat, lon int) bool {
	if lat < 0 || lon < 0 {
		return false
	}
	if lat >= len(model.Map) {
		return false
	}
	if lon >= len(model.Map[0]) {
		return false
	}
	return true
}

// GetTerrain returns the terrain type at the given coordinates
func GetTerrain(model *core.Model, lat, lon int) int {
	if !IsValidPosition(model, lat, lon) {
		return int(InvalidTerrainLow)
	}
	return model.Map[lat][lon]
}

// IsLand checks if the given coordinates are on land
func IsLand(model *core.Model, lat, lon int) bool {
	return GetTerrain(model, lat, lon) == int(Land)
}

// FindLandCells returns all land cell coordinates
func FindLandCells(model *core.Model) [][2]int {
	var landCells [][2]int
	for lat := 0; lat < len(model.Map); lat++ {
		for lon := 0; lon < len(model.Map[lat]); lon++ {
			if model.Map[lat][lon] == int(Land) {
				landCells = append(landCells, [2]int{lat, lon})
			}
		}
	}
	return landCells
}

// GetPixelCoordinates converts map coordinates to pixel coordinates for image saving
func GetPixelCoordinates(model *core.Model, mapLat, mapLon int) (int, int) {
	pixelLat := int(float64(mapLat) * model.Parameters["save_lat_scale"])
	pixelLon := int(float64(mapLon) * model.Parameters["save_lon_scale"])
	return pixelLat, pixelLon
}

// GetMapCoordinates converts pixel coordinates back to map coordinates
func GetMapCoordinates(model *core.Model, pixelLat, pixelLon int) (int, int) {
	mapLat := int(float64(pixelLat) / model.Parameters["save_lat_scale"])
	mapLon := int(float64(pixelLon) / model.Parameters["save_lon_scale"])
	return mapLat, mapLon
}

// GetMapDimensions returns the map width and height
func GetMapDimensions(model *core.Model) (int, int) {
	if len(model.Map) == 0 {
		return 0, 0
	}
	return len(model.Map[0]), len(model.Map)
}

// Wander moves an individual to a random location within a square area
// centered on the current position.
func Wander(model *core.Model, lat, lon int) (int, int) {
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

	if IsLand(model, newLat, newLon) {
		return newLat, newLon
	}
	return lat, lon
}