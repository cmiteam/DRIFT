// Takes three columns of data (lat, lon, terrain) and converts it to a grid of terrain values
// Lattitudes will be inverted to account for maps (high lat is at the top) vs images (high y is at the bottom)

package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
	"strings"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Println("Usage: go run converter.go <input_file.csv> [output_file.csv]")
		fmt.Println("Example: go run converter.go old_map.csv new_map.csv")
		os.Exit(1)
	}

	inputFile := os.Args[1]

	// Default output filename
	outputFile := "converted_map.csv"
	if len(os.Args) >= 3 {
		outputFile = os.Args[2]
	}

	err := convertMap(inputFile, outputFile)
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("Successfully converted %s to %s\n", inputFile, outputFile)
}

func convertMap(inputFile, outputFile string) error {
	// Read the input CSV file
	file, err := os.Open(inputFile)
	if err != nil {
		return fmt.Errorf("failed to open input file: %v", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		return fmt.Errorf("failed to read CSV: %v", err)
	}

	if len(records) == 0 {
		return fmt.Errorf("empty CSV file")
	}

	// Parse the data and find bounds
	type Point struct {
		lat, lon, terrain int
	}

	var points []Point
	minLat, maxLat := 999999, -999999
	minLon, maxLon := 999999, -999999

	// Skip header row (start from index 1)
	startIndex := 0
	if len(records) > 0 {
		// Check if first row looks like headers (non-numeric)
		firstRow := records[0]
		if len(firstRow) >= 3 {
			if _, err := strconv.Atoi(strings.TrimSpace(firstRow[0])); err != nil {
				fmt.Printf("Detected header row: %v\n", firstRow)
				startIndex = 1
			}
		}
	}

	for i := startIndex; i < len(records); i++ {
		record := records[i]
		if len(record) < 3 {
			return fmt.Errorf("row %d has fewer than 3 columns", i+1)
		}

		// Parse lat, lon, terrain
		lat, err := strconv.Atoi(strings.TrimSpace(record[0]))
		if err != nil {
			return fmt.Errorf("invalid latitude '%s' at row %d: %v", record[0], i+1, err)
		}

		lon, err := strconv.Atoi(strings.TrimSpace(record[1]))
		if err != nil {
			return fmt.Errorf("invalid longitude '%s' at row %d: %v", record[1], i+1, err)
		}

		terrain, err := strconv.Atoi(strings.TrimSpace(record[2]))
		if err != nil {
			return fmt.Errorf("invalid terrain '%s' at row %d: %v", record[2], i+1, err)
		}

		points = append(points, Point{lat, lon, terrain})

		// Update bounds
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
	}

	fmt.Printf("Input data bounds: lat[%d,%d], lon[%d,%d]\n", minLat, maxLat, minLon, maxLon)
	fmt.Printf("Grid size: %d x %d\n", maxLat-minLat+1, maxLon-minLon+1)

	// Create the grid
	height := maxLat - minLat + 1
	width := maxLon - minLon + 1
	grid := make([][]int, height)
	for i := range grid {
		grid[i] = make([]int, width)
		// Initialize with default terrain (0 for edge/invalid)
		for j := range grid[i] {
			grid[i][j] = 0
		}
	}

	// Fill the grid
	for _, point := range points {
		gridRow := maxLat - point.lat
		gridCol := point.lon - minLon
		grid[gridRow][gridCol] = point.terrain
	}

	// Write the output CSV
	outFile, err := os.Create(outputFile)
	if err != nil {
		return fmt.Errorf("failed to create output file: %v", err)
	}
	defer outFile.Close()

	writer := csv.NewWriter(outFile)
	defer writer.Flush()

	for _, row := range grid {
		strRow := make([]string, len(row))
		for j, val := range row {
			strRow[j] = strconv.Itoa(val)
		}
		if err := writer.Write(strRow); err != nil {
			return fmt.Errorf("failed to write row: %v", err)
		}
	}

	// Print some statistics
	terrainCounts := make(map[int]int)
	for _, row := range grid {
		for _, terrain := range row {
			terrainCounts[terrain]++
		}
	}

	fmt.Println("Terrain distribution:")
	for terrain, count := range terrainCounts {
		percentage := float64(count) / float64(height*width) * 100
		fmt.Printf("  Type %d: %d cells (%.1f%%)\n", terrain, count, percentage)
	}

	return nil
}
