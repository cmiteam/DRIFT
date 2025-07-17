package animations

import (
	"drift/modules/individual"
	"drift/modules/maploader"
	"drift/types"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/gif"
	"os"
	"path/filepath"
)

// Initialize creates and sets up an animations container
func Initialize(model *types.Model, mapRoot string) (*types.AnimationsContainer, error) {
	// Create the base map
	baseMap, tileSize := CreateBaseMap(model, TerrainColors, mapRoot)
	if baseMap == nil {
		return nil, fmt.Errorf("failed to create base map")
	}

	// Create the animations container
	container := &types.AnimationsContainer{
		BaseMap:    baseMap,
		TileSize:   tileSize,
		Collection: make(map[string]*types.AnimationWriter),
	}

	// Create the standard animations
	CreateAnimation(container, "genetic")
	CreateAnimation(container, "genealogical")

	return container, nil
}

// CreateAnimation adds a new animation to the container
func CreateAnimation(container *types.AnimationsContainer, name string) {
	bounds := container.BaseMap.Bounds()
	width := bounds.Dx()
	height := bounds.Dy()

	container.Collection[name] = &types.AnimationWriter{
		Frames: []*image.Paletted{},
		Delays: []int{},
		Width:  width,
		Height: height,
	}
}

// CreateBaseMap transforms terrain data into a 2D color array

func CreateBaseMap(model *types.Model, terrainColors map[int]color.RGBA, mapRoot string) (*image.RGBA, int) {
	maploader.LoadMap(model, mapRoot)
	mapHeight := len(model.Map)
	mapWidth := len(model.Map[0])
	tileSize := 3600 / mapHeight // Scale to fit 3600 pixel height
	imgWidth := mapWidth * tileSize
	imgHeight := 3610 // 3600 + 10 for progress bar

	// Create the base map image with scaled dimensions
	bounds := image.Rect(0, 0, imgWidth, imgHeight)
	baseMap := image.NewRGBA(bounds)

	// Pre-color the bottom 10 rows with red (for the indicator area)
	for y := bounds.Dy() - 10; y < bounds.Dy(); y++ {
		for x := 0; x < bounds.Dx(); x++ {
			baseMap.SetRGBA(x, y, color.RGBA{255, 0, 0, 255}) // Red background
		}
	}

	// Fill the base map with terrain colors
	for lat := range model.Map {
		for lon := range model.Map[lat] {
			terrainType := model.Map[lat][lon]
			pixelColor, exists := terrainColors[terrainType]
			if !exists {
				pixelColor = color.RGBA{0, 0, 0, 255} // Default color
			}

			// Calculate tile position
			startX := lon * tileSize
			startY := lat * tileSize // Flip Y, account for progress bar and tile height

			// Create rectangle for the terrain tile
			tileRect := image.Rect(startX, startY, startX+tileSize+1, startY+tileSize+1)

			// Check bounds and draw terrain tile
			if startX >= 0 && startY >= 0 && startX+tileSize <= bounds.Dx() && startY+tileSize <= bounds.Dy()-10 {
				uniformColor := &image.Uniform{pixelColor}
				draw.Draw(baseMap, tileRect, uniformColor, image.Point{}, draw.Over)
			}
		}
	}
	fmt.Printf("BaseMap bounds: %v\n", baseMap.Bounds())
	fmt.Printf("Calculated tileSize: %d\n", tileSize)
	fmt.Printf("Map has %d rows, %d columns\n", len(model.Map), len(model.Map[0]))
	fmt.Printf("First few terrain values: %v\n", model.Map[0][:min(5, len(model.Map[0]))])
	return baseMap, tileSize
}

// AddFrame adds a new frame to an animation
func AddFrame(model *types.Model, container *types.AnimationsContainer, animName string, updates map[[2]int]color.RGBA) error {

	// Find the animation
	anim, exists := container.Collection[animName]
	if !exists {
		return fmt.Errorf("animation %s does not exist", animName)
	}

	progressPercent := float64(model.FreeParameters["year"]) / float64(model.Parameters["end_year"])
	model.FreeParameters["progressPercent"] = int(progressPercent * 100)

	// Create a temporary RGBA image with base map
	bounds := container.BaseMap.Bounds()
	rgba := image.NewRGBA(bounds)
	draw.Draw(rgba, bounds, container.BaseMap, image.Point{}, draw.Src)

	// Apply updates with tiles
	for coords, tileColor := range updates {
		mapLat, mapLon := coords[0], coords[1]

		// Convert map coordinates to pixel coordinates
		startX := mapLon * container.TileSize
		startY := mapLat * container.TileSize // No Y-flipping needed

		// Create rectangle for the tile
		tileRect := image.Rect(startX, startY, startX+container.TileSize, startY+container.TileSize)

		// Check if tile is within bounds (leave space for progress bar at bottom)
		if startX >= 0 && startY >= 0 && startX+container.TileSize <= bounds.Dx() && startY+container.TileSize <= bounds.Dy()-10 {
			// Create uniform color image and draw tile
			uniformColor := &image.Uniform{tileColor}
			draw.Draw(rgba, tileRect, uniformColor, image.Point{}, draw.Over)
		}
	}

	barWidth := int(float64(bounds.Dx()) * progressPercent)

	// Draw a 5-pixel thick green progress bar over the red background at the bottom
	for y := bounds.Dy() - 10; y < bounds.Dy(); y++ {
		for x := 0; x < barWidth; x++ {
			rgba.SetRGBA(x, y, color.RGBA{0, 255, 0, 255}) // Bright green
		}
	}

	// Create a palette with your terrain colors plus some extras
	palette := color.Palette{}

	// Add terrain colors first
	for _, c := range TerrainColors {
		palette = append(palette, c)
	}

	// Create a paletted image
	palettedImg := image.NewPaletted(bounds, palette)

	// Use Go's draw package to convert the image
	draw.Draw(palettedImg, bounds, rgba, image.Point{}, draw.Src)

	// Add the frame
	anim.Frames = append(anim.Frames, palettedImg)
	anim.Delays = append(anim.Delays, int(model.Parameters["delay"]))

	return nil
}

// SaveGIF saves the animation as a GIF file
func SaveGIF(container *types.AnimationsContainer, animName string, results string) error {
	// Find the animation
	anim, exists := container.Collection[animName]
	if !exists {
		return fmt.Errorf("animation %s does not exist", animName)
	}

	// Create the output file
	outputPath := filepath.Join(results, animName+".gif")
	f, err := os.Create(outputPath)
	if err != nil {
		return err
	}
	defer f.Close()

	// Create and encode the GIF
	g := &gif.GIF{
		Image:     anim.Frames,
		Delay:     anim.Delays,
		LoopCount: 0, // 0 means loop forever
	}
	return gif.EncodeAll(f, g)
}

// SaveAllGIFs saves all animations in the container
func SaveAllGIFs(container *types.AnimationsContainer, results string) error {
	for name := range container.Collection {
		if err := SaveGIF(container, name, results); err != nil {
			return fmt.Errorf("error saving animation %s: %w", name, err)
		}
	}
	return nil
}

// TerrainColors maps terrain types to colors
var TerrainColors = map[int]color.RGBA{
	0: color.RGBA{0, 0, 0, 255},       // Edge - Black
	1: color.RGBA{76, 153, 0, 255},    // Land - Green
	2: color.RGBA{0, 191, 255, 255},   // CoastalWater - Light blue
	3: color.RGBA{0, 0, 204, 255},     // OpenWater - Deep blue
	4: color.RGBA{255, 255, 255, 255}, // HighMountain - Brown
	5: color.RGBA{255, 204, 102, 255}, // Desert - Sand
	6: color.RGBA{255, 255, 255, 255}, // Ice - White
	7: color.RGBA{0, 0, 0, 255},       // InvalidTerrainHigh - Black
}

func AddAnimationFrames(model *types.Model, pop *types.Pop, container *types.AnimationsContainer) error {
	if model.Parameters["track_map"] != 1 {
		return nil
	}

	// Create genetic frame data
	geneticPoints := createGeneticFrameData(model, pop)
	err := AddFrame(model, container, "genetic", geneticPoints)
	if err != nil {
		return fmt.Errorf("error adding genetic frame: %w", err)
	}

	// Create genealogical frame data (returns map[[2]float64]color.RGBA)
	genealogicalPoints := createGenealogicalFrameData(model, pop)
	err = AddFrame(model, container, "genealogical", genealogicalPoints)
	if err != nil {
		return fmt.Errorf("error adding genealogical frame: %w", err)
	}

	return nil
}

func createGeneticFrameData(model *types.Model, pop *types.Pop) map[[2]int]color.RGBA {
	changePoints := make(map[[2]int]color.RGBA)

	// Create base points for all individuals (white dots)
	for _, indData := range pop.IndData {
		lat := indData[individual.Lat]
		lon := indData[individual.Lon]
		changePoints[[2]int{lat, lon}] = color.RGBA{255, 255, 255, 255}
	}

	// Add genetic markers (red dots) - these will overwrite white dots where present
	for _, indData := range pop.IndData {
		hasMarker := indData[individual.YGens] > 0 ||
			indData[individual.MtGens] > 0 ||
			indData[individual.AlleleCount] > 0 ||
			indData[individual.NumCentromeres] > 0 ||
			indData[individual.NumBlocks] > 0

		if hasMarker {
			lat := indData[individual.Lat] // Pixel coordinates
			lon := indData[individual.Lon] // Pixel coordinates
			changePoints[[2]int{lat, lon}] = color.RGBA{255, 0, 0, 255}
		}
	}

	return changePoints
}

func createGenealogicalFrameData(model *types.Model, pop *types.Pop) map[[2]int]color.RGBA {
	changePoints := make(map[[2]int]color.RGBA)

	// Create base points for all individuals (white dots)
	for _, indData := range pop.IndData {
		lat := indData[individual.Lat]
		lon := indData[individual.Lon]
		changePoints[[2]int{lat, lon}] = color.RGBA{255, 255, 255, 255}
	}

	// Add genealogical markers (red dots)
	for _, indData := range pop.IndData {
		hasMarker := indData[individual.YGens] > 0 ||
			indData[individual.MtGens] > 0 ||
			indData[individual.MinGenealoGens] > 0 ||
			indData[individual.MaxGenealoGens] > 0

		if hasMarker {
			lat := indData[individual.Lat]
			lon := indData[individual.Lon]
			changePoints[[2]int{lat, lon}] = color.RGBA{255, 0, 0, 255}
		}
	}

	return changePoints
}
