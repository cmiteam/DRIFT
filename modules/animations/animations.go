package animations

import (
	"drift/modules/maploader"
	"drift/types"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/gif"
	"math"
	"os"
	"path/filepath"
)

// Initialize creates and sets up an animations container
func Initialize(model *types.Model, mapRoot string) (*types.AnimationsContainer, error) {
	// Create the base map
	baseMap, minLat, minLon, maxLat, maxLon := CreateBaseMap(model, TerrainColors, mapRoot)
	if baseMap == nil {
		return nil, fmt.Errorf("failed to create base map")
	}

	// Create the animations container
	container := &types.AnimationsContainer{
		BaseMap:    baseMap,
		Collection: make(map[string]*types.AnimationWriter),
		MinLat:     minLat,
		MinLon:     minLon,
		MaxLat:     maxLat,
		MaxLon:     maxLon,
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

func CreateBaseMap(model *types.Model, terrainColors map[int]color.RGBA, mapRoot string) (*image.RGBA, int, int, int, int) {
	err, minLat, minLon, maxLat, maxLon := maploader.LoadMap(model, mapRoot)
	if err != nil {
		return nil, 0, 0, 0, 0
	}

	// Get pixel scale from model parameters
	pixelSize := 1
	if size, exists := model.Parameters["pixel_size"]; exists {
		pixelSize = int(size)
		if pixelSize < 1 {
			pixelSize = 1 // Ensure minimum scale is 1
		}
	}

	imgWidth := (maxLon - minLon + 1) * pixelSize
	imgHeight := (maxLat-minLat+1)*pixelSize + 10 // +10 for progress bar

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

			// Transform coordinates
			x, y := TransformCoordinates(lat, lon, minLat, minLon, maxLat, maxLon, pixelSize)

			// Draw the pixel (with size scaling)
			for dy := 0; dy < pixelSize; dy++ {
				for dx := 0; dx < pixelSize; dx++ {
					if x+dx < bounds.Dx() && y+dy < bounds.Dy()-10 {
						baseMap.SetRGBA(x+dx, y+dy, pixelColor)
					}
				}
			}
		}
	}

	return baseMap, minLat, minLon, maxLat, maxLon
}

// AddFrame adds a new frame to an animation
func AddFrame(model *types.Model, container *types.AnimationsContainer, animName string, updates map[[2]int]color.RGBA) error {

	// Find the animation
	anim, exists := container.Collection[animName]
	if !exists {
		return fmt.Errorf("animation %s does not exist", animName)
	}
	//	delay := model.Parameters["delay"]
	pixelSize := int(model.Parameters["pixel_size"])
	progressPercent := float64(model.FreeParameters["year"]) / float64(model.Parameters["end_year"])
	model.FreeParameters["progressPercent"] = int(progressPercent * 100)
	minLat := container.MinLat
	minLon := container.MinLon
	maxLat := container.MaxLat
	maxLon := container.MaxLon
	pointSize := int(model.Parameters["point_size"])

	// Create a temporary RGBA image with base map
	bounds := container.BaseMap.Bounds()
	rgba := image.NewRGBA(bounds)
	draw.Draw(rgba, bounds, container.BaseMap, image.Point{}, draw.Src)

	// Apply updates with scaling
	for coords, pixelColor := range updates {
		lat, lon := coords[0], coords[1]
		x, y := TransformCoordinates(lat, lon, minLat, minLon, maxLat, maxLon, pixelSize)
		//		fmt.Printf("  Drawing from lat/lon (%.1f, %.1f) to pixel x/y (%d, %d)\n", lat, lon, x, y)

		// Check bounds and draw each update
		if x >= 0 && x < bounds.Dx() && y >= 0 && y < bounds.Dy() {
			for dy := 0; dy < pointSize; dy++ {
				for dx := 0; dx < pointSize; dx++ {
					newX := x + dx
					newY := y + dy
					if newX >= 0 && newX < bounds.Dx() && newY >= 0 && newY < bounds.Dy() {
						rgba.SetRGBA(newX, newY, pixelColor)
					}
				}
			}
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

	// Add common colors
	for _, c := range []color.RGBA{
		{0, 0, 0, 255},       // Black
		{255, 255, 255, 255}, // White
		{255, 0, 0, 255},     // Red
		{0, 255, 0, 255},     // Green
		{0, 0, 255, 255},     // Blue
		{255, 255, 0, 255},   // Yellow
		{255, 0, 255, 255},   // Magenta
		{0, 255, 255, 255},   // Cyan
		{128, 128, 128, 255}, // Gray
	} {
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

// Transform from geographic coordinates to image pixel coordinates
func TransformCoordinates(lat, lon, minLat, minLon, maxLat, maxLon, pixelSize int) (x, y int) {

	latRange := maxLat - minLat
	if latRange <= 0 {
		latRange = 1
	}
	lonRange := maxLon - minLon
	if lonRange <= 0 {
		lonRange = 1
	}
	// Calculate normalized position within bounds (0.0 to 1.0)
	normalizedX := float64(lon-minLon) / float64(lonRange)
	normalizedY := float64(lat-minLat) / float64(latRange)

	// Calculate image dimensions
	imgWidth := lonRange * pixelSize
	imgHeight := latRange * pixelSize

	// Convert to pixel coordinates and scale by pixelSize
	x = int(math.Round(normalizedX * float64(imgWidth)))
	y = int(math.Round((1.0 - normalizedY) * float64(imgHeight)))
	return x, y
}

// TerrainColors maps terrain types to colors
var TerrainColors = map[int]color.RGBA{
	1: color.RGBA{76, 153, 0, 255},    // Land - Green
	2: color.RGBA{0, 191, 255, 255},   // CoastalWater - Light blue
	3: color.RGBA{0, 0, 204, 255},     // OpenWater - Deep blue
	4: color.RGBA{255, 255, 255, 255}, // HighMountain - Brown
	5: color.RGBA{255, 204, 102, 255}, // Desert - Sand
	6: color.RGBA{255, 255, 255, 255}, // Ice - White
	0: color.RGBA{0, 0, 0, 255},       // InvalidTerrainLow - Black
	7: color.RGBA{0, 0, 0, 255},       // InvalidTerrainHigh - Black
}
