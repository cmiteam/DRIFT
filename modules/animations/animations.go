package animations

import (
	"drift/types"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/gif"
	"os"
	"path/filepath"
)

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

// CreateBaseMap transforms terrain data into a 2D color array
func CreateBaseMap(model *types.Model, terrainColors map[int]color.RGBA) [][]color.RGBA {
	// Get map dimensions from model
	maxLat := model.FreeParameters["maxLat"]
	maxLon := model.FreeParameters["maxLon"]

	// Create the base map array
	baseMap := make([][]color.RGBA, maxLat+1)
	for i := range baseMap {
		baseMap[i] = make([]color.RGBA, maxLon+1)
	}

	// Fill the base map with terrain colors
	for lat := 0; lat <= maxLat; lat++ {
		for lon := 0; lon <= maxLon; lon++ {
			// Get terrain type from model.Map
			terrainType := model.Map[lat][lon]

			// Get color from the provided mapping
			pixelColor, exists := terrainColors[terrainType]
			if !exists {
				// Default color for unknown terrain types
				pixelColor = color.RGBA{0, 0, 0, 255} // Black
			}

			// Set the color in the base map
			baseMap[maxLat-lat][lon] = pixelColor
		}
	}

	return baseMap
}

// Initialize creates and sets up an animations container
func Initialize(model *types.Model) *types.AnimationsContainer {
	// Create the base map
	baseMap := CreateBaseMap(model, TerrainColors)

	// Create the animations container
	pixelSize := 10
	container := &types.AnimationsContainer{
		BaseMap:    baseMap,
		Collection: make(map[string]*types.AnimationWriter),
		PixelSize:  pixelSize,
	}

	// Create the standard animations
	CreateAnimation(container, "genetic")
	CreateAnimation(container, "genealogical")

	return container
}

// CreateAnimation adds a new animation to the container
func CreateAnimation(container *types.AnimationsContainer, name string) {
	height := len(container.BaseMap)
	width := 0
	if height > 0 {
		width = len(container.BaseMap[0])
	}

	container.Collection[name] = &types.AnimationWriter{
		Frames: []*image.Paletted{},
		Delays: []int{},
		Width:  width,
		Height: height,
	}
}

// AddFrame adds a new frame to an animation
func AddFrame(container *types.AnimationsContainer, animName string, updates map[[2]int]color.RGBA, delay int) error {
	// Find the animation
	anim, exists := container.Collection[animName]
	if !exists {
		return fmt.Errorf("animation %s does not exist", animName)
	}

	// Create a temporary RGBA image with base map
	bounds := image.Rect(0, 0, anim.Width, anim.Height)
	rgba := image.NewRGBA(bounds)

	// Copy base map colors
	for y := 0; y < len(container.BaseMap); y++ {
		for x := 0; x < len(container.BaseMap[0]); x++ {
			rgba.SetRGBA(x, y, container.BaseMap[y][x])
		}
	}

	// Apply updates
	pixelSize := container.PixelSize
	for coords, pixelColor := range updates {
		lat, lon := coords[0], coords[1]
		if lat >= 0 && lat < len(container.BaseMap) &&
			lon >= 0 && lon < len(container.BaseMap[0]) {
			// Draw a larger pixel (pixelSize x pixelSize square)
			for dy := 0; dy < pixelSize; dy++ {
				for dx := 0; dx < pixelSize; dx++ {
					// Check boundaries to avoid out-of-bounds errors
					newLat := lat + dy
					newLon := lon + dx
					if newLat < len(container.BaseMap) && newLon < len(container.BaseMap[0]) {
						rgba.SetRGBA(newLon, newLat, pixelColor)
					}
				}
			}
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
	anim.Delays = append(anim.Delays, delay)

	return nil
}

// SaveGIF saves the animation as a GIF file
func SaveGIF(container *types.AnimationsContainer, animName string) error {
	// Find the animation
	anim, exists := container.Collection[animName]
	if !exists {
		return fmt.Errorf("animation %s does not exist", animName)
	}

	// Create the output file
	outputPath := filepath.Join("Results", animName+".gif")
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
func SaveAllGIFs(container *types.AnimationsContainer) error {
	for name := range container.Collection {
		if err := SaveGIF(container, name); err != nil {
			return fmt.Errorf("error saving animation %s: %w", name, err)
		}
	}
	return nil
}
