package save

import (
	"drift/modules/animations"
	"drift/types"
	"encoding/csv"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"os"
	"strings"
)

// SaveHeaders creates a CSV file with the headers for the results
func SaveHeaders(modelName string) error {
	filename := fmt.Sprintf("results/%s_results.csv", modelName)
	file, err := os.OpenFile(filename, os.O_TRUNC|os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		return fmt.Errorf("failed to open file: %v", err)
	}
	defer file.Close()
	writer := csv.NewWriter(file)
	defer writer.Flush()
	headers := []string{
		"run", "year", "n", "marrs", "births", "randDs", "cullDs", "GenetDes", "GeneaDes", "YDes", "MtDes",
		"nCents", "nAlleles", "nBlocks", "TotFitness", "nMuts", "PercSeedGenoRet", "AvSeedGenoCov", "AvHet",
	}
	if err := writer.Write(headers); err != nil {
		return fmt.Errorf("failed to write headers: %v", err)
	}
	return nil
}

// Save writes the current simulation state to a CSV file
func Save(model *types.Model, pop *types.Pop, animManager *types.AnimationsContainer) {
	fmt.Printf("   Year: %d  n: %d  b: %d  m: %d c: %d\n",
		model.FreeParameters["year"],
		model.FreeParameters["last_pop_size"],
		pop.Tracking["births"],
		pop.Tracking["marriages"],
		pop.Tracking["cull_deaths"],
	)

	filename := fmt.Sprintf("results/%s_results.csv", model.ModelName)
	numInds := len(pop.IndData)
	var YDescends, mtDescends, genealoDescends, geneticDescends, numAlleles, numBlocks, numCentromeres, numMutations, popFitness int
	var percSeedGenomeRetained, avSeedGenomeCoverage float64
	var totHet, totHomMin, totHomMaj, numbitsRetained int

	if model.Parameters["track_DNA"] == 1 {
		YDescends, mtDescends, genealoDescends, geneticDescends, numAlleles, numBlocks, numCentromeres = calculateMiscStats(pop.IndData)
		numbitsRetained, totHet, totHomMin, totHomMaj = seedCounts(model, pop)
		percSeedGenomeRetained = float64(numbitsRetained) / float64(model.FreeParameters["NumBits"]) * 100
		avSeedGenomeCoverage = 0
	}

	if model.Parameters["track_mutations"] == 1 {
		numMutations, popFitness = calculateFitnessStats(model, pop)
	}

	if model.Parameters["track_map"] == 1 {

		allPoints := make(map[[2]float64]color.RGBA)
		freqMap := make(map[float64]map[float64]int)
		changePoints := make(map[[2]float64]color.RGBA)
		progressPercent := float64(model.FreeParameters["year"]) / float64(model.Parameters["end_year"])
		model.FreeParameters["progressPercent"] = int(progressPercent * 100)

		// Count individuals with genetic markers at each location
		for _, individual := range pop.IndData {
			lat := float64(individual["lat"]) / 10.0
			lon := float64(individual["lon"]) / 10.0
			allPoints[[2]float64{lat, lon}] = color.RGBA{255, 255, 255, 255}

			// Check if this individual has any of the genetic markers we're tracking
			hasMarker := individual["Y_gens"] > 0 ||
				individual["mt_gens"] > 0 ||
				individual["allele_count"] > 0 ||
				individual["num_centromeres"] > 0 ||
				individual["num_blocks"] > 0

			if hasMarker {
				// Initialize the inner map if needed
				if _, exists := freqMap[lat]; !exists {
					freqMap[lat] = make(map[float64]int)
				}
				// Increment the counter
				freqMap[lat][lon]++
			}
		}
		for coords, color := range allPoints {
			changePoints[coords] = color
		}

		// Add your specific data points based on population data
		//totalPopSize := model.FreeParameters["last_pop_size"]
		for lat, innerMap := range freqMap {
			for lon, count := range innerMap {
				if count > 0 {
					changePoints[[2]float64{lat, lon}] = color.RGBA{255, 0, 0, 255}
					// Calculate color based on frequency
					//ratio := float64(count) / float64(totalPopSize)
					//					minRedness := 50
					//intensity := 255 - int(ratio*255)
					//colorValue := (0xFF << 16) | (intensity << 8) | intensity
					//r := uint8((colorValue >> 16) & 0xFF)
					//g := uint8((colorValue >> 8) & 0xFF)
					//b := uint8(colorValue & 0xFF)
					// Add this point to our changed points
					//changePoints[[2]float64{lat, lon}] = color.RGBA{r, g, b, 255}
				}
			}
		}

		err := animations.AddFrame(model, animManager, "genetic", changePoints)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error adding frame to Genetic animation: %v\n", err)
			os.Exit(1)
		}

		// Count individuals with genealogical markers at each location
		freqMap = make(map[float64]map[float64]int)
		changePoints = make(map[[2]float64]color.RGBA)
		for coords, color := range allPoints {
			changePoints[coords] = color
		}

		for _, individual := range pop.IndData {
			// Check if this individual has any of the genealogical markers we're tracking
			hasMarker := individual["Y_gens"] > 0 ||
				individual["mt_gens"] > 0 ||
				individual["min_genealo_gens"] > 0 ||
				individual["max_genealo_gens"] > 0

			if hasMarker {
				lat := float64(individual["lat"] / 10.0)
				lon := float64(individual["lon"] / 10.0)

				// Initialize the inner map if needed
				if _, exists := freqMap[lat]; !exists {
					freqMap[lat] = make(map[float64]int)
				}

				// Increment the counter
				freqMap[lat][lon]++
			}
		}

		// Prepare data points for genealogical animation
		//totalPopSize = model.FreeParameters["last_pop_size"]
		for lat, innerMap := range freqMap {
			for lon, count := range innerMap {
				if count > 0 {
					changePoints[[2]float64{lat, lon}] = color.RGBA{255, 0, 0, 255}
					// Calculate color based on frequency
					//ratio := float64(count) / float64(totalPopSize)
					//intensity := 255 - int(ratio*255)
					//colorValue := (0xFF << 16) | (intensity << 8) | intensity
					//r := uint8((colorValue >> 16) & 0xFF)
					//g := uint8((colorValue >> 8) & 0xFF)
					//b := uint8(colorValue & 0xFF)

					// Add this point to our changed points
					//changePoints[[2]float64{lat, lon}] = color.RGBA{r, g, b, 255}
				}
			}
		}
		// For testing purposes, set a bright red pixel at (year, year)
		//		xy := float64(model.FreeParameters["year"]) / 100.0
		//		changePoints[[2]float64{xy, xy}] = color.RGBA{255, 0, 0, 255}

		err = animations.AddFrame(model, animManager, "genealogical", changePoints)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error adding frame to Genetic animation: %v\n", err)
		}
	}

	file, _ := os.OpenFile(filename, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	defer file.Close()
	writer := csv.NewWriter(file)
	defer writer.Flush()

	data := []string{
		fmt.Sprintf("%d", model.FreeParameters["run"]),
		fmt.Sprintf("%d", model.FreeParameters["year"]),
		fmt.Sprintf("%d", numInds),
		fmt.Sprintf("%d", pop.Tracking["marriages"]),
		fmt.Sprintf("%d", pop.Tracking["births"]),
		fmt.Sprintf("%d", pop.Tracking["random_deaths"]),
		fmt.Sprintf("%d", pop.Tracking["cull_deaths"]),
		fmt.Sprintf("%d", geneticDescends),
		fmt.Sprintf("%d", genealoDescends),
		fmt.Sprintf("%d", YDescends),
		fmt.Sprintf("%d", mtDescends),
		fmt.Sprintf("%d", numCentromeres),
		fmt.Sprintf("%d", numAlleles),
		fmt.Sprintf("%d", numBlocks),
		fmt.Sprintf("%d", popFitness),
		fmt.Sprintf("%d", numMutations),
		fmt.Sprintf("%.1f", percSeedGenomeRetained),
		fmt.Sprintf("%.1f", avSeedGenomeCoverage),
		fmt.Sprintf("%.2f", totHet),
		fmt.Sprintf("%.2f", totHomMin),
		fmt.Sprintf("%.2f", totHomMaj),
	}
	writer.Write(data)

	pop.Tracking["births"] = 0
	pop.Tracking["deaths"] = 0
	pop.Tracking["marriages"] = 0
	pop.Tracking["random_deaths"] = 0
	pop.Tracking["cull_deaths"] = 0
}

// personDataString formats detailed individual data with state information
func personDataString(model *types.Model, pop *types.Pop, ind int, state string) string {
	var info strings.Builder
	indInfo, exists := pop.IndData[ind]
	if exists {
		fields := []string{
			"birth_year", "sex", "dad", "mom", "lifespan",
			"lat", "lon", "marriage_state", "numbirths",
			"Y_gens", "mt_gens", "min_genealo_gens", "max_genealo_gens",
			"allele_count", "num_blocks", "centromeres", "fitness", "mutations",
		}
		info.WriteString(fmt.Sprintf("%d,%d,%d,", ind, getOrDefault(indInfo, "birth_year", -1), model.FreeParameters["year"]))
		for _, field := range fields {
			info.WriteString(fmt.Sprintf("%d,", getOrDefault(indInfo, field, -1)))
		}
		info.WriteString(state) // Append state ('R' for removed, 'A' for alive)
		info.WriteString("\n")
	}
	return info.String()
}

// getOrDefault retrieves a value from a map with a default fallback
func getOrDefault(data map[string]int, key string, defaultVal int) int {
	if val, ok := data[key]; ok {
		return val
	}
	return defaultVal
}

// writeToFile writes content to a file
func writeToFile(filename string, content string) error {
	file, err := os.OpenFile(filename, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		return fmt.Errorf("could not open file %s: %w", filename, err)
	}
	defer file.Close()

	_, err = file.WriteString(content)
	if err != nil {
		return fmt.Errorf("could not write to file %s: %w", filename, err)
	}
	return err
}

// SaveGenomeMap saves the chromosomes data as an image with rows representing individuals and columns as bit positions.
func SaveGenomeMap(
	chromosomes map[int][][]uint64,
	chromosomeArms map[int]map[int][]int,
	fileName string,
	pixelSize int,
	numbits int,
) error {
	nIndividuals := len(chromosomes)
	imgWidth := numbits*pixelSize + len(chromosomeArms)*4 + 500
	imgHeight := nIndividuals * 2 * pixelSize

	// Create a blank image
	img := image.NewRGBA(image.Rect(0, 0, imgWidth, imgHeight))
	red := color.RGBA{255, 0, 0, 255}
	green := color.RGBA{0, 255, 0, 255}
	black := color.RGBA{0, 0, 0, 255}
	centromereColor := color.RGBA{100, 100, 100, 255}
	chromosomeBoundaryColor := green
	color := red
	counter := 0
	spacer := 0
	for _, chromosomeData := range chromosomes {
		counter += 1
		for genomecopy := 0; genomecopy < 2; genomecopy++ {
			spacer = 0
			sequenceData := chromosomeData[genomecopy]
			startY := (counter*2 + genomecopy) * pixelSize
			for chromosome := 1; chromosome < len(chromosomeArms); chromosome++ {
				// Draw p arm bits
				pStart := chromosomeArms[chromosome][0][0]
				pLength := chromosomeArms[chromosome][0][1]
				for bitPos := pStart; bitPos < pStart+pLength; bitPos++ {
					bitValue := (sequenceData[bitPos/64] >> (bitPos % 64)) & 1
					if bitValue == 0 {
						color = black
					} else {
						color = red
					}
					for y := 0; y < pixelSize; y++ {
						for x := 0; x < pixelSize; x++ {
							img.Set(bitPos*pixelSize+x+spacer, startY+y, color)
						}
					}
				}

				// Draw centromere spacer
				centromereStartX := (pStart + pLength) * pixelSize
				for y := 0; y < pixelSize; y++ {
					for x := 0; x < pixelSize; x++ {
						img.Set(centromereStartX+x+spacer, startY+y, centromereColor)
					}
				}
				spacer += pixelSize

				// Draw q arm bits
				qStart := chromosomeArms[chromosome][1][0]
				qLength := chromosomeArms[chromosome][1][1]
				for bitPos := qStart; bitPos < qStart+qLength; bitPos++ {
					bitValue := (sequenceData[bitPos/64] >> (bitPos % 64)) & 1
					if bitValue == 0 {
						color = black
					} else {
						color = red
					}
					for y := 0; y < pixelSize; y++ {
						for x := 0; x < pixelSize; x++ {
							img.Set(bitPos*pixelSize+x+spacer, startY+y, color)
						}
					}
				}

				// Draw chromosome boundary spacer after q arm
				boundaryStartX := (qStart + qLength) * pixelSize
				for y := 0; y < pixelSize; y++ {
					for i := 0; i < pixelSize; i++ { // Adjust thickness here if needed
						img.Set(boundaryStartX+i+spacer, startY+y, chromosomeBoundaryColor)
					}
				}
				spacer += pixelSize
			}
		}
	}

	// Save the image to file
	file, err := os.Create(fileName)
	if err != nil {
		return err
	}
	defer file.Close()

	err = png.Encode(file, img)
	if err != nil {
		return err
	}

	return nil
}

func calculateMiscStats(indData map[int]map[string]int) (int, int, int, int, int, int, int) {
	var Y, mt, genealo, genetic, alleles, blocks, cents int
	for _, ind := range indData {
		if ind["Y_gens"] > 0 {
			Y += 1
		}
		if ind["mt_gens"] > 0 {
			mt += 1
		}
		if ind["max_genealo_gens"] > -1 {
			genealo++
		}
		if ind["allele_count"] > 0 {
			alleles += ind["allele_count"]
			genetic++
		}
		if ind["num_blocks"] > 0 {
			blocks += ind["num_blocks"]
		}
		if ind["num_centromeres"] > 0 {
			cents += ind["num_centromeres"]
		}
	}
	return Y, mt, genealo, genetic, alleles, blocks, cents
}

func seedCounts(model *types.Model, pop *types.Pop) (int, int, int, int) {

	bitCounts := make([]int, model.FreeParameters["NumBits"])
	totHet, totHomMin, totHomMaj := 0, 0, 0
	seedGenomeRetained := make([]uint64, (model.FreeParameters["NumBits"]+63)/64)

	for _, chromosomePairs := range pop.Chromosomes {
		if len(chromosomePairs) > 0 && len(chromosomePairs[0]) > 0 && len(chromosomePairs[1]) > 0 {
			seedGenomeRetained = bitwiseOR(seedGenomeRetained, chromosomePairs[0])
			seedGenomeRetained = bitwiseOR(seedGenomeRetained, chromosomePairs[1])
			for j := range chromosomePairs[0] {
				b0, b1 := chromosomePairs[0][j], chromosomePairs[1][j]
				bitCounts[j] += countSetBitsSingleVar(b0)
				bitCounts[j] += countSetBitsSingleVar(b1)
				xorBits := b0 ^ b1
				andBits := b0 & b1
				norBits := ^(b0 | b1)
				totHet += countSetBitsSingleVar(xorBits)
				totHomMin += countSetBitsSingleVar(andBits)
				totHomMaj += countSetBitsSingleVar(norBits)
			}
		}
	}
	numbitsRetained := countSetBits(seedGenomeRetained)
	return numbitsRetained, totHet, totHomMin, totHomMaj
}

func calculateFitnessStats(model *types.Model, pop *types.Pop) (numMuts int, totalFitness int) {
	for _, ind := range pop.IndData {
		numMuts += ind["num_mutations"]
		totalFitness += ind["fitness"]
	}
	return numMuts, totalFitness
}

func bitwiseOR(a, b []uint64) []uint64 {
	result := make([]uint64, len(a))
	for i := range a {
		result[i] = a[i] | b[i]
	}
	return result
}

func countSetBitsSingleVar(value uint64) int {
	count := 0
	for value > 0 {
		count += int(value & 1)
		value >>= 1
	}
	return count
}

func countSetBits(data []uint64) int {
	count := 0
	for _, value := range data {
		count += countSetBitsSingleVar(value)
	}
	return count
}
