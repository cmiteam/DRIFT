package simulation

import (
	"drift/modules/coalescence"
	"drift/modules/individual"
	"drift/pkg/utils"
	"drift/pkg/core"
	"drift/pkg/analysis"
	"encoding/csv"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"math/bits"
	"os"
	"sort"
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
		"nCents", "nAlleles", "nBlocks", "TotFitness", "nMuts", "PercSeedGenoRet", "AvSeedGenoCov", "AvHet", "HomMin", "HomMaj",
		"numSegSites", "allelesLost", "allelesFixed", "avMinorFreq", "freqVariance", "TajimaPi", "WattersonsTheta", "Fis", "Ne",
		"numSegSitesSFS", "singletons", "doubletons", "TajimasD", "FuLiD", "FayWuH", "sfsMean", "sfsVar", "sfsSkew", "sfsExtremes",
		"LD_short", "LD_medium", "LD_long", "HighLDPairs", "LDDecaySlope",
		"meanCoalTime", "maxCoalTime", "coalPairs", "foundMRCA",
		"tmrcaMean", "tmrcaMax", "tmrcaVariance", "estimatedNe", "tmrcaSampleSize",
		"yadamID", "yadamGensBack", "yadamBirthYear",
	}
	if err := writer.Write(headers); err != nil {
		return fmt.Errorf("failed to write headers: %v", err)
	}
	return nil
}

// Save writes the current simulation state to a CSV file
func Save(model *core.Model, pop *core.Pop, animManager *core.AnimationsContainer) {

	filename := fmt.Sprintf("results/%s_results.csv", model.ModelName)
	numInds := len(pop.IndData)
	var YDescends, mtDescends, genealoDescends, geneticDescends, numAlleles, numBlocks, numCentromeres, numMutations, popFitness int
	var percSeedGenomeRetained, avSeedGenomeCoverage float64
	var totHet, totHomMin, totHomMaj, numbitsRetained int
	var numSegSites, allelesLost, allelesFixed int
	var avMinorFreq, freqVariance float64
	var tajimaPi, wattersonsTheta, fis, ne float64
	var numSegSitesSFS, singletons, doubletons int
	var tajimasD, fuLiD, fayWuH float64
	var sfsMean, sfsVar, sfsSkew, sfsExtremes float64
	var ldShort, ldMedium, ldLong, ldDecaySlope float64
	var highLDPairs int
	var meanCoalTime, maxCoalTime float64
	var coalescencePairs, foundMRCA int
	var tmrcaMean, tmrcaMax, tmrcaVariance, estimatedNe float64
	var tmrcaSampleSize int
	var yadamID, yadamGensBack, yadamBirthYear int

	if model.Parameters["track_DNA"] == 1 {
		YDescends, mtDescends, genealoDescends, geneticDescends, numAlleles, numBlocks, numCentromeres = calculateMiscStats(&pop.IndData)
		numbitsRetained, totHet, totHomMin, totHomMaj = seedCounts(model, pop)
		percSeedGenomeRetained = float64(numbitsRetained) / float64(model.FreeParameters["genome_bits"]) * 100
		avSeedGenomeCoverage = 0
	}

	if model.Parameters["track_mutations"] == 1 {
		numMutations, popFitness = calculateFitnessStats(model, pop)
	}

	if model.Parameters["track_drift"] == 1 {
		frequencies := trackAlleleFrequencies(model, pop)
		storeAlleleFrequencies(model, pop, frequencies)
		numSegSites, allelesLost, allelesFixed, avMinorFreq, freqVariance = calculateAlleleStats(frequencies)
		tajimaPi, wattersonsTheta, fis = calculatePopGenStats(model, pop, frequencies)
	}

	// Calculate SFS statistics if enabled
	if model.Parameters["track_SFS"] == 1 {
		sfsResult := analysis.CalculateSFS(model, pop, false) // Use false for folded SFS
		numSegSitesSFS, singletons, doubletons, tajimasD, fuLiD, fayWuH = analysis.GetSFSForSave(sfsResult)
		binnedSFS := analysis.CalculateBinnedSFS(model, pop)
		sfsMean, sfsVar, sfsSkew, sfsExtremes = analysis.GetSFSMoments(binnedSFS)

		// Save detailed SFS to separate file if requested
		if model.Parameters["save_detailed_SFS"] == 1 {
			binnedSFS := analysis.CalculateBinnedSFS(model, pop)
			//necalcs.SaveBinnedSFS(model, binnedSFS)
			if pop.SFSHistory == nil {
				pop.SFSHistory = make(map[int][]int)
			}
			pop.SFSHistory[model.FreeParameters["year"]] = binnedSFS
		}
	}
	if model.Parameters["track_LD"] == 1 {
		ldSummary := analysis.CalculateLDSummary(model, pop)
		ldShort, ldMedium, ldLong, highLDPairs, ldDecaySlope = analysis.GetLDForSave(ldSummary)
	}

	if model.Parameters["track_coalescence"] == 1 {
		yadamResult := coalescence.FindYAdam(model, pop)
		mteveResult := coalescence.FindMtEve(model, pop)

		if yadamResult.FoundYAdam {
			fmt.Printf("Y-Adam found: ID %d, born year %d\n",
				yadamResult.YAdamID, yadamResult.YAdamBirthYear)
		}

		if mteveResult.FoundMtEve {
			fmt.Printf("Mt-Eve found: ID %d, born year %d\n",
				mteveResult.MtEveID, mteveResult.MtEveBirthYear)
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
		fmt.Sprintf("%d", numSegSites),
		fmt.Sprintf("%d", allelesLost),
		fmt.Sprintf("%d", allelesFixed),
		fmt.Sprintf("%.4f", avMinorFreq),
		fmt.Sprintf("%.6f", freqVariance),
		fmt.Sprintf("%.4f", tajimaPi),
		fmt.Sprintf("%.4f", wattersonsTheta),
		fmt.Sprintf("%.4f", fis),
		fmt.Sprintf("%.2f", ne),

		fmt.Sprintf("%d", numSegSitesSFS),
		fmt.Sprintf("%d", singletons),
		fmt.Sprintf("%d", doubletons),
		fmt.Sprintf("%.4f", tajimasD),
		fmt.Sprintf("%.4f", fuLiD),
		fmt.Sprintf("%.4f", fayWuH),
		fmt.Sprintf("%.4f", sfsMean),
		fmt.Sprintf("%.4f", sfsVar),
		fmt.Sprintf("%.4f", sfsSkew),
		fmt.Sprintf("%.4f", sfsExtremes),

		fmt.Sprintf("%.4f", ldShort),
		fmt.Sprintf("%.4f", ldMedium),
		fmt.Sprintf("%.4f", ldLong),
		fmt.Sprintf("%d", highLDPairs),
		fmt.Sprintf("%.6f", ldDecaySlope),

		fmt.Sprintf("%.2f", meanCoalTime),
		fmt.Sprintf("%.2f", maxCoalTime),
		fmt.Sprintf("%d", coalescencePairs),
		fmt.Sprintf("%d", foundMRCA),
		fmt.Sprintf("%.2f", tmrcaMean),
		fmt.Sprintf("%.2f", tmrcaMax),
		fmt.Sprintf("%.4f", tmrcaVariance),
		fmt.Sprintf("%.1f", estimatedNe),
		fmt.Sprintf("%d", tmrcaSampleSize),
		// Y-Adam analysis fields
		fmt.Sprintf("%d", yadamID),
		fmt.Sprintf("%d", yadamGensBack),
		fmt.Sprintf("%d", yadamBirthYear),
	}
	writer.Write(data)

	fmt.Printf("\n   Year: %d  n: %d  b: %d  m: %d c: %d g: %d g: %d",
		model.FreeParameters["year"],
		model.FreeParameters["last_pop_size"],
		pop.Tracking["births"],
		pop.Tracking["marriages"],
		pop.Tracking["cull_deaths"],
		genealoDescends,
		geneticDescends,
	)

	pop.Tracking["births"] = 0
	pop.Tracking["deaths"] = 0
	pop.Tracking["marriages"] = 0
	pop.Tracking["random_deaths"] = 0
	pop.Tracking["cull_deaths"] = 0
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

func calculateMiscStats(indData *map[int][]int) (int, int, int, int, int, int, int) {
	var Y, mt, genealo, genetic, alleles, blocks, cents int
	for _, ind := range *indData {
		if ind[individual.YGens] > 0 {
			Y += 1
		}
		if ind[individual.MtGens] > 0 {
			mt += 1
		}
		if ind[individual.MaxGenealoGens] > -1 {
			genealo++
		}
		if ind[individual.AlleleCount] > 0 {
			alleles += ind[individual.AlleleCount]
			genetic++
		}
		if ind[individual.NumBlocks] > 0 {
			blocks += ind[individual.NumBlocks]
		}
		if ind[individual.NumCentromeres] > 0 {
			cents += ind[individual.NumCentromeres]
		}
	}
	return Y, mt, genealo, genetic, alleles, blocks, cents
}

func seedCounts(model *core.Model, pop *core.Pop) (int, int, int, int) {

	bitCounts := make([]int, model.FreeParameters["genome_bits"])
	totHet, totHomMin, totHomMaj := 0, 0, 0
	seedGenomeRetained := make([]uint64, (model.FreeParameters["genome_bits"]+63)/64)

	for _, chromosomePairs := range pop.Chromosomes {
		if len(chromosomePairs) > 0 && len(chromosomePairs[0]) > 0 && len(chromosomePairs[1]) > 0 {
			seedGenomeRetained = utils.BitwiseOR(seedGenomeRetained, chromosomePairs[0])
			seedGenomeRetained = utils.BitwiseOR(seedGenomeRetained, chromosomePairs[1])
			for j := range chromosomePairs[0] {
				b0, b1 := chromosomePairs[0][j], chromosomePairs[1][j]
				bitCounts[j] += utils.CountSetBitsSingleVar(b0)
				bitCounts[j] += utils.CountSetBitsSingleVar(b1)
				xorBits := b0 ^ b1
				andBits := b0 & b1
				norBits := ^(b0 | b1)
				totHet += utils.CountSetBitsSingleVar(xorBits)
				totHomMin += utils.CountSetBitsSingleVar(andBits)
				totHomMaj += utils.CountSetBitsSingleVar(norBits)
			}
		}
	}
	numbitsRetained := utils.CountSetBits(seedGenomeRetained)
	return numbitsRetained, totHet, totHomMin, totHomMaj
}

func calculateFitnessStats(model *core.Model, pop *core.Pop) (numMuts int, totalFitness int) {
	for _, ind := range pop.IndData {
		numMuts += ind[individual.NumMutations]
		totalFitness += ind[individual.Fitness]
	}
	return numMuts, totalFitness
}

func trackAlleleFrequencies(model *core.Model, pop *core.Pop) []float64 {
	totalBits := model.FreeParameters["genome_bits"]
	numWords := (totalBits + 63) / 64
	popSize := len(pop.IndData)

	// Create bit count array for each bit position
	bitCounts := make([]int, totalBits)

	// Process each word position across all individuals
	for wordIdx := 0; wordIdx < numWords; wordIdx++ {
		// For each bit in this word
		for bitIdx := 0; bitIdx < 64; bitIdx++ {
			bitPos := wordIdx*64 + bitIdx
			if bitPos >= totalBits {
				break
			}

			mask := uint64(1) << bitIdx
			count := 0

			// Count this specific bit across all chromosomes
			for id := range pop.IndData {
				if chromosomes, exists := pop.Chromosomes[id]; exists {
					for copy := 0; copy < 2 && copy < len(chromosomes); copy++ {
						if wordIdx < len(chromosomes[copy]) {
							if chromosomes[copy][wordIdx]&mask != 0 {
								count++
							}
						}
					}
				}
			}
			bitCounts[bitPos] = count
		}
	}

	// Convert to frequencies
	frequencies := make([]float64, totalBits)
	totalChromosomes := float64(popSize * 2)

	for i, count := range bitCounts {
		frequencies[i] = float64(count) / totalChromosomes
	}

	return frequencies
}

func storeAlleleFrequencies(model *core.Model, pop *core.Pop, frequencies []float64) {
	year := model.FreeParameters["year"]
	scaled := make([]int16, len(frequencies))

	for i, freq := range frequencies {
		scaled[i] = int16(freq * 10000) // Scale and convert
	}

	pop.AlleleFreqs[year] = scaled
}

func calculateNe(pop *core.Pop, startYear, endYear int) float64 {
	deltaT := float64(endYear - startYear)

	freqs1 := pop.AlleleFreqs[startYear]
	freqs2 := pop.AlleleFreqs[endYear]

	var variance float64
	validLoci := 0

	for i := range freqs1 {
		p1 := float64(freqs1[i]) / 10000.0
		p2 := float64(freqs2[i]) / 10000.0
		deltaP := p2 - p1
		variance += deltaP * deltaP
		validLoci++
	}

	variance /= float64(validLoci)
	return deltaT / (2.0 * variance)
}

func calculatePopGenStats(model *core.Model, pop *core.Pop, frequencies []float64) (float64, float64, float64) {
	n := 2.0 * float64(len(pop.IndData)) // Number of chromosomes (diploid)

	// Fast Tajima's π and segregating sites calculation
	tajimaPi := 0.0
	segregatingSites := 0

	for _, p := range frequencies {
		if p > 0 && p < 1 {
			segregatingSites++
			tajimaPi += 2 * p * (1 - p) // This is the per-site π
		}
	}

	// Average over all sites
	if len(frequencies) > 0 {
		tajimaPi /= float64(len(frequencies))
	}

	// Watterson's θ
	harmonicSum := 0.0
	for i := 1; i < int(n); i++ {
		harmonicSum += 1.0 / float64(i)
	}

	wattersonsTheta := 0.0
	if harmonicSum > 0 && len(frequencies) > 0 {
		wattersonsTheta = float64(segregatingSites) / harmonicSum / float64(len(frequencies))
	}

	// F_IS calculation using optimized heterozygosity
	observedHet := calculateObservedHeterozygosity(model, pop)
	expectedHet := tajimaPi // Expected het = Tajima's π

	fis := 0.0
	if expectedHet > 0 {
		fis = (expectedHet - observedHet) / expectedHet
	}

	return tajimaPi, wattersonsTheta, fis
}

func calculateObservedHeterozygosity(model *core.Model, pop *core.Pop) float64 {
	totalBits := model.FreeParameters["genome_bits"]
	numWords := (totalBits + 63) / 64
	totalIndividuals := len(pop.IndData)

	if totalIndividuals == 0 {
		return 0.0
	}

	totalHetBits := 0

	for id := range pop.IndData {
		if chromosomes, exists := pop.Chromosomes[id]; exists && len(chromosomes) >= 2 {
			for wordIdx := 0; wordIdx < numWords && wordIdx < len(chromosomes[0]) && wordIdx < len(chromosomes[1]); wordIdx++ {
				// XOR gives us heterozygous positions in one operation
				hetWord := chromosomes[0][wordIdx] ^ chromosomes[1][wordIdx]

				// Count set bits in the XOR result
				if wordIdx == numWords-1 {
					// Handle partial last word
					bitsInLastWord := totalBits % 64
					if bitsInLastWord == 0 {
						bitsInLastWord = 64
					}
					mask := (uint64(1) << bitsInLastWord) - 1
					hetWord &= mask
				}

				totalHetBits += bits.OnesCount64(hetWord)
			}
		}
	}

	return float64(totalHetBits) / (float64(totalIndividuals) * float64(totalBits))
}

func calculateAlleleStats(frequencies []float64) (int, int, int, float64, float64) {
	numSegSites := 0
	allelesLost := 0
	allelesFixed := 0
	minorFreqSum := 0.0
	freqVariance := 0.0
	mean := 0.0

	// Single pass through frequencies
	for _, p := range frequencies {
		mean += p
		if p > 0.0 && p < 1.0 {
			numSegSites++
			// Minor allele frequency
			minorFreq := p
			if p > 0.5 {
				minorFreq = 1.0 - p
			}
			minorFreqSum += minorFreq
		} else if p == 0.0 {
			allelesLost++
		} else if p == 1.0 {
			allelesFixed++
		}
	}

	mean /= float64(len(frequencies))

	// Calculate variance in second pass
	for _, p := range frequencies {
		diff := p - mean
		freqVariance += diff * diff
	}
	freqVariance /= float64(len(frequencies))

	avMinorFreq := 0.0
	if numSegSites > 0 {
		avMinorFreq = minorFreqSum / float64(numSegSites)
	}

	return numSegSites, allelesLost, allelesFixed, avMinorFreq, freqVariance
}

func calculateNeFromHeterozygosityDecline(model *core.Model, currentHet float64) float64 {
	if len(model.HetHistory) < 3 {
		return -1
	}

	// Calculate slope of heterozygosity over time using linear regression
	n := len(model.HetHistory)
	if n < 5 { // Need at least 5 time points for stability
		return -1
	}

	// Use recent time points only (last 10 or all if less than 10)
	startIdx := 0
	if n > 10 {
		startIdx = n - 10
	}

	sumX := 0.0
	sumY := 0.0
	sumXY := 0.0
	sumX2 := 0.0
	count := 0.0

	for i := startIdx; i < n; i++ {
		x := float64(model.TimeHistory[i] - model.TimeHistory[startIdx]) // Relative time
		y := model.HetHistory[i]

		sumX += x
		sumY += y
		sumXY += x * y
		sumX2 += x * x
		count++
	}

	if count < 3 || sumX2 <= 0 {
		return -1
	}

	// Linear regression slope: (n*ΣXY - ΣX*ΣY) / (n*ΣX² - (ΣX)²)
	denominator := count*sumX2 - sumX*sumX
	if math.Abs(denominator) < 1e-12 {
		return -1
	}

	slope := (count*sumXY - sumX*sumY) / denominator

	if slope >= 0 { // No decline in heterozygosity
		return -1
	}

	// Ne ≈ -He / (2 × slope) but with correction for discrete generations
	ne := -currentHet / (2.0 * slope)

	// Apply reasonable bounds
	if ne < 10 || ne > 50000 {
		return -1
	}

	return ne
}

func calculateNeIfPossibleImproved(pop *core.Pop) float64 {
	if len(pop.AlleleFreqs) < 2 {
		return -1 // Not enough data points
	}

	// Get all available years
	var years []int
	for year := range pop.AlleleFreqs {
		years = append(years, year)
	}
	sort.Ints(years)

	// Try different time intervals for stability
	bestNe := -1.0

	// Try intervals of 1, 5, 10, and 20 years if available
	intervals := []int{1, 5, 10, 20}

	for _, interval := range intervals {
		if len(years) < interval+1 {
			continue
		}

		// Use the most recent data points with this interval
		endYear := years[len(years)-1]
		startYear := endYear - interval

		// Find the closest year to startYear
		closestStart := -1
		minDiff := math.MaxInt32
		for _, year := range years {
			diff := utils.Abs(year - startYear)
			if diff < minDiff {
				minDiff = diff
				closestStart = year
			}
		}

		if closestStart != -1 {
			ne := calculateNeWithFiltering(pop, closestStart, endYear)
			if ne > 0 && ne < 1000000 { // Reasonable bounds
				if bestNe == -1.0 || ne < bestNe { // Prefer smaller, more conservative estimates
					bestNe = ne
				}
			}
		}
	}

	return bestNe
}

func calculateObservedHeterozygosityOptimized(model *core.Model, pop *core.Pop) float64 {
	totalBits := model.FreeParameters["genome_bits"]
	numWords := (totalBits + 63) / 64
	totalIndividuals := len(pop.IndData)

	if totalIndividuals == 0 {
		return 0.0
	}

	totalHetBits := 0

	for id := range pop.IndData {
		if chromosomes, exists := pop.Chromosomes[id]; exists && len(chromosomes) >= 2 {
			for wordIdx := 0; wordIdx < numWords && wordIdx < len(chromosomes[0]) && wordIdx < len(chromosomes[1]); wordIdx++ {
				// XOR gives us heterozygous positions in one operation
				hetWord := chromosomes[0][wordIdx] ^ chromosomes[1][wordIdx]

				// Handle partial last word
				if wordIdx == numWords-1 {
					bitsInLastWord := totalBits % 64
					if bitsInLastWord == 0 {
						bitsInLastWord = 64
					}
					mask := (uint64(1) << bitsInLastWord) - 1
					hetWord &= mask
				}

				totalHetBits += bits.OnesCount64(hetWord)
			}
		}
	}

	return float64(totalHetBits) / (float64(totalIndividuals) * float64(totalBits))
}

func calculateNeWithFiltering(pop *core.Pop, startYear, endYear int) float64 {
	deltaT := float64(endYear - startYear)
	if deltaT <= 0 {
		return -1
	}

	freqs1 := pop.AlleleFreqs[startYear]
	freqs2 := pop.AlleleFreqs[endYear]

	var totalVariance float64
	validLoci := 0
	minFreqThreshold := 0.01 // Only consider loci with MAF > X%
	maxVariance := 0.01      // Cap individual locus variance

	for i := range freqs1 {
		p1 := float64(freqs1[i]) / 10000.0
		p2 := float64(freqs2[i]) / 10000.0

		// Only include loci that are polymorphic and not too rare
		minFreq1 := math.Min(p1, 1-p1)
		minFreq2 := math.Min(p2, 1-p2)

		if minFreq1 > minFreqThreshold && minFreq2 > minFreqThreshold {
			deltaP := p2 - p1
			variance := deltaP * deltaP

			// Cap variance to avoid extreme values
			if variance > maxVariance {
				variance = maxVariance
			}

			totalVariance += variance
			validLoci++
		}
	}

	if validLoci < 10 { // Need sufficient loci
		return -1
	}

	meanVariance := totalVariance / float64(validLoci)

	if meanVariance <= 0 {
		return -1
	}

	ne := deltaT / (2.0 * meanVariance)

	// Apply reasonable bounds
	if ne < 10 || ne > 100000 {
		return -1
	}

	return ne
}

func getCoalescenceForSave(results []*core.CoalescenceResult) (float64, float64, int, int) {
	if results == nil || len(results) == 0 {
		return 0.0, 0.0, 0, 0
	}

	var totalTime float64
	var maxTime float64
	foundCount := 0

	for _, result := range results {
		if result.FoundMRCA && result.CoalescenceTime > 0 {
			totalTime += float64(result.CoalescenceTime)
			foundCount++

			if float64(result.CoalescenceTime) > maxTime {
				maxTime = float64(result.CoalescenceTime)
			}
		}
	}

	meanTime := 0.0
	if foundCount > 0 {
		meanTime = totalTime / float64(foundCount)
	}

	return meanTime, maxTime, len(results), foundCount
}
