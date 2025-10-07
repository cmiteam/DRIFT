package necalcs

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"sort"
)

// CalculateSFS generates the site frequency spectrum from allele frequencies
// This function takes the current population and calculates both folded and unfolded SFS
func CalculateSFS(model *core.Model, pop *core.Pop, useAncestralState bool) *core.SFSResult {
	if len(pop.IndData) == 0 {
		return &core.SFSResult{}
	}

	numIndividuals := len(pop.IndData)
	numChromosomes := numIndividuals * 2 // Diploid
	totalBits := model.FreeParameters["genome_bits"]

	// Initialize SFS bins
	foldedSFS := make([]int, (numChromosomes/2)+1) // 0 to n/2
	unfoldedSFS := make([]int, numChromosomes+1)   // 0 to n

	segregatingSites := 0
	singletons := 0
	doubletons := 0
	thetaPiSum := 0.0

	// Count alleles at each position
	for bitPos := 0; bitPos < totalBits; bitPos++ {
		alleleCount := 0

		// Count derived alleles across all individuals
		for id := range pop.IndData {
			if chromosomes, exists := pop.Chromosomes[id]; exists {
				wordIdx := bitPos / 64
				bitIdx := bitPos % 64

				// Count bits from both chromosomes
				if wordIdx < len(chromosomes[0]) {
					if (chromosomes[0][wordIdx]>>bitIdx)&1 == 1 {
						alleleCount++
					}
				}
				if wordIdx < len(chromosomes[1]) {
					if (chromosomes[1][wordIdx]>>bitIdx)&1 == 1 {
						alleleCount++
					}
				}
			}
		}

		// Skip monomorphic sites (all 0 or all 1)
		if alleleCount > 0 && alleleCount < numChromosomes {
			segregatingSites++

			// Add to unfolded SFS
			unfoldedSFS[alleleCount]++

			// Add to folded SFS (use minor allele frequency)
			minorCount := alleleCount
			if alleleCount > numChromosomes/2 {
				minorCount = numChromosomes - alleleCount
			}
			foldedSFS[minorCount]++

			// Count special categories
			if alleleCount == 1 || alleleCount == numChromosomes-1 {
				singletons++
			}
			if alleleCount == 2 || alleleCount == numChromosomes-2 {
				doubletons++
			}

			// Accumulate for Tajima's pi calculation
			p := float64(alleleCount) / float64(numChromosomes)
			thetaPiSum += 2.0 * p * (1.0 - p)
		}
	}

	// Calculate theta statistics
	thetaPi := thetaPiSum / float64(totalBits)
	thetaW := calculateThetaWatterson(segregatingSites, numChromosomes) / float64(totalBits)

	// Calculate Tajima's D
	tajimasD := calculateTajimasD(thetaPi, thetaW, segregatingSites, numChromosomes)

	// Calculate Fu & Li's D*
	fuLiD := calculateFuLiD(singletons, segregatingSites, numChromosomes)

	// Calculate Fay & Wu's H (placeholder - needs ancestral state info)
	fayWuH := 0.0
	if useAncestralState {
		fayWuH = calculateFayWuH(unfoldedSFS, thetaPi, numChromosomes)
	}

	return &core.SFSResult{
		FoldedSFS:   foldedSFS,
		UnfoldedSFS: unfoldedSFS,
		TotalSites:  segregatingSites,
		Singletons:  singletons,
		Doubletons:  doubletons,
		TajimasD:    tajimasD,
		FuLiD:       fuLiD,
		FayWuH:      fayWuH,
		ThetaW:      thetaW,
		ThetaPi:     thetaPi,
		NumSamples:  numChromosomes,
	}
}

// GetSFSForSave extracts key SFS statistics for inclusion in main results file
func GetSFSForSave(sfsResult *core.SFSResult) (int, int, int, float64, float64, float64) {
	if sfsResult == nil {
		return 0, 0, 0, 0.0, 0.0, 0.0
	}

	return sfsResult.TotalSites, sfsResult.Singletons, sfsResult.Doubletons,
		sfsResult.TajimasD, sfsResult.FuLiD, sfsResult.FayWuH
}

// SaveDetailedSFS saves the complete SFS data to a separate file
func SaveDetailedSFS(model *core.Model, sfsResult *core.SFSResult) error {
	if sfsResult == nil || sfsResult.TotalSites == 0 {
		return nil
	}

	filename := fmt.Sprintf("results/%s_sfs_run%d_year%d.csv",
		model.ModelName, model.FreeParameters["run"], model.FreeParameters["year"])

	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("failed to create SFS file: %v", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write header
	header := []string{"FrequencyBin", "FoldedCount", "UnfoldedCount", "Proportion"}
	if err := writer.Write(header); err != nil {
		return err
	}

	// Write folded SFS data
	for i, count := range sfsResult.FoldedSFS {
		if count > 0 {
			proportion := float64(count) / float64(sfsResult.TotalSites)
			unfoldedCount := 0
			if i < len(sfsResult.UnfoldedSFS) {
				unfoldedCount = sfsResult.UnfoldedSFS[i]
			}

			row := []string{
				fmt.Sprintf("%d", i),
				fmt.Sprintf("%d", count),
				fmt.Sprintf("%d", unfoldedCount),
				fmt.Sprintf("%.6f", proportion),
			}
			if err := writer.Write(row); err != nil {
				return err
			}
		}
	}

	return nil
}

// calculateThetaWatterson calculates Watterson's estimator of theta
func calculateThetaWatterson(segregatingSites int, numChromosomes int) float64 {
	if numChromosomes <= 1 {
		return 0.0
	}

	// Calculate harmonic sum: Σ(1/i) for i = 1 to n-1
	harmonicSum := 0.0
	for i := 1; i < numChromosomes; i++ {
		harmonicSum += 1.0 / float64(i)
	}

	if harmonicSum == 0 {
		return 0.0
	}

	return float64(segregatingSites) / harmonicSum
}

// calculateTajimasD calculates Tajima's D neutrality test statistic
func calculateTajimasD(thetaPi, thetaW float64, segregatingSites int, numChromosomes int) float64 {
	if segregatingSites == 0 || numChromosomes <= 1 {
		return 0.0
	}

	n := float64(numChromosomes)
	S := float64(segregatingSites)

	// Calculate variance of Tajima's D
	a1 := 0.0
	a2 := 0.0
	for i := 1; i < numChromosomes; i++ {
		fi := float64(i)
		a1 += 1.0 / fi
		a2 += 1.0 / (fi * fi)
	}

	b1 := (n + 1.0) / (3.0 * (n - 1.0))
	b2 := (2.0 * (n*n + n + 3.0)) / (9.0 * n * (n - 1.0))
	c1 := b1 - (1.0 / a1)
	c2 := b2 - ((n + 2.0) / (a1 * n)) + (a2 / (a1 * a1))
	e1 := c1 / a1
	e2 := c2 / (a1*a1 + a2)

	variance := e1*S + e2*S*(S-1.0)

	if variance <= 0 {
		return 0.0
	}

	return (thetaPi - thetaW) / math.Sqrt(variance)
}

// calculateFuLiD calculates Fu & Li's D* statistic
func calculateFuLiD(singletons int, segregatingSites int, numChromosomes int) float64 {
	if segregatingSites == 0 || numChromosomes <= 1 {
		return 0.0
	}

	n := float64(numChromosomes)
	S := float64(segregatingSites)
	eta := float64(singletons)

	// Calculate a1 (harmonic sum)
	a1 := 0.0
	for i := 1; i < numChromosomes; i++ {
		a1 += 1.0 / float64(i)
	}

	if a1 == 0 {
		return 0.0
	}

	// Expected number of singletons under neutrality
	expectedSingletons := S / a1

	// Simplified variance calculation for Fu & Li's D*
	vD := (n/(n-1.0))*a1 - 1.0
	uD := vD + ((n-1.0)/n)*a1*a1

	if uD <= 0 {
		return 0.0
	}

	return (eta - expectedSingletons) / math.Sqrt(uD*S)
}

// calculateFayWuH calculates Fay & Wu's H statistic (requires unfolded SFS)
func calculateFayWuH(unfoldedSFS []int, thetaPi float64, numChromosomes int) float64 {
	if len(unfoldedSFS) == 0 || numChromosomes <= 1 {
		return 0.0
	}

	n := float64(numChromosomes)
	thetaH := 0.0

	// Calculate theta_H = Σ(i²ξᵢ) * 2/(n(n-1))
	for i := 1; i < len(unfoldedSFS); i++ {
		if i < numChromosomes {
			thetaH += float64(i*i) * float64(unfoldedSFS[i])
		}
	}
	thetaH *= 2.0 / (n * (n - 1.0))

	// H = θ_π - θ_H
	return thetaPi - thetaH
}

func CalculateBinnedSFS(model *core.Model, pop *core.Pop) []int {
	if len(pop.IndData) == 0 {
		return make([]int, 101) // Return empty bins
	}

	numIndividuals := len(pop.IndData)
	numChromosomes := numIndividuals * 2
	totalBits := model.FreeParameters["genome_bits"]

	// Initialize 101 bins for frequencies 0.00, 0.01, 0.02, ... 1.00
	binnedSFS := make([]int, 101)

	// Count alleles at each position and bin by frequency
	for bitPos := 0; bitPos < totalBits; bitPos++ {
		alleleCount := 0

		// Count derived alleles (1-bits) across all individuals
		for id := range pop.IndData {
			if chromosomes, exists := pop.Chromosomes[id]; exists {
				wordIdx := bitPos / 64
				bitIdx := bitPos % 64

				// Count bits from both chromosomes
				if wordIdx < len(chromosomes[0]) {
					if (chromosomes[0][wordIdx]>>bitIdx)&1 == 1 {
						alleleCount++
					}
				}
				if wordIdx < len(chromosomes[1]) {
					if (chromosomes[1][wordIdx]>>bitIdx)&1 == 1 {
						alleleCount++
					}
				}
			}
		}

		// Calculate frequency and determine bin
		frequency := float64(alleleCount) / float64(numChromosomes)

		// Determine which bin this frequency falls into
		// frequency 0.00-0.004999 -> bin 0
		// frequency 0.005-0.014999 -> bin 1
		// frequency 0.015-0.024999 -> bin 2
		// etc.
		binIndex := int(frequency*100 + 0.5) // Round to nearest percent
		if binIndex > 100 {
			binIndex = 100 // Cap at 100 for frequency = 1.0
		}

		binnedSFS[binIndex]++
	}

	return binnedSFS
}

// SaveBinnedSFS saves the binned SFS data to a standardized format file
func SaveBinnedSFS(model *core.Model, binnedSFS []int) error {
	filename := fmt.Sprintf("results/%s_binned_sfs_run%d_year%d.csv",
		model.ModelName, model.FreeParameters["run"], model.FreeParameters["year"])

	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("failed to create binned SFS file: %v", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write header
	header := []string{"Frequency", "Count", "Proportion"}
	if err := writer.Write(header); err != nil {
		return err
	}

	// Calculate total sites for proportion calculation
	totalSites := 0
	for _, count := range binnedSFS {
		totalSites += count
	}

	// Write binned data
	for i, count := range binnedSFS {
		frequency := float64(i) / 100.0 // Convert bin index to frequency
		proportion := 0.0
		if totalSites > 0 {
			proportion = float64(count) / float64(totalSites)
		}

		row := []string{
			fmt.Sprintf("%.2f", frequency),
			fmt.Sprintf("%d", count),
			fmt.Sprintf("%.6f", proportion),
		}
		if err := writer.Write(row); err != nil {
			return err
		}
	}

	return nil
}

// GetSFSMoments calculates summary statistics from binned SFS for main results
func GetSFSMoments(binnedSFS []int) (float64, float64, float64, float64) {
	if len(binnedSFS) == 0 {
		return 0.0, 0.0, 0.0, 0.0
	}

	totalSites := 0
	for _, count := range binnedSFS {
		totalSites += count
	}

	if totalSites == 0 {
		return 0.0, 0.0, 0.0, 0.0
	}

	// Calculate moments of the frequency distribution
	mean := 0.0
	variance := 0.0
	skewness := 0.0

	// First pass: calculate mean
	for i, count := range binnedSFS {
		frequency := float64(i) / 100.0
		weight := float64(count) / float64(totalSites)
		mean += frequency * weight
	}

	// Second pass: calculate variance and skewness
	for i, count := range binnedSFS {
		frequency := float64(i) / 100.0
		weight := float64(count) / float64(totalSites)
		diff := frequency - mean
		variance += diff * diff * weight
		if variance > 0 {
			skewness += diff * diff * diff * weight
		}
	}

	// Normalize skewness
	if variance > 0 {
		skewness = skewness / (variance * math.Sqrt(variance))
	}

	// Calculate proportion of sites at extremes (< 5% or > 95%)
	extremesProportion := float64(binnedSFS[0]+binnedSFS[1]+binnedSFS[2]+binnedSFS[3]+binnedSFS[4]) / float64(totalSites)
	extremesProportion += float64(binnedSFS[96]+binnedSFS[97]+binnedSFS[98]+binnedSFS[99]+binnedSFS[100]) / float64(totalSites)

	return mean, variance, skewness, extremesProportion
}

// SaveSFSTimeSeries saves SFS data with frequency bins as rows and years as columns
func SaveSFSTimeSeries(model *core.Model, pop *core.Pop) error {
	if pop.SFSHistory == nil || len(pop.SFSHistory) == 0 {
		return nil
	}

	filename := fmt.Sprintf("results/%s_sfs_timeseries_run%d.csv",
		model.ModelName, model.FreeParameters["run"])

	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("failed to create SFS time series file: %v", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Get sorted years for consistent column order
	years := make([]int, 0, len(pop.SFSHistory))
	for year := range pop.SFSHistory {
		years = append(years, year)
	}
	sort.Ints(years)

	// Create header: Frequency, Year1, Year2, Year3, ...
	header := []string{"Frequency"}
	for _, year := range years {
		header = append(header, fmt.Sprintf("Year%d", year))
	}
	if err := writer.Write(header); err != nil {
		return err
	}

	// Write data rows - one row per frequency bin
	for freqBin := 0; freqBin <= 100; freqBin++ {
		frequency := float64(freqBin) / 100.0
		row := []string{fmt.Sprintf("%.2f", frequency)}

		// Add count for this frequency bin across all years
		for _, year := range years {
			count := 0
			if sfs, exists := pop.SFSHistory[year]; exists && freqBin < len(sfs) {
				count = sfs[freqBin]
			}
			row = append(row, fmt.Sprintf("%d", count))
		}

		if err := writer.Write(row); err != nil {
			return err
		}
	}

	return nil
}
