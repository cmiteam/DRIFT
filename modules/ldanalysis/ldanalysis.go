package ldanalysis

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"math"
	"math/rand"
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

// CalculateBinnedSFS generates SFS in standardized frequency bins (0.00, 0.01, 0.02, ... 1.00)
// This makes it easy to compare across different population sizes and runs
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

// CalculatePairwiseLD calculates linkage disequilibrium between two genomic positions
func CalculatePairwiseLD(model *core.Model, pop *core.Pop, site1, site2 int) *core.LDResult {
	if len(pop.IndData) == 0 || site1 == site2 {
		return nil
	}

	// Count haplotype frequencies: 00, 01, 10, 11
	var n00, n01, n10, n11 int
	totalHaplotypes := 0

	// For each individual, examine both chromosomes
	for id := range pop.IndData {
		if chromosomes, exists := pop.Chromosomes[id]; exists {
			// Get alleles at both sites for both chromosomes
			alleles1 := make([]int, 2)
			alleles2 := make([]int, 2)

			for chrom := 0; chrom < 2; chrom++ {
				wordIdx1 := site1 / 64
				bitIdx1 := site1 % 64
				wordIdx2 := site2 / 64
				bitIdx2 := site2 % 64

				// Get allele at site1
				if wordIdx1 < len(chromosomes[chrom]) {
					alleles1[chrom] = int((chromosomes[chrom][wordIdx1] >> bitIdx1) & 1)
				}

				// Get allele at site2
				if wordIdx2 < len(chromosomes[chrom]) {
					alleles2[chrom] = int((chromosomes[chrom][wordIdx2] >> bitIdx2) & 1)
				}

				// Count haplotype combinations
				if alleles1[chrom] == 0 && alleles2[chrom] == 0 {
					n00++
				} else if alleles1[chrom] == 0 && alleles2[chrom] == 1 {
					n01++
				} else if alleles1[chrom] == 1 && alleles2[chrom] == 0 {
					n10++
				} else if alleles1[chrom] == 1 && alleles2[chrom] == 1 {
					n11++
				}
				totalHaplotypes++
			}
		}
	}

	if totalHaplotypes == 0 {
		return nil
	}

	// Calculate haplotype frequencies
	//f00 := float64(n00) / float64(totalHaplotypes)
	f01 := float64(n01) / float64(totalHaplotypes)
	f10 := float64(n10) / float64(totalHaplotypes)
	f11 := float64(n11) / float64(totalHaplotypes)

	// Calculate allele frequencies
	p1 := f10 + f11 // Frequency of allele 1 at site 1
	p2 := f01 + f11 // Frequency of allele 1 at site 2
	q1 := 1.0 - p1  // Frequency of allele 0 at site 1
	q2 := 1.0 - p2  // Frequency of allele 0 at site 2

	// Calculate linkage disequilibrium coefficient D
	D := f11 - (p1 * p2)

	// Calculate D' (normalized D)
	var DPrime float64
	if D >= 0 {
		Dmax := math.Min(p1*q2, q1*p2)
		if Dmax > 0 {
			DPrime = D / Dmax
		}
	} else {
		Dmax := math.Min(p1*p2, q1*q2)
		if Dmax > 0 {
			DPrime = -D / Dmax
		}
	}

	// Calculate r² (correlation coefficient)
	var rSquared float64
	denominator := p1 * q1 * p2 * q2
	if denominator > 0 {
		rSquared = (D * D) / denominator
	}

	// Calculate chi-square test statistic
	var chiSq float64
	expected00 := float64(totalHaplotypes) * q1 * q2
	expected01 := float64(totalHaplotypes) * q1 * p2
	expected10 := float64(totalHaplotypes) * p1 * q2
	expected11 := float64(totalHaplotypes) * p1 * p2

	if expected00 > 0 {
		chiSq += math.Pow(float64(n00)-expected00, 2) / expected00
	}
	if expected01 > 0 {
		chiSq += math.Pow(float64(n01)-expected01, 2) / expected01
	}
	if expected10 > 0 {
		chiSq += math.Pow(float64(n10)-expected10, 2) / expected10
	}
	if expected11 > 0 {
		chiSq += math.Pow(float64(n11)-expected11, 2) / expected11
	}

	// Simple p-value approximation (chi-square with 1 df)
	pValue := math.Exp(-chiSq / 2)

	distance := int(math.Abs(float64(site2 - site1)))

	return &core.LDResult{
		Site1:    site1,
		Site2:    site2,
		Distance: distance,
		RSquared: rSquared,
		DPrime:   math.Abs(DPrime), // Take absolute value
		ChiSq:    chiSq,
		PValue:   pValue,
	}
}

// LDSummary contains summary statistics for LD analysis
type LDSummary struct {
	MeanLD_0_10    float64 // Mean r² for distances 0-10 bp
	MeanLD_11_50   float64 // Mean r² for distances 11-50 bp
	MeanLD_51_100  float64 // Mean r² for distances 51-100 bp
	MeanLD_101_500 float64 // Mean r² for distances 101-500 bp
	MeanLD_500plus float64 // Mean r² for distances >500 bp
	HighLDPairs    int     // Number of pairs with r² > 0.8
	TotalPairs     int     // Total pairs analyzed
	LDDecaySlope   float64 // Slope of LD decay curve
}

// CalculateLDSummary analyzes linkage disequilibrium patterns with distance binning
func CalculateLDSummary(model *core.Model, pop *core.Pop) *LDSummary {
	if len(pop.IndData) == 0 {
		return &LDSummary{}
	}

	// Get list of segregating sites to analyze
	segregatingSites := getSegregatingSites(model, pop)
	if len(segregatingSites) < 2 {
		return &LDSummary{}
	}

	// Sample pairs to avoid O(n²) explosion
	maxPairs := 10000 // Limit total pairs analyzed
	sampleRate := 1.0
	totalPossiblePairs := len(segregatingSites) * (len(segregatingSites) - 1) / 2
	if totalPossiblePairs > maxPairs {
		sampleRate = float64(maxPairs) / float64(totalPossiblePairs)
	}

	// Distance bins for LD analysis
	var sum0_10, sum11_50, sum51_100, sum101_500, sum500plus float64
	var count0_10, count11_50, count51_100, count101_500, count500plus int
	var highLDPairs int
	var totalPairs int

	// Analyze sampled pairs
	for i := 0; i < len(segregatingSites); i++ {
		for j := i + 1; j < len(segregatingSites); j++ {
			// Sample pairs based on sampling rate
			if sampleRate < 1.0 && rand.Float64() > sampleRate {
				continue
			}

			site1 := segregatingSites[i]
			site2 := segregatingSites[j]

			ldResult := CalculatePairwiseLD(model, pop, site1, site2)
			if ldResult == nil {
				continue
			}

			totalPairs++
			distance := ldResult.Distance
			rSquared := ldResult.RSquared

			// Count high LD pairs
			if rSquared > 0.8 {
				highLDPairs++
			}

			// Bin by distance
			if distance <= 10 {
				sum0_10 += rSquared
				count0_10++
			} else if distance <= 50 {
				sum11_50 += rSquared
				count11_50++
			} else if distance <= 100 {
				sum51_100 += rSquared
				count51_100++
			} else if distance <= 500 {
				sum101_500 += rSquared
				count101_500++
			} else {
				sum500plus += rSquared
				count500plus++
			}
		}
	}

	// Calculate means for each distance bin
	summary := &LDSummary{
		HighLDPairs: highLDPairs,
		TotalPairs:  totalPairs,
	}

	if count0_10 > 0 {
		summary.MeanLD_0_10 = sum0_10 / float64(count0_10)
	}
	if count11_50 > 0 {
		summary.MeanLD_11_50 = sum11_50 / float64(count11_50)
	}
	if count51_100 > 0 {
		summary.MeanLD_51_100 = sum51_100 / float64(count51_100)
	}
	if count101_500 > 0 {
		summary.MeanLD_101_500 = sum101_500 / float64(count101_500)
	}
	if count500plus > 0 {
		summary.MeanLD_500plus = sum500plus / float64(count500plus)
	}

	// Calculate LD decay slope (simple linear regression)
	summary.LDDecaySlope = calculateLDDecaySlope([]float64{summary.MeanLD_0_10, summary.MeanLD_11_50, summary.MeanLD_51_100, summary.MeanLD_101_500, summary.MeanLD_500plus})

	return summary
}

// getSegregatingSites returns positions of all segregating sites
func getSegregatingSites(model *core.Model, pop *core.Pop) []int {
	totalBits := model.FreeParameters["genome_bits"]
	numChromosomes := len(pop.IndData) * 2
	var segregatingSites []int

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

		// Include if segregating (not fixed at 0 or 1)
		if alleleCount > 0 && alleleCount < numChromosomes {
			segregatingSites = append(segregatingSites, bitPos)
		}
	}

	return segregatingSites
}

// calculateLDDecaySlope estimates slope of LD decay with distance
func calculateLDDecaySlope(means []float64) float64 {
	if len(means) < 2 {
		return 0.0
	}

	// Simple slope calculation: (last - first) / (bins - 1)
	validMeans := 0
	for _, mean := range means {
		if mean > 0 {
			validMeans++
		}
	}

	if validMeans < 2 {
		return 0.0
	}

	// Calculate slope from first to last valid mean
	firstValid := -1
	lastValid := -1
	for i, mean := range means {
		if mean > 0 {
			if firstValid == -1 {
				firstValid = i
			}
			lastValid = i
		}
	}

	if firstValid == lastValid {
		return 0.0
	}

	return (means[lastValid] - means[firstValid]) / float64(lastValid-firstValid)
}

// GetLDForSave extracts key LD statistics for main results file
func GetLDForSave(ldSummary *LDSummary) (float64, float64, float64, int, float64) {
	if ldSummary == nil {
		return 0.0, 0.0, 0.0, 0, 0.0
	}

	return ldSummary.MeanLD_0_10, ldSummary.MeanLD_51_100, ldSummary.MeanLD_500plus,
		ldSummary.HighLDPairs, ldSummary.LDDecaySlope
}
