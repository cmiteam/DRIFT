package analysis

import (
	"drift/pkg/core"
	"math/rand"
)

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
