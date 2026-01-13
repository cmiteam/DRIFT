package analysis

import (
	"drift/pkg/core"
	"math"
)

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

	// Calculate r-squared (correlation coefficient)
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

func CalculateThetaWatterson(segregatingSites int, numChromosomes int) float64 {
	if numChromosomes <= 1 {
		return 0.0
	}

	// Calculate harmonic sum: S(1/i) for i = 1 to n-1
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

func CalculateTajimasD(thetaPi, thetaW float64, segregatingSites int, numChromosomes int) float64 {
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

func CalculateFuLiD(singletons int, segregatingSites int, numChromosomes int) float64 {
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

func CalculateFayWuH(unfoldedSFS []int, thetaPi float64, numChromosomes int) float64 {
	if len(unfoldedSFS) == 0 || numChromosomes <= 1 {
		return 0.0
	}

	n := float64(numChromosomes)
	thetaH := 0.0

	// Calculate theta_H = S(i-squared??) * 2/(n(n-1))
	for i := 1; i < len(unfoldedSFS); i++ {
		if i < numChromosomes {
			thetaH += float64(i*i) * float64(unfoldedSFS[i])
		}
	}
	thetaH *= 2.0 / (n * (n - 1.0))

	// H = ?_p - ?_H
	return thetaPi - thetaH
}
