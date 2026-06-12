package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"fmt"
	"math"
	"math/rand"
)

func SeedingPopulation(model *core.Model, pop *core.Pop) {
	hetLevel := model.Parameters["init_heterozygosity"]
	if hetLevel <= 0 {
		fmt.Println("   No genomic seeding (init_heterozygosity = 0)")
		return
	}
	model.FreeParameters["seed"] = 1

	totalBits := model.FreeParameters["genome_bits"]
	variableSites := 0

	// Initialize all chromosomes first
	for id := range pop.IndData {
		pop.Chromosomes[id] = [][]uint64{
			make([]uint64, (totalBits+63)/64),
			make([]uint64, (totalBits+63)/64),
		}
		// Initialize centromeres
		var blankCentromeres [2][]uint64
		blankCentromeres[0] = make([]uint64, (len(model.ChromosomeArms)+63)/64)
		blankCentromeres[1] = make([]uint64, (len(model.ChromosomeArms)+63)/64)
		pop.Centromeres[id] = blankCentromeres

		// Initialize tracking
		pop.IndData[id][individual.YGens] = 0
		pop.IndData[id][individual.MtGens] = 0
		pop.IndData[id][individual.MinGenealoGens] = 0
		pop.IndData[id][individual.MaxGenealoGens] = 0
		pop.IndData[id][individual.AlleleCount] = 0
		pop.IndData[id][individual.NumCentromeres] = 0
	}

	// For each bit position in genome
	for bitPosition := 0; bitPosition < totalBits; bitPosition++ {
		// Is that bit going to be variable?
		if rand.Float64() < hetLevel {
			variableSites++

			// Pick the frequency of the alt allele (exponential function, >0 to <= 0.5)
			altFreq := -0.05 * math.Log(rand.Float64()) * 0.1 // Exponential decay
			//	if altFreq > 0.5 {
			//		altFreq = 1 - altFreq // Cap at 50%
			//	}
			if altFreq < 0.005 {
				altFreq = 0.005 // Minimum 0.5%
			}

			// Scan through each individual
			for id := range pop.IndData {
				// For each chromosome
				for chromosome := 0; chromosome < 2; chromosome++ {
					// If rand < probability of bit being set, set the bit to 1
					if rand.Float64() < altFreq {
						wordIdx := bitPosition / 64
						bitIdx := bitPosition % 64
						pop.Chromosomes[id][chromosome][wordIdx] |= (1 << bitIdx)
					}
				}
			}
		}
	}

	// Update allele counts for each individual
	for id := range pop.IndData {
		count := utils.CountSetBits(pop.Chromosomes[id][0]) + utils.CountSetBits(pop.Chromosomes[id][1])
		pop.IndData[id][individual.AlleleCount] = count
	}
}
