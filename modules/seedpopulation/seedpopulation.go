package seedpopulation

import (
	"drift/modules/individual"
	"drift/types"
	"fmt"
	"math"
	"math/bits"
	"math/rand"
)

// SeedThePopulation chooses a random seed individual and sets up their genetic data

func SeedThePopulation(model *types.Model, pop *types.Pop) {
	seedStyle := int(model.Parameters["seed_style"])
	switch seedStyle {
	case 0:
		seedSingleIndividual(model, pop)
	case 1:
		seedPopulation(model, pop)
		//	case 2:
		//		seedMaxHet(model, pop)
	default:
		fmt.Printf("Unknown seed style %d, using single seed\n", seedStyle)
		seedSingleIndividual(model, pop)
	}
}

func seedSingleIndividual(model *types.Model, pop *types.Pop) {
	seed := chooseRandomSeed(model, pop)

	if seed == -1 {
		return
	}

	model.FreeParameters["seed"] = seed
	//	fmt.Println("   Seed:", model.FreeParameters["seed"])

	// Create chromosomes
	pop.Chromosomes[seed] = [][]uint64{
		make([]uint64, (model.FreeParameters["genome_bits"]+63)/64),
		make([]uint64, (model.FreeParameters["genome_bits"]+63)/64),
	}

	// Set all bits to 1 in chromosomes
	for i := range pop.Chromosomes[seed][0] {
		pop.Chromosomes[seed][0][i] = ^uint64(0)
	}
	for i := range pop.Chromosomes[seed][1] {
		pop.Chromosomes[seed][1][i] = ^uint64(0)
	}
	// Create centromeres
	var blankCentromeres [2][]uint64
	blankCentromeres[0] = make([]uint64, len(model.ChromosomeArms)+63/64)
	blankCentromeres[1] = make([]uint64, len(model.ChromosomeArms)+63/64)
	pop.Centromeres[seed] = blankCentromeres

	// Set all bits to 1 in centromeres
	for chrom := 0; chrom < len(model.ChromosomeArms); chrom++ {
		idx := chrom / 64
		bit := chrom % 64
		pop.Centromeres[seed][0][idx] |= (1 << bit)
		pop.Centromeres[seed][1][idx] |= (1 << bit)
	}
	pop.IndData[seed][individual.YGens] = 0
	pop.IndData[seed][individual.MtGens] = 0
	pop.IndData[seed][individual.MaxGenealoGens] = 0
	pop.IndData[seed][individual.MinGenealoGens] = 0
	pop.IndData[seed][individual.AlleleCount] = model.FreeParameters["genome_bits"] * 2
	pop.IndData[seed][individual.NumCentromeres] = countSetBits(pop.Centromeres[seed][0])
}

func chooseRandomSeed(model *types.Model, pop *types.Pop) int {
	year := model.FreeParameters["year"]
	matureMales := []int{}
	for id, data := range pop.IndData {
		age := year - data[individual.BirthYear]
		if data[individual.Sex] == 0 && age >= int(model.Parameters["maturity"]) {
			matureMales = append(matureMales, id)
		}
	}
	if len(matureMales) == 0 {
		return -1
	}
	return matureMales[rand.Intn(len(matureMales))]
}

func seedPopulation(model *types.Model, pop *types.Pop) {
	hetLevel := model.Parameters["init_heterozygosity"]
	if hetLevel <= 0 {
		fmt.Println("   No genomic seeding (init_heterozygosity = 0)")
		return
	}

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
			altFreq := -0.05 * math.Log(rand.Float64()) // Exponential decay
			if altFreq > 0.5 {
				altFreq = 1 - altFreq // Cap at 50%
			}
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
		count := countSetBits(pop.Chromosomes[id][0]) + countSetBits(pop.Chromosomes[id][1])
		pop.IndData[id][individual.AlleleCount] = count
	}

}

func countSetBits(words []uint64) int {
	count := 0
	for _, word := range words {
		count += bits.OnesCount64(word)
	}
	return count
}
