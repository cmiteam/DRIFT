package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
)

func SeedingMaxHet(model *core.Model, pop *core.Pop) {
	model.FreeParameters["seed"] = 1

	totalBits := model.FreeParameters["genome_bits"]

	// Initialize all chromosomes first
	for id := range pop.IndData {
		pop.Chromosomes[id] = [][]uint64{
			make([]uint64, (totalBits+63)/64),
			make([]uint64, (totalBits+63)/64),
		}
		for i := range pop.Chromosomes[id][0] {
			pop.Chromosomes[id][0][i] = ^uint64(0)
		}

		// Initialize centromeres
		var blankCentromeres [2][]uint64
		blankCentromeres[0] = make([]uint64, (len(model.ChromosomeArms)+63)/64)
		blankCentromeres[1] = make([]uint64, (len(model.ChromosomeArms)+63)/64)
		pop.Centromeres[id] = blankCentromeres
		for chrom := 0; chrom < len(model.ChromosomeArms); chrom++ {
			idx := chrom / 64
			bit := chrom % 64
			pop.Centromeres[id][0][idx] |= (1 << bit)
			pop.Centromeres[id][1][idx] |= (1 << bit)
		}

		// Initialize tracking
		pop.IndData[id][individual.YGens] = 0
		pop.IndData[id][individual.MtGens] = 0
		pop.IndData[id][individual.MinGenealoGens] = 0
		pop.IndData[id][individual.MaxGenealoGens] = 0
		pop.IndData[id][individual.AlleleCount] = 0
		pop.IndData[id][individual.NumCentromeres] = 0
	}

	// For each bit position in genome

	// Update allele counts for each individual
	for id := range pop.IndData {
		count := utils.CountSetBits(pop.Chromosomes[id][0]) + utils.CountSetBits(pop.Chromosomes[id][1])
		pop.IndData[id][individual.AlleleCount] = count
	}
}
