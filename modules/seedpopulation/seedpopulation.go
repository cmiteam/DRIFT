package seedpopulation

import (
	"drift/modules/individual"
	"drift/types"
	"fmt"
	"math/rand"
)

// SeedThePopulation chooses a random seed individual and sets up their genetic data

func SeedThePopulation(model *types.Model, pop *types.Pop) {
	seed := chooseRandomSeed(model, pop)

	if seed == -1 {
		return
	}

	model.FreeParameters["seed"] = seed
	fmt.Println("   Seed:", model.FreeParameters["seed"])

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
	pop.IndData[seed][individual.Lat] = 0
	pop.IndData[seed][individual.Lon] = 0

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

func countSetBitsSingleVar(n uint64) int {
	// Counts the number of bits set to 1 in n
	count := 0
	for n > 0 {
		count += int(n & 1)
		n >>= 1
	}
	return count
}

func countSetBits(bits []uint64) int {
	count := 0
	for _, b := range bits {
		count += countSetBitsSingleVar(b)
	}
	return count
}
