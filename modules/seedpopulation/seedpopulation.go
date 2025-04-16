package seedpopulation

import (
	"drift/types"
	"fmt"
	"math/rand"
)

// SeedThePopulation chooses a random seed individual and sets up their genetic data

func SeedThePopulation(model *types.Model, pop *types.Pop, year int) {
	seed := chooseRandomSeed(model, pop, year)

	if seed == -1 {
		return
	}

	model.FreeParameters["seed"] = seed
	fmt.Println("   Seed:", model.FreeParameters["seed"])

	// Create chromosomes
	pop.Chromosomes[seed] = [][]uint64{
		make([]uint64, (model.FreeParameters["NumBits"]+63)/64),
		make([]uint64, (model.FreeParameters["NumBits"]+63)/64),
	}

	// Set all bits to 1 in chromosomes
	for i := range pop.Chromosomes[seed][0] {
		pop.Chromosomes[seed][0][i] = ^uint64(0)
	}
	for i := range pop.Chromosomes[seed][1] {
		pop.Chromosomes[seed][1][i] = ^uint64(0)
	}
	// Create centromeres
	pop.Centromeres[seed] = []uint64{0, 0}
	// Set all bits to 1 in centromeres
	for i := 1; i <= len(model.ChromosomeArms); i++ {
		pop.Centromeres[seed][0] = setBit(pop.Centromeres[seed][0], i)
		pop.Centromeres[seed][1] = setBit(pop.Centromeres[seed][1], i)
	}
	pop.IndData[seed]["Y_gens"] = 0
	pop.IndData[seed]["mt_gens"] = 0
	pop.IndData[seed]["max_genealo_gens"] = 0
	pop.IndData[seed]["min_genealo_gens"] = 0
	pop.IndData[seed]["allele_count"] = model.FreeParameters["NumBits"] * 2
	pop.IndData[seed]["num_centomeres"] = countSetBitsSingleVar(pop.Centromeres[seed][0])
	pop.IndData[seed]["num_centomeres"] += countSetBitsSingleVar(pop.Centromeres[seed][1])

}

func chooseRandomSeed(model *types.Model, pop *types.Pop, year int) int {

	matureMales := []int{}
	for id, data := range pop.IndData {
		age := year - data["birth_year"]
		if data["sex"] == 0 && age >= int(model.Parameters["maturity"]) {
			matureMales = append(matureMales, id)
		}
	}
	if len(matureMales) == 0 {
		return -1
	}
	return matureMales[rand.Intn(len(matureMales))]
}

func setBit(value uint64, bitPosition int) uint64 {
	return value | (1 << bitPosition)
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
