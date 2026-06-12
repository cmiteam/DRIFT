package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"math/rand"
)

func SeedingSingle(model *core.Model, pop *core.Pop) {
	seed := chooseRandomSeed(model, pop)

	if seed == -1 {
		return
	}

	model.FreeParameters["seed"] = seed

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
	pop.IndData[seed][individual.NumCentromeres] = utils.CountSetBits(pop.Centromeres[seed][0])
}

func chooseRandomSeed(model *core.Model, pop *core.Pop) int {
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
