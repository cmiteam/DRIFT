package mutation

import (
	"drift/types"
	"gonum.org/v1/gonum/stat/distuv"
	"math"
	"math/rand"
	"time"
)

func InheritMutations(pop *types.Pop, genomemask []uint64, parent int, child int, copy int) {
	if _, exists := pop.IndMutations[child]; !exists {
		pop.IndMutations[child] = map[int][]int{
			0: {},
			1: {},
		}
	}
	for _, mutationID := range pop.IndMutations[parent][0] {
		mutation := pop.MutationPool[mutationID]
		// determine the bit position of the mutation
		// TO DO: don't hard code "1000000", this applies to the default genome only and will cause problems when multiplier != 1
		mutationBin := mutation.Position / 1000000
		inheritedStrand := (genomemask[mutationBin/64] >> (mutationBin % 64)) & 1
		if int(inheritedStrand) == 0 {
			pop.IndMutations[child][copy] = append(pop.IndMutations[child][copy], mutationID)
			pop.MutationCount++
			pop.MutationPool[mutationID] = mutation
		}
	}
	for _, mutationID := range pop.IndMutations[parent][1] {
		mutation := pop.MutationPool[mutationID]
		mutationBin := mutation.Position / 1000000 // See comment above
		inheritedStrand := (genomemask[mutationBin/64] >> (mutationBin % 64)) & 1
		if int(inheritedStrand) == 1 {
			pop.IndMutations[child][copy] = append(pop.IndMutations[child][copy], mutationID)
			mutation.Count++
			pop.MutationPool[mutationID] = mutation
		}
	}
}

func GenerateNewMutations(model *types.Model, pop *types.Pop, ind int) {

	rand.Seed(time.Now().UnixNano())
	poisson := distuv.Poisson{Lambda: model.Parameters["mu"]}
	numNewMutations := int(poisson.Rand())

	for i := 0; i < numNewMutations; i++ {
		model.FreeParameters["mutID"]++
		mutationID := model.FreeParameters["mutID"]
		position := rand.Intn(int(model.FreeParameters["numbits"]))
		mutationEffect := 0.0
		isMutationNonNeutral := rand.Float64()
		if isMutationNonNeutral >= model.Parameters["f_neutral"] {
			mutationEffect = weibullRandom(model.Parameters["shape"], model.Parameters["scale"]) / model.Parameters["Weibull_adj"]
			isMutationDeleterious := rand.Float64()
			if isMutationDeleterious > model.Parameters["f_beneficial"] {
				mutationEffect = -mutationEffect
			}
		}
		pop.MutationHist[int(mutationEffect*model.Parameters["mu_scale_factor"])]++

		strand := rand.Intn(2)
		if pop.IndMutations[ind] == nil {
			pop.IndMutations[ind] = make(map[int][]int)
		}
		if pop.IndMutations[ind][strand] == nil {
			pop.IndMutations[ind][strand] = []int{}
		}
		pop.IndMutations[ind][strand] = append(pop.IndMutations[ind][strand], mutationID)
		pop.MutationPool[mutationID] = types.Mutation{
			Id:        mutationID,
			Position:  position,
			Effect:    mutationEffect,
			Origin:    ind,
			Count:     1,
			Dominance: 0,
		}
	}
}

func weibullRandom(shape, scale float64) float64 {
	u := rand.Float64()
	return scale * math.Pow(-math.Log(u), 1/shape)
}

func addMutation(pop *types.Pop, ind int, position int, copy int, value int) {
	if pop.IndMutations[ind] == nil {
		pop.IndMutations[ind] = make(map[int][]int)
	}
	if pop.IndMutations[ind][position] == nil {
		pop.IndMutations[ind][position] = []int{}
	}
	pop.IndMutations[ind][position][copy] = value
}

func CountFitnessAndMutations(pop *types.Pop, child int) (int, float64) {
	numMutations := 0
	fitnessEffect := 0.0
	if len(pop.IndMutations[child]) > 0 {
		for mutationID, _ := range pop.IndMutations[child] {
			numMutations++
			if mutation, found := pop.MutationPool[mutationID]; found {
				fitnessEffect += mutation.Effect
			}
		}
	}
	//    print (" FE:", fitnessEffect, ":", numMutations)
	return numMutations, fitnessEffect
}
