package utils

import (
	"drift/pkg/core"
)

func InheritMutations(pop *core.Pop, genomemask []uint64, parent int, child int, copy int) {
	if _, exists := pop.IndMutations[child]; !exists {
		pop.IndMutations[child] = map[int][]int{0: {}, 1: {}}
	}

	for _, mutationID := range pop.IndMutations[parent][0] {
		mutation := pop.MutationPool[mutationID]

		// Check which strand this bit comes from in the genomemask
		wordIdx := mutation.Position / 64
		bitIdx := mutation.Position % 64
		inheritedStrand := (genomemask[wordIdx] >> bitIdx) & 1

		if int(inheritedStrand) == 0 {
			pop.IndMutations[child][copy] = append(pop.IndMutations[child][copy], mutationID)
			pop.MutationCount++
			pop.MutationPool[mutationID] = mutation
		}
	}

	for _, mutationID := range pop.IndMutations[parent][1] {
		mutation := pop.MutationPool[mutationID]
		wordIdx := mutation.Position / 64
		bitIdx := mutation.Position % 64
		inheritedStrand := (genomemask[wordIdx] >> bitIdx) & 1

		if int(inheritedStrand) == 1 {
			pop.IndMutations[child][copy] = append(pop.IndMutations[child][copy], mutationID)
			mutation.Count++
			pop.MutationPool[mutationID] = mutation
		}
	}
}

func GenerateNewMutations(model *core.Model, pop *core.Pop, ind int) {

	numNewMutations := RandPoisson(model.Parameters["mu"])

	for i := 0; i < numNewMutations; i++ {
		model.FreeParameters["mutID"]++
		mutationID := model.FreeParameters["mutID"]
		position := RandIntn(int(model.FreeParameters["genome_bits"]))
		mutationEffect := 0.0
		isMutationNonNeutral := RandFloat64()
		if isMutationNonNeutral >= model.Parameters["f_neutral"] {
			mutationEffect = WeibullRandom(model.Parameters["shape"], model.Parameters["scale"]) / model.Parameters["Weibull_adj"]
			isMutationDeleterious := RandFloat64()
			if isMutationDeleterious > model.Parameters["f_beneficial"] {
				mutationEffect = -mutationEffect
			}
		}
		pop.MutationHist[int(mutationEffect*model.Parameters["mu_scale_factor"])]++

		strand := RandIntn(2)
		if pop.IndMutations[ind] == nil {
			pop.IndMutations[ind] = make(map[int][]int)
		}
		if pop.IndMutations[ind][strand] == nil {
			pop.IndMutations[ind][strand] = []int{}
		}
		pop.IndMutations[ind][strand] = append(pop.IndMutations[ind][strand], mutationID)
		pop.MutationPool[mutationID] = core.Mutation{
			Id:        mutationID,
			Position:  position,
			Effect:    mutationEffect,
			Origin:    ind,
			Count:     1,
			Dominance: 0,
		}
	}
}

func addMutation(pop *core.Pop, ind int, position int, copy int, value int) {
	if pop.IndMutations[ind] == nil {
		pop.IndMutations[ind] = make(map[int][]int)
	}
	if pop.IndMutations[ind][position] == nil {
		pop.IndMutations[ind][position] = []int{}
	}
	pop.IndMutations[ind][position][copy] = value
}

func CountFitnessAndMutations(pop *core.Pop, child int) (int, float64) {
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
