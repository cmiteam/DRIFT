package utils

import (
	"drift/pkg/core"
	"math"
	"sort"
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

	// Dominance coefficient h for new mutations. Default 0.5 (additive) when the
	// param is absent — models created before fitness_dominance existed don't
	// carry the key, and a bare map lookup would silently yield h=0 (recessive).
	dominancePct := 50
	if h, ok := model.Parameters["fitness_dominance"]; ok {
		dominancePct = int(math.Round(h * 100))
	}

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
			Id:       mutationID,
			Position: position,
			Effect:   mutationEffect,
			Origin:   ind,
			Count:    1,
			// Dominance coefficient h, stored as an int percentage (h*100): 0 =
			// fully recessive, 100 = fully dominant, 50 = additive/codominant. It
			// weights a mutation's Effect in the heterozygous state (see
			// CountFitnessAndMutations). A single global h (fitness_dominance) is
			// used for now; a per-mutation DFE-linked h would just vary this line.
			Dominance: dominancePct,
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

// CountFitnessAndMutations computes an individual's mutation copy count and its
// genotype fitness load with per-locus dominance.
//
// pop.IndMutations[ind] is keyed by inherited strand (0/1); each value is the
// list of mutation ids on that strand. Every distinct mutation id is treated as
// one biallelic locus (infinite-sites): a mutation carried on BOTH strands is
// homozygous derived and contributes its full Effect; carried on ONE strand it
// is heterozygous and contributes h*Effect, where h is the mutation's dominance
// coefficient (Mutation.Dominance/100; h=0 recessive, 0.5 additive, 1 dominant).
// Per-locus contributions are combined per the fitness_model param (see
// CombineFitness, fitness.go). Homozygosity here means the SAME mutation
// on both strands (identity by descent) — the mechanism that makes recessive
// load express under bottlenecks/founder events, when IBD rises.
//
// Returns the total number of mutation copies carried (both strands) and the
// combined relative fitness (per fitness_model). Iteration is sorted (strands, then mutation
// ids) so the float sum is reproducible under a fixed seed.
func CountFitnessAndMutations(pop *core.Pop, ind int, model *core.Model) (int, float64) {
	strands := pop.IndMutations[ind]
	if len(strands) == 0 {
		return 0, 0.0
	}

	// Count copies (1 = heterozygous, 2 = homozygous) of each mutation id across
	// the up-to-two strands, and the total copies carried.
	copies := make(map[int]int)
	numCopies := 0
	strandKeys := make([]int, 0, len(strands))
	for s := range strands {
		strandKeys = append(strandKeys, s)
	}
	sort.Ints(strandKeys)
	for _, s := range strandKeys {
		for _, mutationID := range strands[s] {
			copies[mutationID]++
			numCopies++
		}
	}

	// Sum per-locus fitness contributions in sorted mutation-id order (float
	// addition is not associative; keep it deterministic under a fixed seed).
	ids := make([]int, 0, len(copies))
	for id := range copies {
		ids = append(ids, id)
	}
	sort.Ints(ids)

	contribs := make([]float64, 0, len(ids))
	for _, id := range ids {
		mutation, found := pop.MutationPool[id]
		if !found {
			continue
		}
		if copies[id] >= 2 {
			// Homozygous derived: full effect.
			contribs = append(contribs, mutation.Effect)
		} else {
			// Heterozygous: effect weighted by the dominance coefficient h.
			contribs = append(contribs, float64(mutation.Dominance)/100.0*mutation.Effect)
		}
	}
	return numCopies, CombineFitness(contribs, model)
}
