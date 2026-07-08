package utils

import (
	"drift/pkg/core"
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

// GenerateNewMutations adds de-novo mutations to individual ind, drawing each
// mutation class in the configured spectrum (roadmap §1) independently.
//
// For each class it draws Poisson(class.Rate) new mutations and, per mutation,
// the same sequence as the historical single-class kernel: position, a
// neutral/non-neutral coin, (if non-neutral) a Weibull magnitude + a
// deleterious/beneficial coin, then a strand. The mutation records its class
// index and target size so downstream stats can partition by class.
//
// NEUTRAL NO-OP (preserved for the §6h -validate harness): with mutation_classes
// unset the spectrum is a single point class whose rate and DFE are the model's
// global mu / Weibull / dominance params, so the loop draws exactly one
// Poisson(mu) and then the identical per-mutation RNG stream as before the
// spectrum existed — strictly-neutral runs stay byte-identical.
func GenerateNewMutations(model *core.Model, pop *core.Pop, ind int) {
	classes := MutationClasses(model)
	muScale := model.Parameters["mu_scale_factor"]
	genomeBits := int(model.FreeParameters["genome_bits"])

	for classIdx, class := range classes {
		numNewMutations := RandPoisson(class.Rate)
		for i := 0; i < numNewMutations; i++ {
			model.FreeParameters["mutID"]++
			mutationID := model.FreeParameters["mutID"]
			position := RandIntn(genomeBits)
			mutationEffect := 0.0
			isMutationNonNeutral := RandFloat64()
			if isMutationNonNeutral >= class.FNeutral {
				mutationEffect = WeibullRandom(class.Shape, class.Scale) / class.WeibullAdj
				isMutationDeleterious := RandFloat64()
				if isMutationDeleterious > class.FBeneficial {
					mutationEffect = -mutationEffect
				}
			}
			pop.MutationHist[int(mutationEffect*muScale)]++

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
				// CountFitnessAndMutations). Inherited from the class (which defaults
				// to the global fitness_dominance); a per-mutation DFE-linked h would
				// just vary this line.
				Dominance: class.Dominance,
				// Class spectrum bookkeeping (roadmap §1): the class index and its
				// target size, so downstream analysis can filter by mutation class.
				Class: classIdx,
				Size:  class.Size,
			}
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
