package utils

import (
	"drift/pkg/core"
	"sort"
)

// InheritMutations passes a parent's de-novo mutations to a child gamete through
// the same per-birth recombination mask (genomemask) that meiosis applies to the
// founder Chromosomes bitfield, so mutations that share an arm co-segregate and
// recombination breaks their linkage with distance (roadmap §1, linkage-aware
// positions). The genomemask bit at a mutation's Position decides whether it is
// inherited on this gamete.
//
// LINKAGE MODEL (opt-in, default byte-identical). Historically this used the
// OPPOSITE mask polarity from meiosis: a strand-0 mutation was inherited where
// mask==0, whereas meiosis takes the copy-0 bitfield allele where mask==1. So a
// de-novo mutation and a founder allele at the SAME position anti-segregated —
// the pool and the bitfield were inconsistent users of the same coordinate space
// (the §6h "de-novo mutations bypass the bitfield" finding). With
// linkage_model="arm" the polarity is aligned with meiosis, so a mutation at
// position P now follows the founder bit at P and the two systems co-segregate.
//
// The default linkage_model="legacy" preserves the original (inverted) polarity
// exactly — InheritMutations consumes no RNG, so both paths leave the RNG stream
// untouched and the default run stays byte-identical (the §6h -validate gate).
//
// MUTATION-POOL REFERENCE COUNTING (mutation_count_model; DEFAULT "refcount";
// roadmap §6f). Mutation.Count is the reference count death.go uses to garbage-
// collect a lineage (Count-- per dead copy; delete at Count<=0). The original
// ("legacy") increment was ASYMMETRIC: only the strand-1 branch below bumped
// Mutation.Count (the strand-0 branch bumped the vestigial global pop.MutationCount
// instead), so Count under-counted true copies and the COMMON/FIXED lineages —
// those with the most inheritance/death events — were GC'd first, WHILE STILL
// CARRIED. That had two downstream effects: (1) the engine's own fitness under-
// counted load, because CountFitnessAndMutations skips ids missing from the pool,
// so fixed load / the Muller's ratchet were invisible (§6f finding P4); and (2) a
// still-carried but GC'd mutation, when inherited here, resolved pop.MutationPool
// [id] to the zero Mutation (Position 0) and was therefore placed by the mask bit
// at position 0 rather than its true locus — which also DISTORTED the neutral SFS
// (artificial position-0 linkage; §6f characterization: legacy Tajima's D −0.89 /
// Ne/N 0.19 / SFS χ²/dof 18.9 vs refcount −0.66 / 0.31 / 4.8).
//
// The DEFAULT mutation_count_model="refcount" makes BOTH branches increment
// Mutation.Count, so Count equals the living-copy count (creation=1, +1 per
// inherited copy, -1 per dead copy) and a lineage is deleted only when genuinely
// lost — fixed load is observable and inherited positions stay correct. This is the
// corrected neutral baseline the §6h -validate gate now characterizes.
//
// COMPATIBILITY. InheritMutations consumes NO RNG under either model, and the
// demography (births/deaths/N) is unchanged, so a fixed seed stays reproducible.
// The opt-out "legacy" reproduces the original asymmetry exactly (for byte-for-byte
// comparison against pre-§6f runs); it is strictly opt-in now that refcount is the
// default.
func InheritMutations(model *core.Model, pop *core.Pop, genomemask []uint64, parent int, child int, copy int) {
	if _, exists := pop.IndMutations[child]; !exists {
		pop.IndMutations[child] = map[int][]int{0: {}, 1: {}}
	}

	// Which mask bit value keeps a strand-0 / strand-1 mutation on this gamete.
	// Legacy: strand0 where mask==0, strand1 where mask==1 (inverted vs meiosis).
	// Arm: swap them so mutation@P tracks the founder bitfield allele@P (meiosis
	// takes copy0 where mask==1), reconciling the pool with the bitfield.
	var keep0, keep1 uint64 = 0, 1
	if model.StringParam("linkage_model", "legacy") == "arm" {
		keep0, keep1 = 1, 0
	}
	// Under "refcount" (the DEFAULT) the strand-0 branch also increments
	// Mutation.Count so the reference count is symmetric (see the doc comment
	// above); "legacy" opts back into the original asymmetry. Only the explicit
	// "legacy" opts out — an unset/unrecognized value is the corrected default.
	refcount := model.StringParam("mutation_count_model", "refcount") != "legacy"

	for _, mutationID := range pop.IndMutations[parent][0] {
		mutation := pop.MutationPool[mutationID]

		// Check which strand this bit comes from in the genomemask
		wordIdx := mutation.Position / 64
		bitIdx := mutation.Position % 64
		inheritedStrand := (genomemask[wordIdx] >> bitIdx) & 1

		if inheritedStrand == keep0 {
			pop.IndMutations[child][copy] = append(pop.IndMutations[child][copy], mutationID)
			pop.MutationCount++
			if refcount {
				mutation.Count++
			}
			pop.MutationPool[mutationID] = mutation
		}
	}

	for _, mutationID := range pop.IndMutations[parent][1] {
		mutation := pop.MutationPool[mutationID]
		wordIdx := mutation.Position / 64
		bitIdx := mutation.Position % 64
		inheritedStrand := (genomemask[wordIdx] >> bitIdx) & 1

		if inheritedStrand == keep1 {
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
// The position is drawn uniformly by default, or — when the class carries a
// regional rate map (the `mutation_rate_map` param, or a per-class `ratemap=`
// override; roadmap §1) — in proportion to per-region rate multipliers so
// mutations concentrate in hotspots. The uniform path is the strict no-op.
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
			// Position draw: uniform by default (strict no-op — the identical single
			// RandIntn(genome_bits) draw as before variable rates existed), or a
			// rate-weighted draw when this class carries a regional rate map (§1).
			var position int
			if class.RateMap == nil {
				position = RandIntn(genomeBits)
			} else {
				position = class.RateMap.Draw()
			}
			mutationEffect := 0.0
			isMutationNonNeutral := RandFloat64()
			if isMutationNonNeutral >= class.FNeutral {
				// Magnitude draw. "weibull" (default) is the byte-identical historical
				// path (WeibullRandom / Weibull_adj). "gamma" (opt-in, roadmap §6f)
				// draws from Gamma(shape, mean/shape) so the mean |effect| is GammaMean
				// and the shape matches the human deleterious DFE (alpha ~ 0.2). Only
				// this branch differs; the neutral coin above and the beneficial coin
				// below are drawn identically, so a non-gamma class is unchanged.
				if class.DFEModel == "gamma" {
					mutationEffect = GammaRandom(class.GammaShape, class.GammaMean/class.GammaShape)
				} else {
					mutationEffect = WeibullRandom(class.Shape, class.Scale) / class.WeibullAdj
				}
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
