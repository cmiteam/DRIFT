package events

import (
	"drift/pkg/core"
	"drift/pkg/methods"
	"drift/pkg/modules"
	"math/bits"
)

// Register the seeding modules into the central registry (pkg/modules).
func init() {
	modules.RegisterSeed("single", methods.SeedingSingle)
	modules.RegisterSeed("population", methods.SeedingPopulation)
	modules.RegisterSeed("max_het", methods.SeedingMaxHet)
	// Created founder alleles (TMR4A.md W4): builds A created haplotypes diverged by
	// founder_allele_divergence and gives each founder a germline pool of them.
	// Has no legacy integer code — seed_style="created" selects it by name.
	modules.RegisterSeed("created", methods.SeedingCreated)
}

// Seed dispatches to the seeding module selected by the seed_style parameter
// (legacy integer codes still resolve: 0=single, 1=population, 2=max_het).
func Seed(model *core.Model, pop *core.Pop) { modules.DispatchSeed(model, pop) }

func countSetBits(words []uint64) int {
	count := 0
	for _, word := range words {
		count += bits.OnesCount64(word)
	}
	return count
}
