package events

import (
	"drift/pkg/core"
	"drift/pkg/methods"
	"fmt"
	"math/bits"
)

func Seed(model *core.Model, pop *core.Pop) {
	seedStyle := int(model.Parameters["seed_style"])
	switch seedStyle {
	case core.SeedingSingle:
		methods.SeedingSingle(model, pop)
	case core.SeedingPopulation:
		methods.SeedingPopulation(model, pop)
	case core.SeedingMaxHet:
		methods.SeedingMaxHet(model, pop)
	default:
		fmt.Printf("Unknown seed style %d, using single seed\n", seedStyle)
		methods.SeedingSingle(model, pop)
	}
}

func countSetBits(words []uint64) int {
	count := 0
	for _, word := range words {
		count += bits.OnesCount64(word)
	}
	return count
}
