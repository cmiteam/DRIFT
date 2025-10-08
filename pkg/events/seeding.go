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
	case 0:
		methods.SeedingSingle(model, pop)
	case 1:
		methods.SeedingPopulation(model, pop)
	case 2:
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
