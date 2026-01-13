// File: pkg/config/initializepop.go
package config

import (
	"drift/pkg/core"
	"drift/pkg/methods"
	"fmt"
)

// InitializePop creates and initializes a population based on the scenario parameter
func InitializePop(model *core.Model) *core.Pop {
	// Create base population structure
	pop := createEmptyPop(model)

	// Call the appropriate setup method
	var err error
	switch model.Scenario {
	case "default":
		err = methods.SetupPopDefault(model, pop)
	case "flood":
		err = methods.SetupPopFlood(model, pop)
	default:
		// If scenario not recognized, use default
		print("Setting up default pop")
		err = methods.SetupPopDefault(model, pop)
	}

	if err != nil {
		panic(fmt.Sprintf("Failed to initialize population: %v", err))
	}

	return pop
}

// createEmptyPop creates the base population structure
func createEmptyPop(model *core.Model) *core.Pop {
	pop := &core.Pop{
		IndData:      make(map[int][]int),
		Chromosomes:  make(map[int][][]uint64),
		Centromeres:  make(map[int][2][]uint64),
		IndMutations: make(map[int]map[int][]int),
		MutationPool: make(map[int]core.Mutation),
		MutationHist: make(map[int]int),
		Tracking:     make(map[string]int),
		AlleleFreqs:  make(map[int][]int16),
		MaleDB:       make(map[int]core.Ancestor),
		FemaleDB:     make(map[int]core.Ancestor),
	}

	// Reset run-specific parameters
	model.FreeParameters["seed"] = -1
	model.FreeParameters["last_pop_size"] = 0

	// Reset tracking counters
	pop.Tracking["births"] = 0
	pop.Tracking["deaths"] = 0
	pop.Tracking["marriages"] = 0
	pop.Tracking["random_deaths"] = 0
	pop.Tracking["cull_deaths"] = 0

	return pop
}
