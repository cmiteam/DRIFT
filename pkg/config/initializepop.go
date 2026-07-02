// File: pkg/config/initializepop.go
package config

import (
	"drift/pkg/core"
	"drift/pkg/methods"
	"drift/pkg/modules"
	"fmt"
)

// Register the population-setup modules into the central registry (pkg/modules).
// These are keyed by scenario name; setup dispatch is scenario-driven by default
// but can be overridden with a setup_style parameter.
func init() {
	modules.RegisterSetup("default", methods.SetupPopDefault)
	modules.RegisterSetup("flood", methods.SetupPopFlood)
	modules.RegisterSetup("creation", methods.SetupPopCreation)
}

// InitializePop creates and initializes a population using the selected setup
// module (scenario-driven, or the setup_style parameter if set).
func InitializePop(model *core.Model) *core.Pop {
	pop := createEmptyPop(model)
	if err := modules.DispatchSetup(model, pop); err != nil {
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
