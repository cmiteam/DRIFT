// File: pkg/config/initializepop.go
package config

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/methods"
	"drift/pkg/modules"
	"fmt"
	"sort"
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
	seedPedigreeFounders(model, pop)
	return pop
}

// seedPedigreeFounders gives every member of the initial population a pedigree node
// (TMR4A.md W1). Done here, once, after setup dispatch, so it covers every setup module
// (default / flood / creation / any future Eden template) without each having to know
// about the pedigree. Founder nodes carry no parents, which is the substrate for the
// walker's `censored_at_founding` outcome: the backward walk stops at creation and says
// so rather than forcing a root (TMR4A.md §0).
//
// Sorted iteration because Go randomizes map order and this must be reproducible; it
// draws no RNG, so the run stays byte-identical either way, but the resulting map is
// built in a fixed order.
func seedPedigreeFounders(model *core.Model, pop *core.Pop) {
	if model.Parameters["track_pedigree"] != 1 {
		return
	}
	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	sort.Ints(ids)
	for _, id := range ids {
		pop.PedEnsure(id,
			pop.IndData[id][individual.BirthYear],
			pop.IndData[id][individual.Sex])
	}
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
