package initializemodel

import (
	"drift/modules/actuarialloader"
	"drift/modules/chromosomeloader"
	"drift/modules/paramloader"
	"drift/modules/save"
	"drift/types"
	"math"
)

// Initializes the model based on the configuration files.
func InitializeModel(configRoot string) (*types.Model, error) {
	model := &types.Model{
		Parameters:     make(map[string]float64),
		PlotFlags:      make(map[string]bool),
		ChromosomeArms: make(map[int]map[int][]int),
		DeathRisk:      make(map[int]float64),
		CumulativeProb: make(map[int]float64),
		FreeParameters: make(map[string]int),
		Map:            nil, // this will be set later, once the dimensions are known
	}

	// Attempt to load each config file. Failure will be fatal.
	err := paramloader.LoadParameters(model, configRoot)
	if err != nil {
		return nil, err
	}
	err = chromosomeloader.LoadChromosomeArms(model, configRoot)
	if err != nil {
		return nil, err
	}
	err = actuarialloader.LoadActuarialTable(model, configRoot)
	if err != nil {
		return nil, err
	}

	// Calculate derived values
	model.Parameters["mu_sig_figs"] = math.Pow(1, model.Parameters["mu_sig_figs"])

	// Prepare output files
	save.SaveHeaders(model.ModelName)

	// Initialize free parameters
	model.FreeParameters["indID"] = 0         // Starting ID for individuals
	model.FreeParameters["seed"] = -1         // No seed initially
	model.FreeParameters["last_pop_size"] = 0 // Required for growth rate calculations
	model.FreeParameters["mutID"] = 0         // Starting ID for mutations

	return model, nil
}
