package initializemodel

import (
	"drift/modules/actuarialloader"
	"drift/modules/chromosomeloader"
	"drift/modules/paramloader"
	"drift/modules/save"
	"drift/types"
	"fmt"
	"math"
)

func InitializeModel() *types.Model {
	model := &types.Model{
		Parameters:     make(map[string]float64),
		PlotFlags:      make(map[string]bool),
		ChromosomeArms: make(map[int]map[int][]int),
		DeathRisk:      make(map[int]float64),
		CumulativeProb: make(map[int]float64),
		FreeParameters: make(map[string]int),
	}

	// Load files
	paramloader.LoadParameters(model)
	chromosomeloader.LoadChromosomeArmsFromCSV(model)
	actuarialloader.LoadActuarialTable(model)

	// Calculate derived values
	model.ModelName = getModelName(model.Parameters)
	model.Parameters["mu_sig_figs"] = math.Pow(1, model.Parameters["mu_sig_figs"])

	// Prepare output files
	save.SaveHeaders(model.ModelName)

	// Initialize free parameters
	model.FreeParameters["indID"] = 0         // Starting ID for individuals
	model.FreeParameters["seed"] = -1         // No seed initially
	model.FreeParameters["last_pop_size"] = 0 // Required for growth rate calculations
	model.FreeParameters["mutID"] = 0         // Starting ID for mutations

	//PrintModel(model)                       // For doublechecking purposes

	return model

}

func getModelName(params map[string]float64) string {
	if id, ok := params["model_id"]; ok {
		return fmt.Sprintf("model%.0f", id)
	}
	return "default_model"
}

func PrintModel(model *types.Model) {
	fmt.Println("Model Name:", model.ModelName)
	fmt.Println("Parameters:", model.Parameters)
	fmt.Println("Free Parameters:", model.FreeParameters)
	fmt.Println("Plot Flags:", model.PlotFlags)
	fmt.Println("Chromosome Arms:", model.ChromosomeArms)
	fmt.Println("Death Risk:", model.DeathRisk)
	fmt.Println("Cumulative Probability:", model.CumulativeProb)
}
