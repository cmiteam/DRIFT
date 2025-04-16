package initializemodel

import (
	"drift/modules/actuarialloader"
	"drift/modules/chromosomeloader"
	"drift/modules/paramloader"
	"drift/modules/save"
	"drift/types"
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"path/filepath"
)

type configLoader struct {
	Fn       func(*types.Model, [][]string) error
	FileName string
}

func loadConfig(configLoader configLoader,
	model *types.Model,
	configRoot string) error {
	configPath := filepath.Join(configRoot, configLoader.FileName)
	file, err := os.Open(configPath)
	if err != nil {
		return err
	}
	defer file.Close()
	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		return err
	}

	// Call the function with the model and records.
	// The function will modify the model in place.
	// We don't bother to check for a nonzero number of records
	// because it's conceivable that some config files might be
	// intentionally empty.
	return configLoader.Fn(model, records)
}

func InitializeModel(configRoot string) (*types.Model, error) {
	model := &types.Model{
		Parameters:     make(map[string]float64),
		PlotFlags:      make(map[string]bool),
		ChromosomeArms: make(map[int]map[int][]int),
		DeathRisk:      make(map[int]float64),
		CumulativeProb: make(map[int]float64),
		FreeParameters: make(map[string]int),
	}

	// Load files
	loaders := []configLoader{
		{paramloader.LoadParameters, "parameter_defaults.csv"},
		{chromosomeloader.LoadChromosomeArmsFromCSV, "chromosome_data.csv"},
		{actuarialloader.LoadActuarialTable, "actuarial_table.csv"},
	}

	for _, loader := range loaders {
		err := loadConfig(loader, model, configRoot)
		if err != nil {
			return nil,
				fmt.Errorf("Error loading %s: %v", loader.FileName, err.Error)
		}
	}

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

	return model, nil
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
