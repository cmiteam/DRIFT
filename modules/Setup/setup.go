package setup

import (
	"drift/modules/actuarialloader"
	"drift/modules/chromosomeloader"
	"drift/modules/paramloader"
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
)

type FreeParams struct {
	Innoculated map[string]int
	IndID       int
	LastPopSize int
	MutID       int
	Seed        []int
}

func SetupFreeParams(modelparameters map[string]float64) FreeParams {
	freeParams := FreeParams{
		Innoculated: make(map[string]int),
		IndID:       int(modelparameters["start_pop_size"]) - 1,
		LastPopSize: int(modelparameters["start_pop_size"]),
		MutID:       0,
		Seed:        []int{},
	}
	return freeParams
}

func createCSV(filename string, headers []string) error {

	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()
	writer := csv.NewWriter(file)
	defer writer.Flush()
	err = writer.Write(headers)
	if err != nil {
		return err
	}
	return nil

}

func InitializeOutputFiles(modelParameters map[string]float64, run int, resultsDirectory string) error {
	modelID := "Model 1"
	runStr := strconv.Itoa(run)

	filename := filepath.Join(resultsDirectory, fmt.Sprintf("%s-%s_results.csv", modelID, runStr))
	headers := []string{"run", "year", "numinds", "marriages", "births", "random_deaths", "culled_deaths", "genetic", "genealo", "Y", "mt", "cents", "numblocks", "avblocksize", "sdblocksize", "AvFitPerInd", "AvBinFit", "NumMuts", "AvMutsPerInd", "AvMutsPerBin"}
	err := createCSV(filename, headers)
	if err != nil {
		return err
	}

	if modelParameters["track_dead"] == 1 {
		filename := filepath.Join(resultsDirectory, fmt.Sprintf("%s-%s_deaths.csv", modelID, runStr))
		headers := []string{"ID", "birthyear", "deathyear", "sex", "father", "mother", "lifespan", "lat", "lon", "married", "numbirths", "Ygens", "MTgens", "MinGenealGens", "MaxGenealGens", "SeedAlleles", "CentromereCount", "blocks", "fitness", "NumMuts", "CauseOfDeath"}
		err := createCSV(filename, headers)
		if err != nil {
			return err
		}
	}

	if modelParameters["mutation_hist"] == 1 {
		filename := filepath.Join(resultsDirectory, fmt.Sprintf("%s-%s_mutation_histogram.csv", modelID, runStr))
		headers := []string{"Effect", "All", "EndOfRun"}
		err := createCSV(filename, headers)
		if err != nil {
			return err
		}
	}

	return nil
}

func Initialize(modelParametersPath, chromosomeDataPath, actuarialTablePath string) (map[string]int, map[string]bool, []chromosomeloader.ChromosomeArm, map[int]float64, FreeParams, error) {
	modelParametersPath := "C:\\Go\\Programs\\Drift\\static\\parameter_defaults.csv"
	chromosomeDataPath := "C:\\Go\\Programs\\Drift\\static\\chromosome_data.csv"
	actuarialTablePath := "C:\\Go\\Programs\\Drift\\static\\actuarial_table.csv"

	modelParameters, plotParameters, err := paramloader.LoadParametersFromCSV(modelParametersPath)
	if err != nil {
		return nil, nil, nil, nil, FreeParams{}, fmt.Errorf("Error loading parameters: %v", err)
	}

	// Load chromosome arms data
	chromosomeArms, err := chromosomeloader.LoadChromosomeArmsFromCSV(chromosomeDataPath, int(modelParameters["multiplier"]))
	if err != nil {
		return nil, nil, nil, nil, FreeParams{}, fmt.Errorf("Error loading chromosome arms: %v", err)
	}

	// Load actuarial table
	deathRisk, err := actuarialloader.LoadActuarialTable(actuarialTablePath)
	if err != nil {
		return nil, nil, nil, nil, FreeParams{}, fmt.Errorf("Error loading actuarial table: %v", err)
	}

	// Setup free params
	freeParams := SetupFreeParams(modelParameters)

	// Initialize output files
	run := 1
	resultsDirectory := "./results"
	err = InitializeOutputFiles(modelParameters, run, resultsDirectory)
	if err != nil {
		return nil, nil, nil, nil, FreeParams{}, fmt.Errorf("Error loading actuarial table: %v", err)
	}

	return modelParameters, plotParameters, chromosomeArms, deathRisk, freeParams, nil
}
