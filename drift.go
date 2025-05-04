// go run drift3.go

package main

import (
	"drift/modules/animations"
	"drift/modules/birth"
	"drift/modules/death"
	"drift/modules/initializemodel"
	"drift/modules/initializepop"
	"drift/modules/marriage"
	"drift/modules/save"
	"drift/modules/seedpopulation"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	"time"
)

// Default value for the config-root parameter, relative path to config files
const defaultConfigRoot = "static"

// Default value for the map-root parameter, relative path to map files
const defaultMapRoot = "maps"

// Default value for the results parameter, relative path to results files
const defaultResults = "results"

// Main function does the following:
// 1. Parses command-line arguments
// 2. Initializes the model
// 3. Runs the model for the specified number of iterations
// 4. Saves the results
// 5. Prints the execution time
func main() {
	// Start the timer for execution time
	starttime := time.Now()

	// Define a command-line parameter (e.g., for a config file path)
	configRootArg := flag.String("config-root",
		defaultConfigRoot,
		"path to directory containing configuration files")
	mapRootArg := flag.String("map-root",
		defaultMapRoot,
		"path to directory containing map files")
	resultsArg := flag.String("results",
		defaultResults,
		"path to directory containing results files")

	// Add more parameters as needed

	// Parse the command-line arguments
	flag.Parse()

	// Initialize the model
	model, err := initializemodel.InitializeModel(*configRootArg)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error initializing model: %v\n", err)
		os.Exit(1)
	}

	// Create the results directory if it doesn't exist
	if _, err := os.Stat(*resultsArg); os.IsNotExist(err) {
		err = os.Mkdir(*resultsArg, 0755)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating results directory: %v\n", err)
		os.Exit(1)
	}

	// Initialize the animations
	animContainer, err := animations.Initialize(model, *mapRootArg)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error initializing animations: %v\n", err)
		os.Exit(1)
	}

	// Start CPU profiling
	cpuFile, err := os.Create("cpu.prof")
	if err != nil {
		panic(err)
	}
	defer cpuFile.Close()
	pprof.StartCPUProfile(cpuFile)

	// Loop over the number of model runs
	for run := 1; run <= int(model.Parameters["num_runs"]); run++ {
		print("\nRun ", run, "\n")
		model.FreeParameters["run"] = run
		pop := initializepop.InitializePop(model)

		// Loop over the number years in each model run
		for year := 0; year <= int(model.Parameters["end_year"]); year++ {
			model.FreeParameters["year"] = year
			if year >= int(model.Parameters["seed_year"]) &&
				model.FreeParameters["seed"] == -1 {
				seedpopulation.SeedThePopulation(model, pop)
			}

			birth.Birth(model, pop)
			marriage.Marriage(model, pop)
			death.Death(model, pop)
			model.FreeParameters["last_pop_size"] = len(pop.IndData) // save pop size for future growth rate calculations
			if len(pop.IndData) <= 1 {                               // Save and quit if population extinct
				save.Save(model, pop, animContainer)
				break
			}
			if year%int(model.Parameters["save_interval"]) == 0 {
				save.Save(model, pop, animContainer)
			}
		}

		pprof.StopCPUProfile()

		// Things to do at the end of a model run
		if model.Parameters["track_DNA"] == 1 {
			filename := fmt.Sprintf("results/%s genome map.png", model.ModelName)
			pixelSize := 4
			save.SaveGenomeMap(pop.Chromosomes, model.ChromosomeArms, filename, pixelSize, int(model.Parameters["NumBits"]))
		}
		if model.Parameters["track_map"] == 1 {
			print("Saving map...\n")
			err := animations.SaveAllGIFs(animContainer, *resultsArg)
			if err != nil {
				log.Printf("Error saving animation: %v", err)
			}
		}
	}

	// End of model runs
	elapsed := time.Since(starttime)
	fmt.Printf("Execution time: %s\n", elapsed)
	fmt.Print("\a")
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
