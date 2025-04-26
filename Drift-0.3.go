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
	"time"
)

// Default value for the config-root parameter, relative path to config files
const defaultConfigRoot = "static"

// Default value for the map-root parameter, relative path to map files
const defaultMapRoot = "maps"

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
	// Add more parameters as needed

	// Parse the command-line arguments
	flag.Parse()

	// Initialize the model.
	// If there is an error, print it to stderr and exit with a non-zero status code.
	model, err := initializemodel.InitializeModel(*configRootArg)

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error initializing model: %v\n", err)
		os.Exit(1)
	}

	animContainer := animations.Initialize(model, *mapRootArg)

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

		// Things to do at the end of a model run
		if model.Parameters["track_DNA"] == 1 {
			filename := fmt.Sprintf("results/%s genome map.png", model.ModelName)
			pixelSize := 4
			save.SaveGenomeMap(pop.Chromosomes, model.ChromosomeArms, filename, pixelSize, int(model.Parameters["NumBits"]))
		}
		if model.Parameters["track_map"] == 1 {
			print("Saving map...\n")
			err := animations.SaveAllGIFs(animContainer)
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
