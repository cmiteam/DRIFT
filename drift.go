package main

import (
	"drift/modules/animations"
	"drift/modules/birth"
	"drift/modules/coalescence"
	"drift/modules/death"
	"drift/modules/initializedrift"
	"drift/modules/initializemodel"
	"drift/modules/initializepop"
	"drift/modules/necalcs"
	"drift/modules/parsecommands"
	"drift/modules/save"
	"drift/modules/seedpopulation"
	"drift/modules/utils"
	"drift/pkg/simulation"
	"fmt"
	"log"
	"os"
	"time"
)

// Main function does the following:
// 1. Parses command-line arguments
// 2. Initializes the model
// 3. Runs the model for the specified number of iterations
// 4. Saves the results
// 5. Prints the execution time

func main() {
	// Start the timer for execution time
	starttime := time.Now()

	// Parse the command-line arguments
	config := parsecommands.ParseCommandLine()

	// Initialize the model
	model, err := initializemodel.InitializeModel(config.ConfigRoot)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error initializing model: %v\n", err)
		os.Exit(1)
	}

	// Setup directories
	err = initializedrift.SetupDirectories(config.Results)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		os.Exit(1)
	}

	// Initialize animations if enabled
	animContainer, err := animations.InitializeIfEnabled(model, config.MapRoot)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error initializing animations: %v\n", err)
		os.Exit(1)
	}

	// Setup CPU profiling
	cleanup, err := initializedrift.StartCPUProfiling(config.CpuProfile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		os.Exit(1)
	}
	defer cleanup()

	// Loop over the number of model runs
	for run := 1; run <= int(model.Parameters["num_runs"]); run++ {
		print("\nRun ", run, "\n")
		model.FreeParameters["run"] = run

		pop := initializepop.InitializePop(model)

		// Loop over the number years in each model run
		for year := 0; year <= int(model.Parameters["end_year"]); year++ {
			model.FreeParameters["year"] = year
			if year >= int(model.Parameters["seed_year"]) && model.FreeParameters["seed"] == -1 {
				seedpopulation.SeedThePopulation(model, pop)
				model.FreeParameters["seed"] = 1
				fmt.Println("   Seeded population in year", year, "with seed style", int(model.Parameters["seed_style"]))
			}

			birth.Birth(model, pop)
			simulation.Mating(model, pop)
			death.Death(model, pop)
			model.FreeParameters["last_pop_size"] = len(pop.IndData) // save pop size for future growth rate calculations

			// Escape clauses
			// Save and quit if population extinct
			if len(pop.IndData) <= 1 {
				save.Save(model, pop, animContainer)
				break
			}
			// Save and quit if genealo = pop size
			genealo := utils.CountGenealo(&pop.IndData)

			// Save and quit if genealo = 0
			if genealo == 0 {
				save.Save(model, pop, animContainer)
				break
			}
			if model.Parameters["track_map"] == 1 && year%int(model.Parameters["animation_save_interval"]) == 0 {
				err := animations.AddAnimationFrames(model, pop, animContainer)
				if err != nil {
					log.Printf("Error adding animation frames: %v", err)
				}
			}
			// Normal save at save interval
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
			//			print("Saving map...\n")
			err := animations.SaveAllGIFs(animContainer, config.Results)
			if err != nil {
				log.Printf("Error saving animation: %v", err)
			}
		}
		if model.Parameters["save_detailed_SFS"] == 1 {
			necalcs.SaveSFSTimeSeries(model, pop)
		}
		if model.Parameters["track_coalescence"] == 1 {
			yadamResult := coalescence.FindYAdam(model, pop)
			mteveResult := coalescence.FindMtEve(model, pop)
			fmt.Printf("\nY-Adam: ID %d, born year %d, %d generations back\n",
				yadamResult.YAdamID, yadamResult.YAdamBirthYear, yadamResult.GenerationsBack)
			fmt.Printf("Mt-Eve: ID %d, born year %d, %d generations back\n",
				mteveResult.MtEveID, mteveResult.MtEveBirthYear, mteveResult.GenerationsBack)
			fmt.Printf("Final MaleDB size: %d, FemaleDB size: %d\n", len(pop.MaleDB), len(pop.FemaleDB))
			fmt.Printf("Living individuals: %d\n", len(pop.IndData))
		}
	}

	// End of model runs
	elapsed := time.Since(starttime)
	fmt.Printf("Execution time: %s\n", elapsed)
	fmt.Print("\a")
}
