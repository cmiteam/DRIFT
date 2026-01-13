package main

import (
	"drift/pkg/analysis"
	"drift/pkg/config"
	"drift/pkg/events"
	"drift/pkg/simulation"
	"drift/pkg/utils"
	"drift/pkg/visualization"
	"drift/pkg/webserver"
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
	commands := config.ParseCommandLine()

	// If web mode is enabled, start the web server
	if commands.WebMode {
		fmt.Println("Starting Drift web server...")
		err := webserver.Start()
		if err != nil {
			fmt.Fprintf(os.Stderr, "Web server error: %v\n", err)
			os.Exit(1)
		}
		return
	}

	// Initialize the model
	model, err := config.InitializeModel(commands.ConfigRoot)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error initializing model: %v\n", err)
		os.Exit(1)
	}

	// Setup user-specific directories if username is specified
	err = config.SetupUserDirectories(model)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error setting up user directories: %v\n", err)
		os.Exit(1)
	}

	// If output-dir was specified via command line, use it (overrides user directory)
	if commands.OutputDir != "" {
		model.ResultsDir = commands.Results
	}

	// Setup directories
	err = config.SetupDirectories(model.ResultsDir)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		os.Exit(1)
	}

	// Initialize animations if enabled
	animContainer, err := visualization.InitializeIfEnabled(model, commands.MapRoot)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error initializing animations: %v\n", err)
		os.Exit(1)
	}

	// Setup CPU profiling
	cleanup, err := config.StartCPUProfiling(commands.CpuProfile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		os.Exit(1)
	}
	defer cleanup()

	// Loop over the number of model runs
	for run := 1; run <= int(model.Parameters["num_runs"]); run++ {
		print("\nRun ", run, "\n")
		model.FreeParameters["run"] = run

		pop := config.InitializePop(model)

		// Loop over the number years in each model run
		for year := 0; year <= int(model.Parameters["end_year"]); year++ {
			model.FreeParameters["year"] = year
			if year >= int(model.Parameters["seed_year"]) && model.FreeParameters["seed"] == -1 {
				events.Seed(model, pop)
				model.FreeParameters["seed"] = 1
				fmt.Println("   Seeded population in year", year, "with seed style", int(model.Parameters["seed_style"]))
			}

			simulation.Birth(model, pop)
			simulation.Mating(model, pop)
			simulation.Death(model, pop)
			model.FreeParameters["last_pop_size"] = len(pop.IndData) // save pop size for future growth rate calculations

			// Escape clauses
			// Save and quit if population extinct
			if len(pop.IndData) <= 1 {
				simulation.Save(model, pop, animContainer)
				break
			}
			// Save and quit if genealo = pop size
			genealo := utils.CountGenealo(&pop.IndData)

			// Save and quit if genealo = 0
			if genealo == 0 {
				simulation.Save(model, pop, animContainer)
				break
			}
			if model.Parameters["track_map"] == 1 && year%int(model.Parameters["animation_save_interval"]) == 0 {
				err := visualization.AddAnimationFrames(model, pop, animContainer)
				if err != nil {
					log.Printf("Error adding animation frames: %v", err)
				}
			}
			// Normal save at save interval
			if year%int(model.Parameters["save_interval"]) == 0 {
				simulation.Save(model, pop, animContainer)
			}
		}

		// Things to do at the end of a model run
		if model.Parameters["track_DNA"] == 1 {
			filename := fmt.Sprintf("%s/%s genome map.png", model.ResultsDir, model.ModelName)
			pixelSize := 4
			simulation.SaveGenomeMap(pop.Chromosomes, model.ChromosomeArms, filename, pixelSize, int(model.Parameters["NumBits"]))
		}
		if model.Parameters["track_map"] == 1 {
			//			print("Saving map...\n")
			err := visualization.SaveAllGIFs(animContainer, model.ResultsDir)
			if err != nil {
				log.Printf("Error saving animation: %v", err)
			}
		}
		if model.Parameters["save_detailed_SFS"] == 1 {
			analysis.SaveSFSTimeSeries(model, pop)
		}
		if model.Parameters["track_coalescence"] == 1 {
			yadamResult := analysis.FindYAdam(model, pop)
			mteveResult := analysis.FindMtEve(model, pop)
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
