// go run drift3.go

package main

import (
	"drift/modules/birth"
	"drift/modules/death"
	"drift/modules/initializemodel"
	"drift/modules/initializepop"
	"drift/modules/marriage"
	"drift/modules/save"
	"drift/modules/seedpopulation"
	"fmt"
	"time"
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

func main() {

	starttime := time.Now()
	model := initializemodel.InitializeModel()

	// Loop over the number of model runs
	for run := 1; run <= int(model.Parameters["num_runs"]); run++ {
		print("\nRun ", run, "\n")
		pop := initializepop.InitializePop(model)

		// Loop over the number years in each model run
		for year := 0; year <= int(model.Parameters["end_year"]); year++ {
			if year >= int(model.Parameters["seed_year"]) && model.FreeParameters["seed"] == -1 {
				seedpopulation.SeedThePopulation(model, pop, year)
			}

			birth.Birth(model, pop, year)
			marriage.Marriage(model, pop, year)
			death.Death(model, pop, year, run)
			model.FreeParameters["last_pop_size"] = len(pop.IndData) // save pop size for future growth rate calculations
			if year%int(model.Parameters["save_interval"]) == 0 {
				save.Save(model, pop, run, year)
			}
			if len(pop.IndData) <= 1 { // Save and quit if population extinct
				save.Save(model, pop, run, year)
				break
			}
		}

		// Things to do at the end of a model run
		if model.Parameters["track_DNA"] == 1 {
			filename := fmt.Sprintf("results/%s genome map.png", model.ModelName)
			pixelSize := 4
			save.SaveGenomeMap(pop.Chromosomes, model.ChromosomeArms, filename, pixelSize, int(model.Parameters["NumBits"]))
		}
	}
	elapsed := time.Since(starttime)
	fmt.Printf("Execution time: %s\n", elapsed)
	fmt.Print("\a")
}
