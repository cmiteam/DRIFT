package simulation

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"fmt"
	"math"
	"os"
	"sort"
	"strings"
)

// deathStandard is the default death algorithm, registered as "standard" in
// module_registry.go. Callers use the Death dispatcher.
func deathStandard(model *core.Model, pop *core.Pop) int {

	deaths := 0
	var deadPeopleData string
	keyList := generateKeyList(&pop.IndData, -1)

	// Step 1: Random actuarial deaths
	for _, ind := range keyList {
		// People that can potentially live for a long time (e.g., 900 years) need
		// lower death risks or they will NEVER reach that age.
		// Since you cannot test people who can *potentially* live for a long time
		// at the same rate as normal people, we will create a variable called
		// ageGroup. This is an 'effective' age group, proportional to potential
		// lifespan. For example, if a person could potentially live to 850 and the
		// normal lifespan is 85, at age 85 they are only at 1/10 of their potential
		// lifespan. So, at age 85: age group = 85/850 x 85 = 8.5, which rounds down
		// to 5 because the risk factor data are in 5-year increments.
		// The death risk is really high at 85 already, so anyone who makes it to
		// >= 85 has the same risk of dying each year.

		age := model.FreeParameters["year"] - pop.IndData[ind][individual.BirthYear]
		ageGroup := int((float64(age)/float64(pop.IndData[ind][individual.Lifespan]))*model.Parameters["min_lifespan"]/5) * 5
		if ageGroup > 85 {
			ageGroup = 85
		}
		deathrisk := model.DeathRisk[ageGroup]
		die := utils.RandFloat64() // low roll = death
		riskModification := model.Parameters["min_lifespan"] / float64(pop.IndData[ind][individual.Lifespan])
		fitness := 1.0
		if model.Parameters["track_mutations"] == 1 {
			fitness = float64(pop.IndData[ind][individual.Fitness]) / model.Parameters["mu_scale_factor"]
		}
		adjustedDeathRisk := deathrisk * riskModification * fitness
		if die < adjustedDeathRisk {
			if model.FreeParameters["seed"] == ind {
				if age < pop.IndData[ind][individual.Lifespan] {
					break
				}
			}
			if int(model.Parameters["track_dead"]) == 1 {
				deadPersonString := personDataString(model, pop, ind, "R")
				deadPeopleData += deadPersonString
			}
			RIP(ind, pop, model)
			deaths++
			pop.Tracking["random_deaths"]++
		}
	}

	// Steps 2 & 3: enforce carrying capacity. When a demographic scenario is active
	// (roadmap §6c), model.DemeCaps holds a per-deme census target for this year and
	// supersedes the global max_pop_size / max_growth_rate limits — each deme is
	// trimmed to its own cap. Otherwise the original single-population behavior runs
	// unchanged (strict no-op for non-demography runs: byte-identical results).
	if len(model.DemeCaps) > 0 {
		deaths += cullByDeme(model, pop, &deadPeopleData)
	} else {
		// Adjust max population size based on bottleneck
		maxPopSize := int(model.Parameters["max_pop_size"])
		if int(model.Parameters["bottleneck_start"]) <= model.FreeParameters["year"] && int(model.Parameters["bottleneck_end"]) >= model.FreeParameters["year"] {
			maxPopSize = int(model.Parameters["bottleneck_size"])
		}

		// Step 2: Trim excess population by randomly culling individuals
		excess := len(pop.IndData) - maxPopSize
		keyList = generateKeyList(&pop.IndData, model.FreeParameters["seed"])
		keyList = pickVictims(excess, keyList)
		for _, ind := range keyList {
			if model.FreeParameters["seed"] == ind {
				break
			}
			if int(model.Parameters["track_dead"]) == 1 {
				deadPersonString := deadString(model, pop, ind)
				deadPersonString += ",R\n"
				deadPeopleData += deadPersonString
			}
			RIP(ind, pop, model)
			deaths++
			pop.Tracking["cull_deaths"]++
		}

		// Step 3: Tamp down population growth rate by randomly culling individuals.
		// Round up so small populations can grow at all: int(6 * 1.05) = 6 (locks pop),
		// but ceil(6 * 1.05) = 7 (lets it grow by one).
		allowedNumInds := int(math.Ceil(float64(model.FreeParameters["last_pop_size"]) * model.Parameters["max_growth_rate"]))
		if allowedNumInds > int(model.Parameters["max_pop_size"]) {
			allowedNumInds = int(model.Parameters["max_pop_size"])
		}

		diff := len(pop.IndData) - allowedNumInds
		keyList = generateKeyList(&pop.IndData, model.FreeParameters["seed"])
		keyList = pickVictims(diff, keyList)
		for _, ind := range keyList {
			if model.FreeParameters["seed"] == ind {
				break
			}
			if int(model.Parameters["track_dead"]) == 1 {
				deadPersonString := deadString(model, pop, ind)
				deadPersonString += ",R\n"
				deadPeopleData += deadPersonString
			}
			RIP(ind, pop, model)
			deaths++
			pop.Tracking["cull_deaths"]++
		}
	}

	// Step 4: Reduce population to specified number of breeding individuals, if called for, by randomly culling individuals
	if model.Parameters["max_breeding_inds"] > -1 {
		if breeders := countBreedingIndividuals(model, pop); breeders > int(model.Parameters["max_breeding_inds"]) {
			keyList := generateKeyList(&pop.IndData, model.FreeParameters["seed"])
			keyList = pickVictims(breeders-int(model.Parameters["max_breeding_inds"]), keyList)
			for _, ind := range keyList {
				if model.FreeParameters["seed"] == ind {
					break
				}
				if int(model.Parameters["track_dead"]) == 1 {
					deadPersonString := deadString(model, pop, ind)
					deadPersonString += ",R\n"
					deadPeopleData += deadPersonString
				}
				RIP(ind, pop, model)
				deaths++
				pop.Tracking["cull_deaths"]++
				breeders = countBreedingIndividuals(model, pop)
			}
		}
	}

	// Save dead individuals to file
	if int(model.Parameters["track_dead"]) == 1 {
		modelID := fmt.Sprintf("%.0f", model.Parameters["model_id"])
		filename := fmt.Sprintf("results_directory/%s-%d deaths.csv", modelID, model.FreeParameters["run"])
		writeToFile(filename, deadPeopleData)
	}

	return deaths
}

// cullByDeme trims each deme to its per-deme census cap (model.DemeCaps), used when
// a demographic scenario is active (roadmap §6c) instead of the global carrying-
// capacity limits. Demes are processed in ascending id order and individuals culled
// in sorted order (via generateKeyList + pickVictims) so a fixed rng_seed stays
// reproducible. Returns the number of individuals removed.
func cullByDeme(model *core.Model, pop *core.Pop, deadPeopleData *string) int {
	seed := model.FreeParameters["seed"]
	trackDead := int(model.Parameters["track_dead"]) == 1

	// Group living individuals by deme (skipping the protected seed individual).
	byDeme := make(map[int][]int)
	for id := range pop.IndData {
		if id == seed {
			continue
		}
		d := pop.IndData[id][individual.Deme]
		byDeme[d] = append(byDeme[d], id)
	}

	// Iterate demes in ascending id order for deterministic RNG consumption.
	demes := make([]int, 0, len(byDeme))
	for d := range byDeme {
		demes = append(demes, d)
	}
	sort.Ints(demes)

	removed := 0
	for _, d := range demes {
		cap, ok := model.DemeCaps[d]
		if !ok {
			continue // no cap declared for this deme: leave it untouched
		}
		members := byDeme[d]
		sort.Ints(members)
		excess := len(members) - cap
		victims := pickVictims(excess, members)
		for _, ind := range victims {
			if trackDead {
				*deadPeopleData += deadString(model, pop, ind) + ",R\n"
			}
			RIP(ind, pop, model)
			removed++
			pop.Tracking["cull_deaths"]++
		}
	}
	return removed
}

// generateKeyList creates a slice of all individual IDs
func generateKeyList(indData *map[int][]int, seed int) []int {
	keyList := make([]int, 0, len(*indData))
	for key := range *indData {
		if key == seed {
			continue
		}
		keyList = append(keyList, key)
	}
	// Deterministic order before RNG consumption (actuarial rolls, culling), so a
	// fixed rng_seed reproduces the run despite randomized map iteration order.
	sort.Ints(keyList)
	return keyList
}

func pickVictims(excess int, keyList []int) []int {
	if excess <= 0 {
		return []int{}
	}
	victimList := make([]int, excess)
	for i := 0; i < excess; i++ {
		randomIndex := i + utils.RandIntn(len(keyList)-i)
		victimList[i] = keyList[randomIndex]
		keyList[randomIndex] = keyList[i]
	}
	return victimList
}

// RIP removes a deceased individual and updates related data
func RIP(ind int, pop *core.Pop, model *core.Model) {
	if _, exists := pop.IndData[ind]; exists {
		if pop.IndData[ind][individual.MarriageState] > -1 {
			// Clear the dying person's spouse's marriage state
			spouse := pop.IndData[ind][individual.MarriageState]
			if _, spouseExists := pop.IndData[spouse]; spouseExists {
				pop.IndData[spouse][individual.MarriageState] = -1
			} else {
				fmt.Printf("Cannot remove spouse %d of ind %d, already dead\n", spouse, ind)
			}
		}
	} else {
		fmt.Printf("Cannot kill %d, already dead\n", ind)
	}

	// decrement mutation counts
	for strand := 0; strand <= 1; strand++ {
		if mutationIDs, exists := pop.IndMutations[ind][strand]; exists {
			for _, mutationID := range mutationIDs {
				if mutation, exists := pop.MutationPool[mutationID]; exists {
					mutation.Count--
					if mutation.Count <= 0 {
						delete(pop.MutationPool, mutationID)
					} else {
						pop.MutationPool[mutationID] = mutation
					}
				}
			}
		}
	}

	if pop.IndData[ind][individual.Sons] == 0 {
		delete(pop.MaleDB, ind)
	}
	if pop.IndData[ind][individual.Daughters] == 0 {
		delete(pop.FemaleDB, ind)
	}
	delete(pop.IndMutations, ind)
	delete(pop.Chromosomes, ind)
	delete(pop.Centromeres, ind)
	delete(pop.IndData, ind)

}

// deadString formats individual data for death records
func deadString(model *core.Model, pop *core.Pop, ind int) string {
	// keeps track of deceased individuals if they are to be saved
	var info strings.Builder
	indInfo, exists := pop.IndData[ind]
	if exists {
		// Append individual data to the string
		info.WriteString(fmt.Sprintf("%d,%d,%d,%d,%d,%d,%d,",
			ind,
			indInfo[individual.BirthYear],
			model.FreeParameters["year"],
			indInfo[individual.Sex],
			indInfo[individual.Dad],
			indInfo[individual.Mom],
			indInfo[individual.Lifespan],
		))
		info.WriteString(fmt.Sprintf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
			indInfo[individual.Lat],
			indInfo[individual.Lon],
			indInfo[individual.MarriageState],
			indInfo[individual.NumBirths],
			indInfo[individual.YGens],
			indInfo[individual.MtGens],
			indInfo[individual.MinGenealoGens],
			indInfo[individual.MaxGenealoGens],
			indInfo[individual.AlleleCount],
			indInfo[individual.NumBlocks],
			indInfo[individual.NumCentromeres],
			indInfo[individual.Fitness],
			indInfo[individual.NumMutations],
		))
	}
	return info.String()
}

// personDataString formats detailed individual data with state information
func personDataString(model *core.Model, pop *core.Pop, ind int, state string) string {
	var info strings.Builder
	indInfo, exists := pop.IndData[ind]
	if exists {
		fields := []individual.IndData{
			individual.BirthYear, individual.Sex, individual.Dad, individual.Mom, individual.Lifespan,
			individual.Lat, individual.Lon, individual.MarriageState, individual.NumBirths,
			individual.YGens, individual.MtGens, individual.MinGenealoGens, individual.MaxGenealoGens,
			individual.AlleleCount, individual.NumBlocks, individual.NumCentromeres, individual.Fitness,
			individual.NumMutations,
		}
		info.WriteString(fmt.Sprintf("%d,%d,%d,", ind, indInfo[individual.BirthYear], model.FreeParameters["year"]))
		for field := range fields {
			info.WriteString(fmt.Sprintf("%d,", indInfo[field]))
		}
		info.WriteString(state) // Append state ('R' for removed, 'A' for alive)
		info.WriteString("\n")
	}
	return info.String()
}

// getOrDefault retrieves a value from a map with a default fallback
func getOrDefault(data map[string]int, key string, defaultVal int) int {
	if val, ok := data[key]; ok {
		return val
	}
	return defaultVal
}

// countBreedingIndividuals counts individuals of breeding age
func countBreedingIndividuals(model *core.Model, pop *core.Pop) int {
	count := 0
	for _, data := range pop.IndData {
		age := model.FreeParameters["year"] - data[individual.BirthYear]
		if age >= int(model.Parameters["maturity"]) {
			count++
		}
	}
	return count
}

// writeToFile writes content to a file
func writeToFile(filename, content string) error {
	file, err := os.OpenFile(filename, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		return fmt.Errorf("could not open file %s: %w", filename, err)
	}
	defer file.Close()

	_, err = file.WriteString(content)
	if err != nil {
		return fmt.Errorf("could not write to file %s: %w", filename, err)
	}
	return nil
}
