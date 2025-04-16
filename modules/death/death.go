package death

import (
	"drift/types"
	"fmt"
	"math/rand"
	"os"
	"strings"
	"time"
)

func Death(model *types.Model, pop *types.Pop, year int, run int) int {

	rand.Seed(time.Now().UnixNano())
	deaths := 0
	var deadPeopleData string
	keyList := generateKeyList(pop.IndData)

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

		age := year - pop.IndData[ind]["birth_year"]
		ageGroup := int((float64(age)/float64(pop.IndData[ind]["lifespan"]))*model.Parameters["min_lifespan"]/5) * 5
		if ageGroup > 85 {
			ageGroup = 85
		}
		deathrisk := model.DeathRisk[ageGroup]
		die := rand.Float64() // low roll = death
		riskModification := model.Parameters["min_lifespan"] / float64(pop.IndData[ind]["lifespan"])
		fitness := 1.0
		if model.Parameters["track_mutations"] == 1 {
			fitness = float64(pop.IndData[ind]["fitness"]) / model.Parameters["mu_scale_factor"]
		}
		adjustedDeathRisk := deathrisk * riskModification * fitness
		if die < adjustedDeathRisk {
			if int(model.Parameters["track_dead"]) == 1 {
				deadPersonString := personDataString(ind, pop, year, "R")
				deadPeopleData += deadPersonString
			}
			RIP(ind, pop, model)
			deaths++
			model.Parameters["random_deaths"]++
		}
	}

	// Adjust max population size based on bottleneck
	maxPopSize := int(model.Parameters["max_pop_size"])
	if int(model.Parameters["bottleneck_start"]) <= year && int(model.Parameters["bottleneck_end"]) >= year {
		maxPopSize = int(model.Parameters["bottleneck_size"])
	}

	// Step 2: Trim excess population by randomly culling individuals
	excess := len(pop.IndData) - maxPopSize
	for excess > 0 {
		keyList = generateKeyList(pop.IndData)
		randomIndex := rand.Intn(len(keyList))
		ind := keyList[randomIndex]
		if ind == model.FreeParameters["seed"] { // Don't kill off the seed
			continue
		}
		if int(model.Parameters["track_dead"]) == 1 {
			deadPersonString := deadString(ind, pop, year)
			deadPersonString += ",R\n"
			deadPeopleData += deadPersonString
		}
		RIP(ind, pop, model)
		deaths++
		pop.Tracking["cull_deaths"]++
		excess = len(pop.IndData) - maxPopSize
	}

	// Step 3: Tamp down population growth rate by randomly culling individuals
	allowedNumInds := int(float64(model.FreeParameters["last_pop_size"]) * model.Parameters["max_growth_rate"])
	if allowedNumInds > int(model.Parameters["max_pop_size"]) {
		allowedNumInds = int(model.Parameters["max_pop_size"])
	}

	diff := len(pop.IndData) - allowedNumInds
	for diff > 0 {
		keyList = generateKeyList(pop.IndData)
		randomIndex := rand.Intn(len(keyList))
		ind := keyList[randomIndex]
		if ind == model.FreeParameters["seed"] { // Don't kill off the seed
			continue
		}
		if int(model.Parameters["track_dead"]) == 1 {
			deadPersonString := deadString(ind, pop, year)
			deadPersonString += ",R\n"
			deadPeopleData += deadPersonString
		}
		RIP(ind, pop, model)
		deaths++
		pop.Tracking["cull_deaths"]++
		diff = len(pop.IndData) - allowedNumInds
	}

	// Step 4: Reduce population to specified number of breeding individuals, if called for, by randomly culling individuals
	if model.Parameters["max_breeding_inds"] > -1 {
		keyList = generateKeyList(pop.IndData)
		breeders := countBreedingIndividuals(pop, year, model)
		for breeders > int(model.Parameters["max_breeding_inds"]) {
			randomIndex := rand.Intn(len(keyList))
			ind := keyList[randomIndex]
			if ind == model.FreeParameters["seed"] { // Don't kill off the seed
				continue
			}
			if int(model.Parameters["track_dead"]) == 1 {
				deadPersonString := deadString(ind, pop, year)
				deadPersonString += ",R\n"
				deadPeopleData += deadPersonString
			}
			RIP(ind, pop, model)
			deaths++
			pop.Tracking["cull_deaths"]++
			breeders = countBreedingIndividuals(pop, year, model)
		}
	}

	// Save dead individuals to file
	if int(model.Parameters["track_dead"]) == 1 {
		modelID := fmt.Sprintf("%.0f", model.Parameters["model_id"])
		filename := fmt.Sprintf("results_directory/%s-%d deaths.csv", modelID, run)
		writeToFile(filename, deadPeopleData)
	}

	return deaths
}

// generateKeyList creates a slice of all individual IDs
func generateKeyList(indData map[int]map[string]int) []int {
	keyList := make([]int, 0, len(indData))
	for key := range indData {
		keyList = append(keyList, key)
	}
	return keyList
}

// RIP removes a deceased individual and updates related data
func RIP(ind int, pop *types.Pop, model *types.Model) {
	if _, exists := pop.IndData[ind]; exists {
		if pop.IndData[ind]["marriage_state"] > -1 {
			spouse := pop.IndData[ind]["marriage_state"]
			pop.IndData[spouse]["marriage_state"] = -1
		}
	}
	delete(pop.Chromosomes, ind)
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
	delete(pop.IndMutations, ind)
	delete(pop.IndData, ind)
}

// deadString formats individual data for death records
func deadString(ind int, pop *types.Pop, year int) string {
	// keeps track of deceased individuals if they are to be saved
	var info strings.Builder
	indInfo, exists := pop.IndData[ind]
	if exists {
		getValue := func(key string, defaultVal int) int {
			if val, ok := indInfo[key]; ok {
				return val
			}
			return defaultVal
		}
		info.WriteString(fmt.Sprintf("%d,%d,%d,%d,%d,%d,%d,",
			ind,
			getValue("birth_year", -1),
			year,
			getValue("sex", -1),
			getValue("dad", -1),
			getValue("mom", -1),
			getValue("lifespan", -1),
		))
		info.WriteString(fmt.Sprintf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
			getValue("lat", -1),
			getValue("lon", -1),
			getValue("marriage_state", -1),
			getValue("numbirths", -1),
			getValue("Y_gens", -1),
			getValue("mt_gens", -1),
			getValue("min_genealo_gens", -1),
			getValue("max_genealo_gens", -1),
			getValue("allele_count", -1),
			getValue("num_blocks", -1),
			getValue("centromeres", -1),
			getValue("fitness", -1),
			getValue("mutations", -1),
		))
	}
	return info.String()
}

// personDataString formats detailed individual data with state information
func personDataString(ind int, pop *types.Pop, year int, state string) string {
	var info strings.Builder
	indInfo, exists := pop.IndData[ind]
	if exists {
		fields := []string{
			"birth_year", "sex", "dad", "mom", "lifespan",
			"lat", "lon", "marriage_state", "numbirths",
			"Y_gens", "mt_gens", "min_genealo_gens", "max_genealo_gens",
			"allele_count", "num_blocks", "centromeres", "fitness", "mutations",
		}
		info.WriteString(fmt.Sprintf("%d,%d,%d,", ind, getOrDefault(indInfo, "birth_year", -1), year))
		for _, field := range fields {
			info.WriteString(fmt.Sprintf("%d,", getOrDefault(indInfo, field, -1)))
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
func countBreedingIndividuals(pop *types.Pop, year int, model *types.Model) int {
	count := 0
	for _, data := range pop.IndData {
		age := year - data["birth_year"]
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
