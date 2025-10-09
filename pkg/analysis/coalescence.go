package analysis

import (
	"drift/modules/individual"
	"drift/pkg/core"
	"fmt"
)

// FindYAdam identifies Y-chromosome Adam by tracing paternal lineages
func FindYAdam(model *core.Model, pop *core.Pop) *core.YAdamResult {
	// Get all living males
	var livingMales []int
	for id, indData := range pop.IndData {
		if indData[individual.Sex] == 0 { // Male
			livingMales = append(livingMales, id)
		}
	}

	if len(livingMales) == 0 {
		return &core.YAdamResult{
			YAdamID:         -1,
			GenerationsBack: -1,
			YAdamBirthYear:  -1,
			FoundYAdam:      false,
		}
	}

	// Count how many living males trace back to each ancestor
	ancestorCounts := make(map[int]int)

	for _, maleID := range livingMales {
		ancestors := tracePaternalLineage(pop, maleID)
		for _, ancestorID := range ancestors {
			ancestorCounts[ancestorID]++
		}
	}

	// Find ancestors that ALL living males share
	numLivingMales := len(livingMales)
	var commonAncestors []int

	for ancestorID, count := range ancestorCounts {
		if count == numLivingMales {
			commonAncestors = append(commonAncestors, ancestorID)
		}
	}

	// Debug output
	fmt.Printf("Y-Adam analysis: %d living males, %d unique ancestors, %d common ancestors\n",
		numLivingMales, len(ancestorCounts), len(commonAncestors))

	if len(commonAncestors) == 0 {
		// Show how close we are - find the ancestor with highest count
		maxCount := 0
		var bestAncestor int
		for ancestorID, count := range ancestorCounts {
			if count > maxCount {
				maxCount = count
				bestAncestor = ancestorID
			}
		}
		fmt.Printf("Closest to Y-Adam: ancestor %d shared by %d/%d males (%.1f%%)\n",
			bestAncestor, maxCount, numLivingMales, float64(maxCount)/float64(numLivingMales)*100.0)

		return &core.YAdamResult{
			YAdamID:         -1,
			GenerationsBack: -1,
			YAdamBirthYear:  -1,
			FoundYAdam:      false,
		}
	}

	// Find the most recent common ancestor (latest birth year)
	yadamID := commonAncestors[0]
	latestBirthYear := getAncestorBirthYear(pop, yadamID)

	for _, candidateID := range commonAncestors[1:] {
		candidateBirthYear := getAncestorBirthYear(pop, candidateID)
		if candidateBirthYear > latestBirthYear {
			yadamID = candidateID
			latestBirthYear = candidateBirthYear
		}
	}

	generationsBack := model.FreeParameters["year"] - latestBirthYear

	return &core.YAdamResult{
		YAdamID:         yadamID,
		GenerationsBack: generationsBack,
		YAdamBirthYear:  latestBirthYear,
		FoundYAdam:      true,
	}
}

// FindMtEve identifies mitochondrial Eve by tracing maternal lineages
func FindMtEve(model *core.Model, pop *core.Pop) *core.MtEveResult {
	// Get all living individuals (both sexes for maternal lineage)
	var livingIndividuals []int
	for id := range pop.IndData {
		livingIndividuals = append(livingIndividuals, id)
	}

	if len(livingIndividuals) == 0 {
		return &core.MtEveResult{
			MtEveID:         -1,
			GenerationsBack: -1,
			MtEveBirthYear:  -1,
			FoundMtEve:      false,
		}
	}

	// Count how many living individuals trace back to each maternal ancestor
	ancestorCounts := make(map[int]int)

	for _, individualID := range livingIndividuals {
		ancestors := traceMaternalLineage(pop, individualID)
		for _, ancestorID := range ancestors {
			ancestorCounts[ancestorID]++
		}
	}

	// Find ancestors that ALL living individuals share
	numLivingIndividuals := len(livingIndividuals)
	var commonAncestors []int

	for ancestorID, count := range ancestorCounts {
		if count == numLivingIndividuals {
			commonAncestors = append(commonAncestors, ancestorID)
		}
	}

	// Debug output
	fmt.Printf("Mt-Eve analysis: %d living individuals, %d unique ancestors, %d common ancestors\n",
		numLivingIndividuals, len(ancestorCounts), len(commonAncestors))

	if len(commonAncestors) == 0 {
		// Show how close we are
		maxCount := 0
		var bestAncestor int
		for ancestorID, count := range ancestorCounts {
			if count > maxCount {
				maxCount = count
				bestAncestor = ancestorID
			}
		}
		fmt.Printf("Closest to Mt-Eve: ancestor %d shared by %d/%d individuals (%.1f%%)\n",
			bestAncestor, maxCount, numLivingIndividuals, float64(maxCount)/float64(numLivingIndividuals)*100.0)

		return &core.MtEveResult{
			MtEveID:         -1,
			GenerationsBack: -1,
			MtEveBirthYear:  -1,
			FoundMtEve:      false,
		}
	}

	// Find the most recent common ancestor (latest birth year)
	mteveID := commonAncestors[0]
	latestBirthYear := getAncestorBirthYear(pop, mteveID)

	for _, candidateID := range commonAncestors[1:] {
		candidateBirthYear := getAncestorBirthYear(pop, candidateID)
		if candidateBirthYear > latestBirthYear {
			mteveID = candidateID
			latestBirthYear = candidateBirthYear
		}
	}

	generationsBack := model.FreeParameters["year"] - latestBirthYear

	return &core.MtEveResult{
		MtEveID:         mteveID,
		GenerationsBack: generationsBack,
		MtEveBirthYear:  latestBirthYear,
		FoundMtEve:      true,
	}
}

// tracePaternalLineage traces back through fathers using MaleDB
func tracePaternalLineage(pop *core.Pop, startID int) []int {
	var lineage []int
	currentID := startID
	visited := make(map[int]bool)

	for currentID != -1 && !visited[currentID] {
		lineage = append(lineage, currentID)
		visited[currentID] = true

		// Look up father in MaleDB
		if ancestor, exists := pop.MaleDB[currentID]; exists {
			currentID = ancestor.ID
		} else {
			break // No father recorded
		}
	}

	return lineage
}

// traceMaternalLineage traces back through mothers using FemaleDB
func traceMaternalLineage(pop *core.Pop, startID int) []int {
	var lineage []int
	currentID := startID
	visited := make(map[int]bool)

	for currentID != -1 && !visited[currentID] {
		lineage = append(lineage, currentID)
		visited[currentID] = true

		// Look up mother in FemaleDB
		if ancestor, exists := pop.FemaleDB[currentID]; exists {
			currentID = ancestor.ID
		} else {
			break // No mother recorded
		}
	}

	return lineage
}

// getAncestorBirthYear gets birth year for an ancestor using the genealogical databases
func getAncestorBirthYear(pop *core.Pop, ancestorID int) int {
	// First try to find in living population
	if indData, exists := pop.IndData[ancestorID]; exists {
		return indData[individual.BirthYear]
	}

	// If not living, search through MaleDB and FemaleDB to find birth year
	for _, ancestor := range pop.MaleDB {
		if ancestor.ID == ancestorID {
			return ancestor.BirthYear
		}
	}

	for _, ancestor := range pop.FemaleDB {
		if ancestor.ID == ancestorID {
			return ancestor.BirthYear
		}
	}

	// If not found anywhere, return a very old birth year
	return -99999
}
