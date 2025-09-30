package initializepop

import (
	"drift/modules/individual"
	"drift/types"
	"math/rand"
	"time"
)

func InitializePop(model *types.Model) *types.Pop {

	// Create a new population
	pop := &types.Pop{
		IndData:      make(map[int][]int),
		Chromosomes:  make(map[int][][]uint64),
		Centromeres:  make(map[int][2][]uint64),
		IndMutations: make(map[int]map[int][]int),
		MutationPool: make(map[int]types.Mutation),
		MutationHist: make(map[int]int),
		Tracking:     make(map[string]int),
		AlleleFreqs:  make(map[int][]int16),
		MaleDB:       make(map[int]types.Ancestor),
		FemaleDB:     make(map[int]types.Ancestor),
	}

	// Reset run-specific parameters
	model.FreeParameters["seed"] = -1 // No seed initially
	model.FreeParameters["last_pop_size"] = 0

	// Reset tracking counters for this run
	pop.Tracking["births"] = 0
	pop.Tracking["deaths"] = 0
	pop.Tracking["marriages"] = 0
	pop.Tracking["random_deaths"] = 0
	pop.Tracking["cull_deaths"] = 0

	// Initialize the random number generator
	rand.Seed(time.Now().UnixNano())

	// Find all land squares in the map
	var landCoordinates [][2]int
	if model.Parameters["track_map"] == 1 {
		for lat := range model.Map {
			for lon := range model.Map[lat] {
				// Check if this is land (terrain type 1)
				if model.Map[lat][lon] == 1 {
					landCoordinates = append(landCoordinates, [2]int{lat, lon})
				}
			}
		}
		if len(landCoordinates) == 0 {
			landCoordinates = append(landCoordinates, [2]int{0, 0})
		}
	}

	// Set up the individuals
	popSize := int(model.Parameters["start_pop_size"])
	fitness := int(model.Parameters["mu_scale_factor"])

	for i := 0; i < popSize; i++ {
		// assign data to each individual
		age := 0
		r := rand.Float64()
		for a, cumProb := range model.CumulativeProb {
			if r <= cumProb {
				age = a
				break
			}
		}

		pop.IndData[i] = individual.MakeIndData()
		pop.IndData[i][individual.BirthYear] = -age                             // the person was born before the model began to be run
		pop.IndData[i][individual.Lifespan] = int(model.Parameters["lifespan"]) // initial theoretical lifespans
		pop.IndData[i][individual.Sex] = rand.Intn(2)                           // 0 = male, 1 = female
		pop.IndData[i][individual.MarriageState] = -1                           // will be set to the ID # of the spouse
		pop.IndData[i][individual.NumBirths] = 0                                // tracks number of children for females
		pop.IndData[i][individual.LastBirthYear] = 0                            // to allow for spacing between children
		pop.IndData[i][individual.Fitness] = fitness                            // used for survival calculations
		pop.IndData[i][individual.AlleleCount] = 0                              // tracking descent from seed individual(s)
		pop.IndData[i][individual.NumCentromeres] = 0
		pop.IndData[i][individual.YGens] = -1
		pop.IndData[i][individual.MtGens] = -1
		pop.IndData[i][individual.MinGenealoGens] = -1
		pop.IndData[i][individual.MaxGenealoGens] = -1
		// Choose a random land location
		lat := rand.Intn(101)
		lon := rand.Intn(101)

		if model.Parameters["track_map"] == 1 {
			randomIndex := rand.Intn(len(landCoordinates))
			randomLoc := landCoordinates[randomIndex]
			lat = randomLoc[0]
			lon = randomLoc[1]
		}

		pop.IndData[i][individual.Lat] = lat
		pop.IndData[i][individual.Lon] = lon
		pop.IndData[i][individual.Sons] = 0
		pop.IndData[i][individual.Daughters] = 0

		// Create a paired individual and marry them if old enough
		//		pop.IndData[i+1] = pop.IndData[i]
		//		pop.IndData[i+1][individual.Sex] = 1
		//		model.FreeParameters["indID"] = i + 1 // tracks highest ID used so far

		//		if age >= int(model.Parameters["maturity"]) {
		//			pop.IndData[i][individual.MarriageState] = i + 1
		//			pop.IndData[i+1][individual.MarriageState] = i
		//			pop.IndData[i+1][individual.LastBirthYear] = -rand.Intn(int(model.Parameters["spacing"]))
		//		} else {
		//			pop.IndData[i][individual.MarriageState] = -1
		//			pop.IndData[i+1][individual.MarriageState] = -1
		//		}
		//		model.FreeParameters["indID"] = i + 1
	}

	model.FreeParameters["last_pop_size"] = len(pop.IndData) // needed to control population growth
	model.FreeParameters["indID"] = popSize                  // Starting ID for individuals

	return pop
}
