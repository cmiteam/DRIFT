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
		Centromeres:  make(map[int][]uint64),
		IndMutations: make(map[int]map[int][]int),
		MutationPool: make(map[int]types.Mutation),
		MutationHist: make(map[int]int),
		Tracking:     make(map[string]int),
	}
	
	// Reset run-specific parameters
	model.FreeParameters["indID"] = 0 // Starting ID for individuals
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
	var landCoordinates [][2]float64
	for lat := range model.Map {
		for lon := range model.Map[lat] {
			// Check if this is land (terrain type 1)
			if model.Map[lat][lon] == 1 {
				landCoordinates = append(landCoordinates, [2]float64{lat, lon})
				//print("appending land coordinate: ", lat, ",", lon, "\n")
			}
		}
	}
	if len(landCoordinates) == 0 {
		landCoordinates = append(landCoordinates, [2]float64{0, 0})
	}

	// Set up the individuals
	popSize := int(model.Parameters["start_pop_size"])
	fitness := int(model.Parameters["mu_scale_factor"])

	for i := 0; i < popSize; i++ {
		// Choose a random land location
		randomIndex := rand.Intn(len(landCoordinates))
		randomLoc := landCoordinates[randomIndex]

		// Convert to integer coordinates with 10x multiplier for precision
		lat := int(randomLoc[0] * 10.0)
		lon := int(randomLoc[1] * 10.0)

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
		
		pop.IndData[i][individual.BirthYear] = -age // the person was born before the model began to be run
		pop.IndData[i][individual.Lifespan] = int(model.Parameters["lifespan"]) // initial theoretical lifespans
		pop.IndData[i][individual.Sex] = rand.Intn(2) // 0 = male, 1 = female
		pop.IndData[i][individual.MarriageState] = -1 // will be set to the ID # of the spouse
		pop.IndData[i][individual.NumBirths] = 0 // tracks number of children for females
		pop.IndData[i][individual.LastBirthYear] = 0 // to allow for spacing between children
		pop.IndData[i][individual.Fitness] = fitness // used for survival calculations
		pop.IndData[i][individual.AlleleCount] = 0 // tracking descent from seed individual(s)
		pop.IndData[i][individual.Lat] = lat
		pop.IndData[i][individual.Lon] = lon
		
		//print("Individual", i, " position: ", lat, ",", lon, "\n")
		model.FreeParameters["indID"]++ // each ind gets a unique ID
	}

	model.FreeParameters["last_pop_size"] = len(pop.IndData) // needed to control population growth

	return pop

}
