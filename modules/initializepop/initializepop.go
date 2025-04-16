package initializepop

import (
	"drift/types"
	"fmt"
	"math/rand"
	"time"
)

func InitializePop(model *types.Model) *types.Pop {

	// Create a new population
	pop := &types.Pop{
		IndData:      make(map[int]map[string]int),
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

		pop.IndData[i] = map[string]int{
			"dad":              -1,                                // -1 is used often in this program as a placeholder
			"mom":              -1,                                // ditto
			"birth_year":       -age,                              // the person was born before the model began to be run
			"lifespan":         int(model.Parameters["lifespan"]), // initial theoretical lifespans
			"sex":              rand.Intn(2),                      // 0 = male, 1 = female
			"marriage_state":   -1,                                // will be set to the ID # of the spouse
			"num_births":       0,                                 // tracks number of children for females
			"last_birth_year":  0,                                 // to allow for spacing between children
			"fitness":          fitness,                           // used for survival calculations
			"allele_count":     0,                                 // tracking descent from seed individual(s)
			"Y_gens":           -1,                                // generations from male seed
			"mt_gens":          -1,                                // generations from female seed
			"min_genealo_gens": -1,                                // shortest path on family tree to seed
			"max_genealo_gens": -1,                                // longest path on family tree to seed
			"lat":              rand.Intn(1000) - 500,             // for non-random mating or geography
			"lon":              rand.Intn(1000) - 500,             // lat and lon are in a square centered on (0,0)
		}

		model.FreeParameters["indID"]++ // each ind gets a unique ID
	}

	model.FreeParameters["last_pop_size"] = len(pop.IndData) // needed to control population growth

	//PrintPop(pop)  // For doublechecking purposes
	return pop

}

func PrintPop(pop *types.Pop) {
	fmt.Println("Individual Data:", pop.IndData)
	fmt.Println("Chromosomes:", pop.Chromosomes)
	fmt.Println("Centromeres:", pop.Centromeres)
	fmt.Println("Individual Mutations:", pop.IndMutations)
	fmt.Println("Mutation Pool Size:", len(pop.MutationPool))
	fmt.Println("Mutation History Size:", len(pop.MutationHist))
	fmt.Println("Total Mutation Count:", pop.MutationCount)
	fmt.Println("Tracking Information:", pop.Tracking)
}
