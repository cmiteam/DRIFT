package marriage

import (
	"drift/types"
	"math/rand"
	"time"
)

func Marriage(model *types.Model, pop *types.Pop, year int) {
	var availableMen, availableWomen []int

	// Find eligible individuals
	for id, data := range pop.IndData {
		if data["marriage_state"] == -1 && year-data["birth_year"] >= int(model.Parameters["maturity"]) {
			if data["sex"] == 0 {
				availableMen = append(availableMen, id)
			} else if data["sex"] == 1 {
				availableWomen = append(availableWomen, id)
			}
		}
	}

	// Randomize people (note this will create unusual age gaps among married couples, fix?)
	rand.Seed(time.Now().UnixNano())
	rand.Shuffle(len(availableMen), func(i, j int) {
		availableMen[i], availableMen[j] = availableMen[j], availableMen[i]
	})
	rand.Shuffle(len(availableWomen), func(i, j int) {
		availableWomen[i], availableWomen[j] = availableWomen[j], availableWomen[i]
	})

	// Trim male and female lists to the shortest of the two
	if len(availableMen) > len(availableWomen) {
		availableMen = availableMen[:len(availableWomen)]
	} else {
		availableWomen = availableWomen[:len(availableMen)]
	}

	// Assign spouses
	for i := 0; i < len(availableMen); i++ {
		pop.IndData[availableMen[i]]["marriage_state"] = availableWomen[i]
		pop.IndData[availableWomen[i]]["marriage_state"] = availableMen[i]
		pop.Tracking["marriages"]++
		// Uncomment for debugging
		// fmt.Printf(" M: %d : %d\n", availableMen[i], availableWomen[i])
	}
}
