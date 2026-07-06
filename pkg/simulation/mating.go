package simulation

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/modules"
	"drift/pkg/utils"
	"sort"
)

func Mating(model *core.Model, pop *core.Pop) {
	var availableMen, availableWomen []int

	// Find eligible individuals (unmarried and of mature age)
	for id, data := range pop.IndData {
		if data[individual.MarriageState] == -1 &&
			model.FreeParameters["year"]-data[individual.BirthYear] >= int(model.Parameters["maturity"]) {
			if data[individual.Sex] == 0 {
				availableMen = append(availableMen, id)
			} else if data[individual.Sex] == 1 {
				availableWomen = append(availableWomen, id)
			}
		}
	}

	// Sort so the lists have a deterministic order before the mating module
	// consumes randomness (Go map iteration order is randomized; without this a
	// fixed rng_seed would still produce different runs).
	sort.Ints(availableMen)
	sort.Ints(availableWomen)

	// When a demographic scenario is active (roadmap §6c), the scheduler already
	// migrated individuals per the epoch's migration matrix (in Apply, before Birth)
	// and owns population structure, so here we only mate within each deme that is
	// currently present in the population — no scalar island migration.
	if model.DemographyScheduler != nil {
		mateByDeme(model, pop, availableMen, availableWomen)
		return
	}

	numDemes := int(model.Parameters["num_demes"])
	if numDemes <= 1 {
		// Single-population default: identical to pre-deme behavior.
		modules.DispatchMating(model, pop, availableMen, availableWomen)
		return
	}

	// Island model (roadmap §6b): individuals may migrate between demes, then mate
	// only within their (possibly new) deme. This gives the discrete population
	// structure the differentiation statistics (Fst / f-stats) measure.
	migrateDemes(model, pop, numDemes)

	// Mate within each deme independently by dispatching the selected module once
	// per deme. Reuses every existing mating module unchanged (they simply operate
	// on the deme's subset of eligible men/women). Iterate demes in ascending order
	// so RNG is consumed deterministically.
	menByDeme := bucketByDeme(pop, availableMen)
	womenByDeme := bucketByDeme(pop, availableWomen)
	for deme := 0; deme < numDemes; deme++ {
		men := menByDeme[deme]
		women := womenByDeme[deme]
		if len(men) == 0 || len(women) == 0 {
			continue
		}
		modules.DispatchMating(model, pop, men, women)
	}
}

// mateByDeme dispatches the selected mating module once per deme actually present
// in the population (unlike the island path, the set of demes changes over a
// demographic run as splits activate new demes). Demes are iterated in ascending
// order so RNG is consumed deterministically.
func mateByDeme(model *core.Model, pop *core.Pop, availableMen, availableWomen []int) {
	menByDeme := bucketByDeme(pop, availableMen)
	womenByDeme := bucketByDeme(pop, availableWomen)

	demeSet := make(map[int]struct{}, len(menByDeme)+len(womenByDeme))
	for d := range menByDeme {
		demeSet[d] = struct{}{}
	}
	for d := range womenByDeme {
		demeSet[d] = struct{}{}
	}
	demes := make([]int, 0, len(demeSet))
	for d := range demeSet {
		demes = append(demes, d)
	}
	sort.Ints(demes)

	for _, deme := range demes {
		men := menByDeme[deme]
		women := womenByDeme[deme]
		if len(men) == 0 || len(women) == 0 {
			continue
		}
		modules.DispatchMating(model, pop, men, women)
	}
}

// migrateDemes moves each individual to a uniformly-random *other* deme with
// probability deme_migration_rate. This is the island model's gene-flow channel:
// a migrant subsequently mates in its new deme. Individuals are visited in sorted
// order so the RNG stream (hence the whole run) stays reproducible under a fixed
// rng_seed. A no-op when the migration rate is 0.
func migrateDemes(model *core.Model, pop *core.Pop, numDemes int) {
	m := model.Parameters["deme_migration_rate"]
	if m <= 0 || numDemes < 2 {
		return
	}
	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	sort.Ints(ids)
	for _, id := range ids {
		if utils.RandFloat64() >= m {
			continue
		}
		cur := pop.IndData[id][individual.Deme]
		// Draw a different deme uniformly: pick in [0, numDemes-1) then shift past
		// the current one.
		dst := utils.RandIntn(numDemes - 1)
		if dst >= cur {
			dst++
		}
		pop.IndData[id][individual.Deme] = dst
	}
}

// bucketByDeme groups ids by their Deme field, preserving the (already sorted)
// input order within each bucket.
func bucketByDeme(pop *core.Pop, ids []int) map[int][]int {
	buckets := make(map[int][]int)
	for _, id := range ids {
		d := pop.IndData[id][individual.Deme]
		buckets[d] = append(buckets[d], id)
	}
	return buckets
}

// createInfluenceGrid creates a grid where each cell contains IDs of women who can mate there
func createInfluenceGrid(model *core.Model, pop *core.Pop, availableWomen []int, maxDistance int) [][][]int {
	// Get map dimensions
	mapHeight := 101
	mapWidth := 101
	if model.Parameters["track_map"] == 1 {
		mapHeight = len(model.Map)
		mapWidth = len(model.Map[0])
	}
	// Create influence grid - each cell contains list of women who can mate there
	influenceGrid := make([][][]int, mapHeight)
	for i := range influenceGrid {
		influenceGrid[i] = make([][]int, mapWidth)
		for j := range influenceGrid[i] {
			influenceGrid[i][j] = make([]int, 0)
		}
	}

	// For each woman, add her ID to all cells within mating range
	for _, womanID := range availableWomen {
		womanLat := pop.IndData[womanID][individual.Lat]
		womanLon := pop.IndData[womanID][individual.Lon]

		// Add woman's ID to all cells within Manhattan distance
		for lat := womanLat - maxDistance; lat <= womanLat+maxDistance; lat++ {
			for lon := womanLon - maxDistance; lon <= womanLon+maxDistance; lon++ {
				// Check bounds
				if lat >= 0 && lat < mapHeight && lon >= 0 && lon < mapWidth {
					influenceGrid[lat][lon] = append(influenceGrid[lat][lon], womanID)
				}
			}
		}
	}

	return influenceGrid
}

func removeWomanFromGrid(influenceGrid [][][]int, womanID int, womanLat, womanLon, maxDistance int) {
	// Remove woman's ID from all cells within Manhattan distance
	for lat := womanLat - maxDistance; lat <= womanLat+maxDistance; lat++ {
		for lon := womanLon - maxDistance; lon <= womanLon+maxDistance; lon++ {
			// Check bounds
			if lat >= 0 && lat < len(influenceGrid) && lon >= 0 && lon < len(influenceGrid[0]) {
				influenceGrid[lat][lon] = removeFromSlice(influenceGrid[lat][lon], womanID)
			}
		}
	}
}

func removeFromSlice(slice []int, value int) []int {
	for i, v := range slice {
		if v == value {
			// Remove by swapping with last element and truncating
			slice[i] = slice[len(slice)-1]
			return slice[:len(slice)-1]
		}
	}
	return slice // Not found, return unchanged
}
