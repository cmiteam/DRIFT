package simulation

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"math"
	"sort"
)

// Habitat suitability levers (roadmap §3).
//
// A habitat_suitability table (parsed in pkg/core) gives each map cell a suitability
// multiplier. Two death-side levers read it, both gated on model.HabitatActive() so a
// run with no habitat param takes the byte-identical original death path:
//
//  1. Mortality multiplier (death Step 1): a cell's actuarial death risk is scaled by
//     1/suitability, so a harsh cell (suitability 0.25) quadruples the local hazard.
//     This is exactly the shape of the famine event factor — same RNG draw per
//     individual, same sorted-id order, only the death threshold moves.
//
//  2. Soft spatial carrying capacity (death Step 3): the growth-ceiling cull picks its
//     victims weighted toward crowded and/or low-suitability cells instead of uniformly
//     at random. The global logistic / max_pop_size machinery still decides HOW MANY
//     die this year (so the total-N trajectory and every global-cap / num_demes /
//     DemeCaps interaction is unchanged); habitat only decides WHO, driving the
//     equilibrium per-cell occupancy toward a multiple of suitability. This is the
//     "carrying-capacity multiplier" the roadmap names, expressed as a redistribution
//     of survivors rather than a hard per-cell ceiling.
//
// Composition: the mortality lever lives in Step 1, which always runs, so it stacks
// with viability selection, famines, and the bottleneck. The weighted cull lives in
// Step 3, which is skipped entirely when a demographic scenario owns culling
// (model.DemeCaps set -> cullByDeme), so habitat is a no-op under a §6c scenario —
// consistent with how the global logistic carrying_capacity already yields there.

// habitatSuitabilityFloor guards the weight denominator against a zero/near-zero
// suitability (an occupant should never be on a suitability-0 cell, but be safe).
const habitatSuitabilityFloor = 1e-9

// habitatMortalityFactor maps a cell suitability to a multiplier on the actuarial
// death risk: fully suitable (1.0) leaves the risk unchanged, a harsh cell raises it
// (1/suitability), and a super-optimal refugium (suitability > 1) lowers it. A
// non-positive suitability returns 1.0 (no change) as a guard.
func habitatMortalityFactor(suitability float64) float64 {
	if suitability <= 0 {
		return 1.0
	}
	return 1.0 / suitability
}

// buildCellCounts tallies how many living individuals occupy each (Lat, Lon) cell —
// the local crowding used to weight the habitat cull. Ranges over the population map
// only to count (order-independent, no RNG).
func buildCellCounts(pop *core.Pop) map[[2]int]int {
	counts := make(map[[2]int]int, len(pop.IndData))
	for _, data := range pop.IndData {
		counts[[2]int{data[individual.Lat], data[individual.Lon]}]++
	}
	return counts
}

// pickVictimsHabitat selects `count` cull victims from the (sorted) candidate ids,
// weighted toward crowded / low-suitability cells — a weighted sample without
// replacement via the Efraimidis-Spirakis reservoir key (key_i = u_i^(1/w_i); the
// largest keys are selected). Each candidate's weight is crowding/suitability, so an
// individual in a packed, harsh cell is the likeliest victim. Exactly one RandFloat64
// is drawn per candidate, in the candidates' (already sorted) order, so a fixed
// rng_seed reproduces the run and two habitat runs are byte-identical.
//
// Only reached when model.HabitatActive(); the uniform pickVictims path is untouched
// for runs with no habitat table.
func pickVictimsHabitat(model *core.Model, pop *core.Pop, count int, candidates []int) []int {
	if count <= 0 {
		return nil
	}
	if count >= len(candidates) {
		// Everyone is culled. Still consume one draw per candidate so RNG
		// consumption does not depend on where the growth ceiling happens to land.
		for range candidates {
			utils.RandFloat64()
		}
		out := make([]int, len(candidates))
		copy(out, candidates)
		return out
	}

	cellCounts := buildCellCounts(pop)

	type scored struct {
		id    int
		score float64
	}
	items := make([]scored, len(candidates))
	for i, id := range candidates {
		lat := pop.IndData[id][individual.Lat]
		lon := pop.IndData[id][individual.Lon]
		s := model.CellSuitability(lat, lon)
		if s < habitatSuitabilityFloor {
			s = habitatSuitabilityFloor
		}
		crowd := float64(cellCounts[[2]int{lat, lon}])
		if crowd < 1 {
			crowd = 1
		}
		w := crowd / s // higher = likelier victim (packed and/or harsh)

		u := utils.RandFloat64()
		score := math.Inf(-1)
		if u > 0 {
			// ln(u)/w is monotonic in the ES key u^(1/w): a larger weight pushes the
			// (negative) score toward 0, so the crowded/harsh candidates rank highest.
			score = math.Log(u) / w
		}
		items[i] = scored{id: id, score: score}
	}

	sort.Slice(items, func(a, b int) bool {
		if items[a].score != items[b].score {
			return items[a].score > items[b].score
		}
		return items[a].id < items[b].id // deterministic tie-break
	})

	victims := make([]int, count)
	for i := 0; i < count; i++ {
		victims[i] = items[i].id
	}
	return victims
}
