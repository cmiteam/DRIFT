package events

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"fmt"
	"math"
	"sort"
)

// migrationPulse is a one-time gene-flow burst (roadmap §2): in year `year` it
// relabels a `fraction` of deme `from`'s living members into deme `to`. Unlike the
// steady island-model migration (mating.go's per-year deme_migration_rate) or a
// demographic scenario's continuous migration matrix (§6c), this is a single
// discrete pulse — e.g. a founder wave colonizing a new region, or a rescue
// influx after a local crash.
//
// RNG is consumed ONLY in the firing year, and only over the source deme's members
// in sorted-id order (a partial Fisher-Yates), so a fixed rng_seed stays
// reproducible and every other year matches a no-event baseline exactly. It works
// whether demes come from the plain island model (num_demes > 1) or a demographic
// scenario — both track subpopulation membership in the same individual.Deme field.
//
// Spec: migration_pulse year=<year> from=<deme> to=<deme> fraction=<0..1>
//   - year     (required) the single year the pulse fires
//   - from,to  (required) source and destination deme ids (must differ)
//   - fraction (required, in (0,1]) share of the source deme that moves
type migrationPulse struct {
	year, from, to int
	fraction       float64
}

func newMigrationPulse(fields map[string]string) (Event, error) {
	year, err := reqInt(fields, "year")
	if err != nil {
		return nil, fmt.Errorf("migration_pulse: %w", err)
	}
	from, err := reqInt(fields, "from")
	if err != nil {
		return nil, fmt.Errorf("migration_pulse: %w", err)
	}
	to, err := reqInt(fields, "to")
	if err != nil {
		return nil, fmt.Errorf("migration_pulse: %w", err)
	}
	if from == to {
		return nil, fmt.Errorf("migration_pulse: from and to are the same deme (%d)", from)
	}
	fraction, err := reqFloat(fields, "fraction")
	if err != nil {
		return nil, fmt.Errorf("migration_pulse: %w", err)
	}
	if fraction <= 0 || fraction > 1 {
		return nil, fmt.Errorf("migration_pulse: fraction must be in (0,1], got %g", fraction)
	}
	return &migrationPulse{year: year, from: from, to: to, fraction: fraction}, nil
}

// Apply relabels round(fraction * |from|) random members of the source deme into
// the destination deme, but only in the pulse's year. Members are gathered and
// sorted before any RNG is drawn, so the draw sequence is a deterministic function
// of the seed.
func (mp *migrationPulse) Apply(model *core.Model, pop *core.Pop, year int) {
	if year != mp.year {
		return
	}
	members := make([]int, 0)
	for id := range pop.IndData {
		if pop.IndData[id][individual.Deme] == mp.from {
			members = append(members, id)
		}
	}
	sort.Ints(members)

	n := int(math.Round(mp.fraction * float64(len(members))))
	if n <= 0 {
		return
	}
	if n > len(members) {
		n = len(members)
	}
	// Partial Fisher-Yates over the sorted member list: pick n distinct members
	// without replacement and relabel them into the destination deme.
	for i := 0; i < n; i++ {
		j := i + utils.RandIntn(len(members)-i)
		members[i], members[j] = members[j], members[i]
		pop.IndData[members[i]][individual.Deme] = mp.to
	}
	fmt.Printf("migration_pulse: year %d — moved %d individuals from deme %d to deme %d\n", year, n, mp.from, mp.to)
}

func init() { RegisterEvent("migration_pulse", newMigrationPulse) }
