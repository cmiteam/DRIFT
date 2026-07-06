// Package demography implements time-varying demographic scenarios (roadmap §6c):
// per-deme population sizes, split events, growth, and a migration matrix that
// change over epochs. It exists so DRIFT can *run* the canonical out-of-Africa
// models (Gutenkunst et al. 2009 / Gravel et al. 2011) — and, in the same engine
// with the same statistics, alternative (e.g. Babel) histories.
//
// SPECIFICATION FORMAT. A scenario is a JSON file (e.g. a model's demography.json)
// authored in the PUBLISHED units of the demographic-inference literature so it
// reads as "exactly the model we ran": diploid effective sizes (Ne), times in
// GENERATIONS, and per-generation per-lineage migration rates. The file is an
// ordered list of epochs; each epoch declares the demes active during it (with an
// optional exponential growth rate), an optional split that begins it, and the
// pairwise migration matrix in force. See static/basemodels/OoA/demography.json.
//
// FORWARD-TIME MAPPING. DRIFT is forward-time and individual-based, with a loop
// that ticks in calendar YEARS (individuals mature at `maturity`, reproduce every
// year), whereas the published models are backward-time and specified in
// generations of an idealized Wright-Fisher population. Two header knobs bridge
// the two:
//
//   - generation_time (DRIFT years per generation, default 25): a published
//     generation count becomes generation_time DRIFT years.
//   - rescale Q (default 1): the standard simulation rescaling — divide every Ne
//     and every generation count by Q and multiply every rate by Q. This shrinks
//     an ~9000-generation, ~14000-diploid model into something that runs in
//     minutes while keeping drift (∝ generations/Ne), 2Nm, and growth-per-epoch
//     invariant. Q is the single speed⇄fidelity knob; Q=1 runs the model at full
//     published scale (long).
//
// So an epoch of D published generations spans round(D/Q * generation_time) DRIFT
// years; a published Ne becomes a census cap of round(Ne/Q); a per-generation
// migration rate m becomes a per-YEAR per-individual probability m*Q/generation_time.
//
// SPLITS (chosen with Rob: relabel-a-subset). A forward-time split of an existing
// deme carves a child deme out of it by RELABELLING a random subset of the parent's
// living individuals — preserving shared ancestry, so Fst / f-statistics / joint
// SFS reflect the real population tree (seeding fresh founders would falsely look
// like independent origins). The child's founding size is round(child_founders/Q)
// individuals (capped at half the parent so the parent is never emptied).
//
// ENGINE INTEGRATION. Apply is called once per year from the run loop. It (1) fires
// any split due this year, (2) migrates individuals per the epoch matrix, and (3)
// writes the per-deme census caps into model.DemeCaps. Death then culls each deme
// to its cap (pkg/simulation/death.go) and Mating mates within each present deme
// and skips its scalar island-migration (pkg/simulation/mating.go), because the
// scheduler owns migration while active. All randomness is consumed in sorted-id
// order so a fixed rng_seed stays byte-reproducible.
package demography

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"sort"
)

// sizeSpec is one deme's size (and optional growth) during an epoch, in published
// units (diploid Ne; growth is the per-generation exponential rate).
type sizeSpec struct {
	ID     int     `json:"id"`
	Size   float64 `json:"size"`
	Growth float64 `json:"growth,omitempty"`
}

// migSpec is one symmetric pairwise migration rate during an epoch, per published
// generation per lineage.
type migSpec struct {
	A    int     `json:"a"`
	B    int     `json:"b"`
	Rate float64 `json:"rate"`
}

// splitSpec begins an epoch by carving deme Child out of deme Parent, seeding it
// with ChildFounders (published diploid Ne) relabelled individuals.
type splitSpec struct {
	Parent        int     `json:"parent"`
	Child         int     `json:"child"`
	ChildFounders float64 `json:"child_founders"`
}

// Epoch is one interval of constant demographic structure.
type Epoch struct {
	Name        string     `json:"name"`
	DurationGen float64    `json:"duration_gen"`
	Split       *splitSpec `json:"split,omitempty"`
	Demes       []sizeSpec `json:"demes"`
	Migration   []migSpec  `json:"migration,omitempty"`

	startYear int // DRIFT year this epoch begins (inclusive), computed in finalize
	endYear   int // DRIFT year this epoch ends (exclusive), computed in finalize
}

// Scenario is a parsed demographic model. It implements core.DemographyScheduler.
type Scenario struct {
	Name           string  `json:"name"`
	Source         string  `json:"source"`
	GenerationTime float64 `json:"generation_time"`
	Rescale        float64 `json:"rescale"`
	Epochs         []Epoch `json:"epochs"`

	totalYears int // DRIFT year at which the last epoch ends
}

// Load reads and validates a scenario JSON file, computing the DRIFT-year timeline.
func Load(path string) (*Scenario, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("demography: %w", err)
	}
	var s Scenario
	if err := json.Unmarshal(data, &s); err != nil {
		return nil, fmt.Errorf("demography: parsing %s: %w", path, err)
	}
	if err := s.finalize(); err != nil {
		return nil, fmt.Errorf("demography: %s: %w", path, err)
	}
	return &s, nil
}

// finalize fills in defaults and computes each epoch's DRIFT-year span.
func (s *Scenario) finalize() error {
	if len(s.Epochs) == 0 {
		return fmt.Errorf("no epochs")
	}
	if s.GenerationTime <= 0 {
		s.GenerationTime = 25
	}
	if s.Rescale <= 0 {
		s.Rescale = 1
	}
	year := 0
	for i := range s.Epochs {
		e := &s.Epochs[i]
		if e.DurationGen <= 0 {
			return fmt.Errorf("epoch %q has non-positive duration_gen", e.Name)
		}
		if len(e.Demes) == 0 {
			return fmt.Errorf("epoch %q lists no demes", e.Name)
		}
		e.startYear = year
		// Rescale divides generation counts by Q; generation_time converts to years.
		span := int(math.Round(e.DurationGen / s.Rescale * s.GenerationTime))
		if span < 1 {
			span = 1
		}
		year += span
		e.endYear = year
	}
	s.totalYears = year
	return nil
}

// TotalYears is the DRIFT-year length of the scenario (end of the last epoch).
// Callers can use it to set end_year so the run covers the whole history.
func (s *Scenario) TotalYears() int { return s.totalYears }

// DemeRanks returns each deme's serial-founder rank: the number of founder
// bottleneck (split) steps between it and an origin deme. Demes present in the
// first epoch are origins (rank 0); a deme carved off by a split inherits its
// parent's rank + 1. This is DRIFT's intrinsic "distance from origin" axis for
// the serial-founder-effect vs single-origin heterozygosity-gradient comparison
// (roadmap §6c). In the published out-of-Africa analyses (Ramachandran 2005,
// Prugnolle 2005) heterozygosity declines with GEOGRAPHIC distance from Addis
// Ababa, but geography there is just a proxy for exactly this count of cumulative
// founder events — so rank is the causally-correct axis and needs no spatial
// structure, letting an aspatial OoA and an aspatial Babel scenario be scored on
// the same footing. The analysis (pkg/analysis/hetdistance.go) reads this via an
// optional interface, so the core.DemographyScheduler contract is unchanged.
func (s *Scenario) DemeRanks() map[int]int {
	ranks := make(map[int]int)
	if len(s.Epochs) == 0 {
		return ranks
	}
	// Origins: every deme active in the first epoch has rank 0.
	for _, spec := range s.Epochs[0].Demes {
		if _, ok := ranks[spec.ID]; !ok {
			ranks[spec.ID] = 0
		}
	}
	// Splits fire in epoch order (a parent is always founded before its child in a
	// well-formed scenario), so a single forward pass resolves every rank.
	for i := range s.Epochs {
		sp := s.Epochs[i].Split
		if sp == nil {
			continue
		}
		if _, exists := ranks[sp.Child]; exists {
			continue
		}
		parentRank, ok := ranks[sp.Parent]
		if !ok {
			parentRank = 0 // parent unseen (malformed order); treat it as an origin
		}
		ranks[sp.Child] = parentRank + 1
	}
	return ranks
}

// epochAt returns the epoch governing a DRIFT year (the last epoch that has begun);
// years past the end clamp to the final epoch (present-day structure persists).
func (s *Scenario) epochAt(year int) *Epoch {
	idx := 0
	for i := range s.Epochs {
		if year >= s.Epochs[i].startYear {
			idx = i
		}
	}
	return &s.Epochs[idx]
}

// capFor computes a deme's census cap at a given year: base = round(Ne/Q), grown
// exponentially from the epoch start when the deme has a growth rate. Year is
// clamped to the scenario end so growth doesn't run away past the present.
func (s *Scenario) capFor(spec sizeSpec, e *Epoch, year int) int {
	base := spec.Size / s.Rescale
	if spec.Growth != 0 {
		y := year
		if y > s.totalYears {
			y = s.totalYears
		}
		// Convert elapsed DRIFT years back to published generations, then grow.
		gensPublished := float64(y-e.startYear) / s.GenerationTime * s.Rescale
		if gensPublished < 0 {
			gensPublished = 0
		}
		base *= math.Exp(spec.Growth * gensPublished)
	}
	cap := int(math.Round(base))
	if cap < 1 {
		cap = 1
	}
	return cap
}

// perYearMigration converts a published per-generation rate to a per-year
// per-individual migration probability under the current rescaling.
func (s *Scenario) perYearMigration(rate float64) float64 {
	p := rate * s.Rescale / s.GenerationTime
	if p < 0 {
		p = 0
	}
	if p > 1 {
		p = 1
	}
	return p
}

// Apply advances the scenario for the current year: fire a split if one begins this
// year, migrate individuals per the epoch matrix, and publish per-deme caps.
func (s *Scenario) Apply(model *core.Model, pop *core.Pop) {
	year := model.FreeParameters["year"]
	e := s.epochAt(year)

	// A split begins its epoch on the epoch's first year (fires exactly once, since
	// year advances monotonically and resume continues past the boundary).
	if e.Split != nil && year == e.startYear {
		s.doSplit(pop, e.Split)
	}

	// Per-deme census caps for this year (read by Death).
	caps := make(map[int]int, len(e.Demes))
	for _, spec := range e.Demes {
		caps[spec.ID] = s.capFor(spec, e, year)
	}
	model.DemeCaps = caps

	// Migration per the epoch matrix.
	s.migrate(pop, e)
}

// doSplit relabels round(child_founders/Q) random members of the parent deme into
// the child deme (capped at half the parent so it is never emptied). RNG is
// consumed in sorted-id order for reproducibility.
func (s *Scenario) doSplit(pop *core.Pop, sp *splitSpec) {
	members := membersOf(pop, sp.Parent)
	if len(members) == 0 {
		fmt.Printf("demography: split parent deme %d is empty; child %d not founded\n", sp.Parent, sp.Child)
		return
	}
	n := int(math.Round(sp.ChildFounders / s.Rescale))
	if n < 1 {
		n = 1
	}
	if n > len(members)/2 {
		n = len(members) / 2
	}
	if n < 1 {
		n = 1 // parent has a single member; move it so the child at least exists
	}
	// Partial Fisher-Yates: pick n distinct members without replacement.
	for i := 0; i < n && i < len(members); i++ {
		j := i + utils.RandIntn(len(members)-i)
		members[i], members[j] = members[j], members[i]
		pop.IndData[members[i]][individual.Deme] = sp.Child
	}
	fmt.Printf("demography: split — moved %d individuals from deme %d to new deme %d\n", n, sp.Parent, sp.Child)
}

// migrate moves individuals between demes according to the epoch's symmetric
// migration matrix. Each individual is considered once; the first matching pair
// (in file order) that "hits" moves it. Membership is snapshotted via sorted ids
// and moves are applied afterward so within-year order can't cascade.
func (s *Scenario) migrate(pop *core.Pop, e *Epoch) {
	if len(e.Migration) == 0 {
		return
	}
	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		ids = append(ids, id)
	}
	sort.Ints(ids)

	moves := make(map[int]int)
	for _, id := range ids {
		d := pop.IndData[id][individual.Deme]
		for _, ms := range e.Migration {
			var other int
			if ms.A == d {
				other = ms.B
			} else if ms.B == d {
				other = ms.A
			} else {
				continue
			}
			if utils.RandFloat64() < s.perYearMigration(ms.Rate) {
				moves[id] = other
				break
			}
		}
	}
	// Apply in sorted order (deterministic; the map itself is order-independent).
	dests := make([]int, 0, len(moves))
	for id := range moves {
		dests = append(dests, id)
	}
	sort.Ints(dests)
	for _, id := range dests {
		pop.IndData[id][individual.Deme] = moves[id]
	}
}

// membersOf returns the sorted living-individual ids currently in a deme.
func membersOf(pop *core.Pop, deme int) []int {
	var out []int
	for id, data := range pop.IndData {
		if data[individual.Deme] == deme {
			out = append(out, id)
		}
	}
	sort.Ints(out)
	return out
}
