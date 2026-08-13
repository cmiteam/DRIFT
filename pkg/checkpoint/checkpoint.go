// Package checkpoint serializes and restores the full mutable state of a
// simulation run, enabling stop/resume and "fork" experiments: run to year T,
// checkpoint that population as a baseline, then continue it — exactly (same RNG
// stream) or divergently (re-seeded) — while holding the shared history constant.
//
// The blob contains state only: the population, the FreeParameters counters, and
// the RNG position. Model configuration (parameters, chromosomes, map, animation
// state) is reloaded from the model's files on resume, then the checkpoint is
// overlaid.
package checkpoint

import (
	"drift/pkg/core"
	"drift/pkg/utils"
	"encoding/gob"
	"fmt"
	"os"
	"time"
)

// SchemaVersion is bumped whenever the Checkpoint layout changes so that old
// blobs fail loudly instead of decoding into a mismatched struct.
//
// v2 (§6b): added the Deme field to the IndData enum, lengthening each
// individual's serialized []int; v1 saves have shorter slices and would
// mis-index, so they are rejected on load.
const SchemaVersion = 2

// Checkpoint is the serialized state of a run at a point in time.
type Checkpoint struct {
	SchemaVersion  int
	SavedAt        string // RFC3339 timestamp, for provenance
	ModelName      string
	Scenario       string
	RNGSeed        int64  // the seed the run was started from
	RNGState       []byte // exact generator position (see utils.RNGState)
	Year           int    // simulation year at capture
	Run            int    // run index at capture
	FreeParameters map[string]int
	Pop            *core.Pop
}

// Save writes a checkpoint of the current simulation state to path.
func Save(path string, model *core.Model, pop *core.Pop) error {
	cp := &Checkpoint{
		SchemaVersion:  SchemaVersion,
		SavedAt:        time.Now().Format(time.RFC3339),
		ModelName:      model.ModelName,
		Scenario:       model.Scenario,
		RNGSeed:        utils.CurrentSeed(),
		RNGState:       utils.RNGState(),
		Year:           model.FreeParameters["year"],
		Run:            model.FreeParameters["run"],
		FreeParameters: copyIntMap(model.FreeParameters),
		Pop:            pop,
	}

	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create checkpoint %q: %w", path, err)
	}
	defer f.Close()
	if err := gob.NewEncoder(f).Encode(cp); err != nil {
		return fmt.Errorf("encode checkpoint: %w", err)
	}
	return nil
}

// Load reads and validates a checkpoint from path.
func Load(path string) (*Checkpoint, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open checkpoint %q: %w", path, err)
	}
	defer f.Close()

	cp := &Checkpoint{}
	if err := gob.NewDecoder(f).Decode(cp); err != nil {
		return nil, fmt.Errorf("decode checkpoint: %w", err)
	}
	if cp.SchemaVersion != SchemaVersion {
		return nil, fmt.Errorf("checkpoint schema version %d not supported (this build expects %d)",
			cp.SchemaVersion, SchemaVersion)
	}
	if cp.Pop != nil {
		ensurePopMaps(cp.Pop)
	}
	return cp, nil
}

// Restore applies a checkpoint onto a freshly-loaded model: overlays the
// FreeParameters counters and (unless forking) restores the exact RNG position.
// If forkSeed != 0 the RNG is re-seeded with it instead, producing a divergent
// continuation from the same shared history. The caller uses cp.Pop as the live
// population and resumes from cp.Year+1.
func Restore(model *core.Model, cp *Checkpoint, forkSeed int64) error {
	for k, v := range cp.FreeParameters {
		model.FreeParameters[k] = v
	}
	if forkSeed != 0 {
		utils.SeedRNG(forkSeed)
		return nil
	}
	return utils.RestoreRNG(cp.RNGState)
}

// ensurePopMaps initializes any nil maps so a resumed population can be written
// to without a nil-map panic (gob decodes empty maps as nil).
func ensurePopMaps(p *core.Pop) {
	if p.IndData == nil {
		p.IndData = map[int][]int{}
	}
	if p.Chromosomes == nil {
		p.Chromosomes = map[int][][]uint64{}
	}
	if p.Centromeres == nil {
		p.Centromeres = map[int][2][]uint64{}
	}
	if p.IndMutations == nil {
		p.IndMutations = map[int]map[int][]int{}
	}
	if p.MutationPool == nil {
		p.MutationPool = map[int]core.Mutation{}
	}
	if p.MutationHist == nil {
		p.MutationHist = map[int]int{}
	}
	if p.Tracking == nil {
		p.Tracking = map[string]int{}
	}
	if p.AlleleFreqs == nil {
		p.AlleleFreqs = map[int][]int16{}
	}
	if p.SFSHistory == nil {
		p.SFSHistory = map[int][]int{}
	}
	if p.NeHistory == nil {
		p.NeHistory = map[int]*core.NeSnapshot{}
	}
	if p.LoadHistory == nil {
		p.LoadHistory = map[int]*core.LoadSnapshot{}
	}
	if p.MaleDB == nil {
		p.MaleDB = map[int]core.Ancestor{}
	}
	// Pedigree is deliberately NOT force-allocated: nil means "track_pedigree off",
	// and the Ped* methods allocate lazily on first write, so a checkpoint taken with
	// the flag off can be resumed with it on (later births are recorded, earlier
	// ancestry reports as censored) without pretending an empty pedigree is a complete one.
	if p.FemaleDB == nil {
		p.FemaleDB = map[int]core.Ancestor{}
	}
}

func copyIntMap(m map[string]int) map[string]int {
	out := make(map[string]int, len(m))
	for k, v := range m {
		out[k] = v
	}
	return out
}
