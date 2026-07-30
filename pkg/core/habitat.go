package core

import (
	"strconv"
	"strings"
)

// Habitat suitability maps (roadmap §3).
//
// A DRIFT map cell today carries only a terrain code (see pkg/utils: Land=1,
// CoastalWater=2, OpenWater=3, HighMountain=4, Desert=5, Ice=6), and only Land is
// habitable — every accessor treats "on land" as a hard binary. Habitat suitability
// generalizes that binary into a per-cell multiplier in [0, +inf): a suitability of
// 1 is fully hospitable, 0.3 is a harsh-but-livable cell, and 0 is uninhabitable.
//
// The table is authored as a single opt-in param, `habitat_suitability`, mapping
// terrain codes to suitability multipliers, e.g. "1:1.0 4:0.2 5:0.3 6:0.05" (Land
// fully suitable, mountains/desert/ice progressively harsher). Tokens may be
// separated by spaces, ';' or ',' (so it survives a CSV parameter file without
// quoting). A terrain code not listed keeps its default: Land stays fully suitable
// (1.0), every other code stays uninhabitable (0.0). An empty/absent param leaves
// the table inactive, and the accessors then reproduce the original binary
// land/water behavior exactly — so an unset param is byte-identical on every path.
//
// Two levers read the table (both in pkg/simulation, both gated on HabitatActive so
// they never touch a run with no habitat param): a per-cell mortality multiplier on
// the actuarial death risk, and a suitability-weighted growth-ceiling cull (a "soft"
// spatial carrying capacity — the global logistic/cap still decides how many die,
// habitat only decides who). The accessor is deliberately terrain-derived; a
// finer-grained per-cell float grid can slot in behind CellSuitability later without
// touching the callers (mirrors the §6a genetic-map accessor).

// landTerrain is the terrain code for habitable land. It mirrors utils.Land (=1);
// duplicated as a bare const here because pkg/core cannot import pkg/utils (the
// dependency runs the other way). Kept in lockstep with the utils.Terrain enum.
const landTerrain = 1

// HabitatTable is the parsed habitat_suitability param: a terrain-code -> suitability
// multiplier lookup. Active is false when the param was empty or unparseable, in
// which case the cell accessors fall back to binary land/water semantics.
type HabitatTable struct {
	Active bool
	Mult   map[int]float64 // terrain code -> suitability multiplier
}

// parseHabitatTable parses a habitat_suitability spec ("<code>:<mult> ...", tokens
// separated by space/';'/','). Malformed tokens and negative multipliers are
// skipped/clamped rather than failing the load. An empty result is inactive.
func parseHabitatTable(spec string) *HabitatTable {
	tokens := strings.FieldsFunc(spec, func(r rune) bool {
		return r == ' ' || r == '\t' || r == '\n' || r == '\r' || r == ';' || r == ','
	})
	mult := make(map[int]float64)
	for _, tok := range tokens {
		kv := strings.SplitN(tok, ":", 2)
		if len(kv) != 2 {
			continue // not a code:mult token
		}
		code, err1 := strconv.Atoi(strings.TrimSpace(kv[0]))
		val, err2 := strconv.ParseFloat(strings.TrimSpace(kv[1]), 64)
		if err1 != nil || err2 != nil {
			continue
		}
		if val < 0 {
			val = 0 // an uninhabitable cell, not a negative suitability
		}
		mult[code] = val
	}
	if len(mult) == 0 {
		return &HabitatTable{Active: false}
	}
	return &HabitatTable{Active: true, Mult: mult}
}

// ensureHabitat lazily builds the habitat table from the param on first use and
// caches it. Idempotent and cheap (a nil check after the first call), so hot callers
// (Wander, per-individual death) can call it freely. Building from the param here —
// rather than at model load — means every load path (user-model, validation,
// base-model) gets consistent behavior without extra wiring.
func (m *Model) ensureHabitat() {
	if m.Habitat != nil {
		return
	}
	m.Habitat = parseHabitatTable(m.StringParam("habitat_suitability", ""))
}

// HabitatActive reports whether a habitat suitability table is configured. When
// false, the cell accessors reproduce the original binary land/water behavior and
// the §3 death-side levers are skipped entirely (byte-identical).
func (m *Model) HabitatActive() bool {
	m.ensureHabitat()
	return m.Habitat.Active
}

// CellSuitability returns the suitability multiplier for the given map cell:
//   - out of bounds / no map loaded => 0 (uninhabitable)
//   - table inactive => 1.0 on Land, 0.0 elsewhere (the original binary semantics)
//   - table active   => the configured multiplier for that terrain code, defaulting
//     to 1.0 for unlisted Land and 0.0 for any other unlisted terrain.
func (m *Model) CellSuitability(lat, lon int) float64 {
	if lat < 0 || lon < 0 || lat >= len(m.Map) {
		return 0
	}
	if lon >= len(m.Map[lat]) {
		return 0
	}
	terrain := m.Map[lat][lon]
	m.ensureHabitat()
	if m.Habitat.Active {
		if v, ok := m.Habitat.Mult[terrain]; ok {
			return v
		}
	}
	// Unlisted terrain, or inactive table: keep Land fully suitable and everything
	// else uninhabitable — identical to the pre-habitat binary IsLand behavior.
	if terrain == landTerrain {
		return 1.0
	}
	return 0.0
}

// IsHabitable reports whether an individual can occupy the given cell (suitability
// > 0). With no habitat table this is exactly utils.IsLand (Land only); with a table
// it also admits harsh terrains that were given a positive suitability.
func (m *Model) IsHabitable(lat, lon int) bool {
	return m.CellSuitability(lat, lon) > 0
}
