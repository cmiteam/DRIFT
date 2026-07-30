package core

import (
	"strconv"
	"strings"
)

// Movement barriers & corridors (roadmap §3).
//
// DRIFT's only movement primitive is utils.Wander: it draws a random destination
// inside a (2·wander+1)² Chebyshev box centered on the individual and accepts the
// move iff the destination cell is habitable. Because that is a *teleport* — the
// destination is tested but no cell between origin and destination is — a strip of
// impassable terrain narrower than `wander` is simply jumped over: a thin strait or
// mountain ridge does not actually block crossing today. There is no way to make a
// barrier that stops crossing without also making the far side an invalid destination.
//
// Movement barriers close that gap with a *crossing test*. A new opt-in param,
// `movement_barriers`, lists the terrain codes that block crossing (as opposed to
// habitat_suitability, which governs where an individual may *live*). A move is then
// allowed only if the straight cell-line from origin to destination does not pass
// through a barrier cell (line-of-sight; see CanTraverse). This makes a one-cell water
// strip a real wall — to reach the far side an individual must Wander around the tip
// of the barrier (if it is within range) or through a *corridor*: a passable strip the
// map author draws through the barrier (e.g. a land bridge). Corridors therefore need
// no new geometry — any non-barrier cell restores line-of-sight through it.
//
// The same crossing test gates the other place an individual moves across space: the
// spatial mating modules (mating_style distance / age_distance) pair a man with women
// by pure Manhattan proximity through an influence grid that ignores terrain, then move
// the woman to the man's cell. Without the test a barrier narrower than
// max_mating_distance would still let couples form straight through it, so gene flow —
// the whole point of a barrier — would leak across. filterTraversableWomen (pkg/methods)
// therefore drops any woman a barrier separates from the man, the mating analogue of
// Wander's dispersal block. Random (aspatial) mating is left untouched.
//
// In Wander the test is a pure accept/reject on the already-drawn destination, so it
// draws no additional RNG: Wander still makes its two offset draws in the same order,
// and when the table is inactive CanTraverse returns true unconditionally, leaving the
// move-accept predicate — and the whole run — byte-identical. The mating filter is
// likewise a strict no-op when inactive (it returns the candidate slice unchanged).
// Barriers separate the "can I cross here" axis from habitat_suitability's "can I live
// here" axis; the two params compose independently (a wide ocean is both uninhabitable
// and a barrier; an open plain can be uninhabitable-but-passable; a mountain pass
// habitable-and-open).
//
// The accessor is deliberately terrain-derived, mirroring the habitat_suitability and
// §6a genetic-map accessors: a finer-grained per-cell passability grid or explicit
// barrier geometry can slot in behind CanTraverse later without touching the callers.

// barrierPassableTokens are the values that, when a terrain code is written as
// "<code>:<value>", mark the code as passable rather than a barrier. Any other value
// (or a bare "<code>") marks the code as an impassable barrier. Lets a table document
// intent (e.g. "3:block") while a bare "3" means the same thing.
var barrierPassableTokens = map[string]bool{
	"pass": true, "passable": true, "open": true, "none": true, "0": true, "false": true,
}

// BarrierTable is the parsed movement_barriers param: the set of terrain codes that
// block crossing. Active is false when the param was empty or unparseable, in which
// case CanTraverse always returns true (no crossing constraint).
type BarrierTable struct {
	Active bool
	Block  map[int]bool // terrain code -> blocks crossing
}

// parseBarrierTable parses a movement_barriers spec: terrain codes separated by
// space/';'/',', each optionally suffixed ":<value>" (":block" for clarity, or a
// passable token like ":pass" to explicitly mark a code non-blocking). A bare code is
// a barrier. Malformed tokens are skipped rather than failing the load. An empty
// result (no blocking codes) is inactive.
func parseBarrierTable(spec string) *BarrierTable {
	tokens := strings.FieldsFunc(spec, func(r rune) bool {
		return r == ' ' || r == '\t' || r == '\n' || r == '\r' || r == ';' || r == ','
	})
	block := make(map[int]bool)
	for _, tok := range tokens {
		codeStr := tok
		blocks := true
		if i := strings.IndexByte(tok, ':'); i >= 0 {
			codeStr = tok[:i]
			val := strings.ToLower(strings.TrimSpace(tok[i+1:]))
			if barrierPassableTokens[val] {
				blocks = false
			}
		}
		code, err := strconv.Atoi(strings.TrimSpace(codeStr))
		if err != nil {
			continue // not a terrain code
		}
		if blocks {
			block[code] = true
		}
	}
	if len(block) == 0 {
		return &BarrierTable{Active: false}
	}
	return &BarrierTable{Active: true, Block: block}
}

// ensureBarriers lazily builds the barrier table from the param on first use and
// caches it. Idempotent and cheap (a nil check after the first call), so the hot
// caller (Wander, per birth) can call it freely. Building from the param here — rather
// than at model load — means every load path gets consistent behavior without wiring.
func (m *Model) ensureBarriers() {
	if m.Barriers != nil {
		return
	}
	m.Barriers = parseBarrierTable(m.StringParam("movement_barriers", ""))
}

// MovementBarriersActive reports whether a movement-barrier table is configured. When
// false, CanTraverse imposes no constraint and Wander's accept predicate is unchanged.
func (m *Model) MovementBarriersActive() bool {
	m.ensureBarriers()
	return m.Barriers.Active
}

// IsBarrier reports whether the given cell's terrain blocks crossing. Out-of-bounds /
// no-map cells and an inactive table are never barriers.
func (m *Model) IsBarrier(lat, lon int) bool {
	m.ensureBarriers()
	if !m.Barriers.Active {
		return false
	}
	if lat < 0 || lon < 0 || lat >= len(m.Map) {
		return false
	}
	if lon >= len(m.Map[lat]) {
		return false
	}
	return m.Barriers.Block[m.Map[lat][lon]]
}

// CanTraverse reports whether a move from (fromLat, fromLon) to (toLat, toLon) is
// allowed by the barrier layer: the straight cell-line between the two endpoints must
// not pass through a barrier cell (line-of-sight). The two endpoints themselves are
// not tested — the origin is where the individual already lives, and the destination
// is separately gated by IsHabitable — so only the cells *crossed* matter.
//
// With no barrier table this returns true immediately (zero work, byte-identical
// move-accept). The line is walked with an integer Bresenham traversal, at most
// ~2·wander cells, so the per-birth cost is negligible next to meiosis.
func (m *Model) CanTraverse(fromLat, fromLon, toLat, toLon int) bool {
	m.ensureBarriers()
	if !m.Barriers.Active {
		return true
	}
	// Bresenham line from origin to destination; check only the interior cells.
	dLat := abs(toLat - fromLat)
	dLon := abs(toLon - fromLon)
	sLat := sign(toLat - fromLat)
	sLon := sign(toLon - fromLon)
	err := dLat - dLon
	lat, lon := fromLat, fromLon
	for {
		if lat == toLat && lon == toLon {
			return true
		}
		if !(lat == fromLat && lon == fromLon) {
			// Interior cell — a barrier here blocks the crossing.
			if m.IsBarrier(lat, lon) {
				return false
			}
		}
		e2 := 2 * err
		if e2 > -dLon {
			err -= dLon
			lat += sLat
		}
		if e2 < dLat {
			err += dLat
			lon += sLon
		}
	}
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func sign(x int) int {
	switch {
	case x > 0:
		return 1
	case x < 0:
		return -1
	default:
		return 0
	}
}
