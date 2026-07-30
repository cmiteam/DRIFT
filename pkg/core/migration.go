package core

import (
	"sort"
	"strconv"
	"strings"
)

// Deme migration matrix (roadmap §3, "deme / island models with migration matrices").
//
// DRIFT's always-on structured-population channel is the scalar island model: with
// num_demes > 1 and deme_migration_rate = m, migrateDemes (pkg/simulation) moves each
// individual to a UNIFORMLY-RANDOM *other* deme with per-year probability m. Every deme
// pair therefore exchanges migrants at the identical symmetric rate — a Wright island
// model with a single knob. The only per-pair rates that exist otherwise live inside a
// §6c demography.json scenario, and those are symmetric AND tied to the epoch schedule
// (they only run when a DemographyScheduler is attached).
//
// The migration matrix generalizes the scalar island model into a general, ASYMMETRIC,
// per-pair, always-on N-deme matrix that needs no scenario. A new opt-in param,
// `migration_matrix`, lists directional edges "<from>><to>:<rate>", e.g.
// "0>1:0.01 1>0:0.005 1>2:0.02": rate is the per-year per-individual probability that a
// resident of deme <from> migrates to deme <to> (matching the units of the scalar
// deme_migration_rate it generalizes, so it is a drop-in — the island model with rate m
// over N demes is exactly the matrix with every off-diagonal edge = m/(N-1)). This makes
// stepping-stone (only adjacent edges), source→sink (asymmetric edges), and a post-Babel
// radiation (edges out of one origin deme) expressible directly; the resulting migration
// ⇄ drift balance is what the Fst / f-stats / joint-SFS readouts (which already partition
// on the Deme label) measure.
//
// The application (migrateDemesMatrix, pkg/simulation) draws exactly one RandFloat64 per
// individual and walks that deme's outgoing edges in a fixed (Dest-sorted) order, so the
// stream stays byte-reproducible under a fixed rng_seed. Per-source outgoing rates should
// sum to <= 1 (the residual is the stay-put probability); a source whose edges sum above
// 1 simply never stays put (the cumulative walk always hits an edge), which is a config
// error the parser does not reject — it keeps the malformed-token-skipping, never-fail
// parse contract the habitat / barriers accessors use.
//
// When `migration_matrix` is empty/absent the table is inactive and the island path takes
// the original scalar migrateDemes unchanged, so an unset param is byte-identical on every
// path. The accessor is data-only (no RNG) and cached on the model, exactly like the
// habitat, barriers, and §6a genetic-map accessors, so a future richer source (a full
// dense matrix file, distance-decay kernel) can slot in behind MigrationMatrix() later
// without touching the callers.

// MigrationEdge is one directional migration rate: a resident of the source deme moves
// to Dest with per-year probability Rate.
type MigrationEdge struct {
	Dest int
	Rate float64
}

// MigrationMatrix is the parsed migration_matrix param: for each source deme, the sorted
// list of outgoing directional edges. Active is false when the param was empty or held no
// usable edges, in which case the island path falls back to the scalar deme_migration_rate.
type MigrationMatrix struct {
	Active bool
	Out    map[int][]MigrationEdge // source deme -> outgoing edges, sorted by Dest
}

// parseMigrationMatrix parses a migration_matrix spec: directional edges
// "<from>><to>:<rate>" separated by space/';'/','. A rate is the per-year per-individual
// probability of that from→to move. Malformed tokens are skipped (never fail the load),
// negative rates are clamped to 0, self-edges (from == to) are dropped, and duplicate
// from→to edges keep the last value. Edges are sorted by Dest per source so the RNG walk
// is deterministic. An empty result is inactive.
func parseMigrationMatrix(spec string) *MigrationMatrix {
	tokens := strings.FieldsFunc(spec, func(r rune) bool {
		return r == ' ' || r == '\t' || r == '\n' || r == '\r' || r == ';' || r == ','
	})
	// source -> dest -> rate (map dedups repeated edges, last wins)
	tmp := make(map[int]map[int]float64)
	for _, tok := range tokens {
		arrow := strings.IndexByte(tok, '>')
		colon := strings.IndexByte(tok, ':')
		if arrow < 0 || colon < 0 || colon < arrow {
			continue // needs both "<from>><to>" and ":<rate>"
		}
		from, err1 := strconv.Atoi(strings.TrimSpace(tok[:arrow]))
		to, err2 := strconv.Atoi(strings.TrimSpace(tok[arrow+1 : colon]))
		rate, err3 := strconv.ParseFloat(strings.TrimSpace(tok[colon+1:]), 64)
		if err1 != nil || err2 != nil || err3 != nil {
			continue
		}
		if from == to {
			continue // a self-edge is not migration
		}
		if rate < 0 {
			rate = 0
		}
		if tmp[from] == nil {
			tmp[from] = make(map[int]float64)
		}
		tmp[from][to] = rate
	}

	out := make(map[int][]MigrationEdge, len(tmp))
	for from, dests := range tmp {
		edges := make([]MigrationEdge, 0, len(dests))
		for to, rate := range dests {
			if rate <= 0 {
				continue // a zero-rate edge is a no-op; drop it
			}
			edges = append(edges, MigrationEdge{Dest: to, Rate: rate})
		}
		if len(edges) == 0 {
			continue
		}
		sort.Slice(edges, func(i, j int) bool { return edges[i].Dest < edges[j].Dest })
		out[from] = edges
	}
	if len(out) == 0 {
		return &MigrationMatrix{Active: false}
	}
	return &MigrationMatrix{Active: true, Out: out}
}

// ensureMigration lazily builds the migration matrix from the param on first use and
// caches it. Idempotent and cheap after the first call, so the per-year Mating caller can
// consult it freely. Building from the param here — rather than at model load — means
// every load path gets consistent behavior without extra wiring (mirrors ensureHabitat /
// ensureBarriers).
func (m *Model) ensureMigration() {
	if m.Migration != nil {
		return
	}
	m.Migration = parseMigrationMatrix(m.StringParam("migration_matrix", ""))
}

// MigrationMatrix returns the parsed migration matrix (lazily built + cached). When the
// returned matrix is inactive the island path uses the scalar deme_migration_rate instead,
// so an unset param leaves the run byte-identical.
func (m *Model) MigrationMatrix() *MigrationMatrix {
	m.ensureMigration()
	return m.Migration
}
