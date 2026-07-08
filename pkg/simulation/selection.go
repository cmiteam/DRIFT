package simulation

import "drift/pkg/core"

// Selection modes control how an individual's mutation fitness load acts on the
// life cycle (roadmap §2). Set via the string param "selection_mode".
//
//   - none       fitness affects neither birth nor survival (pure drift, even
//     with non-neutral mutations — useful for isolating drift).
//   - fecundity  fitness scales the per-year birth probability only (default).
//   - viability  fitness scales the per-year survival (death hazard) only.
//   - both       fitness acts on both fecundity and viability.
//
// Fitness selection can only act at all when mutations are tracked (otherwise
// there is no per-genotype load), so selectionMode collapses to "none" when
// track_mutations is off. The default is "fecundity", matching the historical
// behavior where fitness gated birth probability. NOTE: a strictly neutral run
// has fitness = 1 for everyone, so every mode yields byte-identical results —
// the mode only bites once mutations carry non-zero effects.
const (
	selNone      = "none"
	selFecundity = "fecundity"
	selViability = "viability"
	selBoth      = "both"
)

// selectionMode returns the normalized selection mode for this model.
func selectionMode(model *core.Model) string {
	if model.Parameters["track_mutations"] != 1 {
		return selNone
	}
	switch m := model.StringParam("selection_mode", selFecundity); m {
	case selNone, selFecundity, selViability, selBoth:
		return m
	default:
		// Unrecognized value: fall back to the historical fecundity behavior
		// rather than silently disabling selection.
		return selFecundity
	}
}

// fecunditySelection reports whether fitness should scale birth probability.
func fecunditySelection(model *core.Model) bool {
	m := selectionMode(model)
	return m == selFecundity || m == selBoth
}

// viabilitySelection reports whether fitness should scale the death hazard.
func viabilitySelection(model *core.Model) bool {
	m := selectionMode(model)
	return m == selViability || m == selBoth
}

// viabilityHazardFactor maps an individual's relative fitness to a multiplicative
// factor on its baseline death risk under viability selection:
//
//	fitness = 1 (no load)          -> 1   (baseline hazard, no effect)
//	fitness < 1 (deleterious load) -> > 1 (higher hazard, dies more)
//	fitness > 1 (beneficial load)  -> < 1 (lower hazard, survives more)
//
// i.e. factor = 2 - fitness = 1 - load, the natural inverse of the fecundity
// side (where birth is accepted with probability = fitness). Floored at 0 so a
// pathological super-deleterious genotype (fitness > 2) can never produce a
// negative risk. This replaces the earlier `deathrisk * fitness`, which was
// inverted (it made less-fit individuals die *less*).
func viabilityHazardFactor(fitness float64) float64 {
	factor := 2.0 - fitness
	if factor < 0 {
		return 0
	}
	return factor
}
