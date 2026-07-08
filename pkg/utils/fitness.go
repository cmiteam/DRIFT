package utils

import (
	"math"

	"drift/pkg/core"
)

// Fitness models (roadmap §1: epistasis / synergistic load) select how the
// per-locus fitness contributions of an individual's mutations combine into a
// single relative fitness. Set via the string param "fitness_model". This is a
// live debate in the genetic-entropy literature — whether load is additive,
// multiplicative, or synergistic materially changes the meltdown prediction — so
// the three models are toggleable and directly comparable in the same engine.
//
// Let c_i be locus i's contribution (Effect if homozygous derived, h*Effect if
// heterozygous; negative for deleterious, from CountFitnessAndMutations) and
// L = Σ c_i the additive load.
//
//   - additive       w = 1 + L                 (default; the historical
//     behavior — a strict no-op equal to the old fitness = 1 + load).
//   - multiplicative w = Π (1 + c_i)           (loci act independently, the
//     classic multiplicative fitness model; for small effects ≈ additive but
//     always ≥ additive for deleterious load, since it can't cross zero).
//   - synergistic    w = 1 + L + β·L·|L|       (deleterious load ACCELERATES —
//     each further mutation costs more than the last, the genetic-entropy
//     "synergistic epistasis" claim). β = epistasis_coefficient; β = 0 reduces
//     it exactly to additive. For a uniform per-locus effect s and n mutations
//     (L = n·s) this is w = 1 + n·s − β·n²·s², the standard quadratic
//     (Crow/Kondrashov) synergistic form, i.e. accelerating with mutation count.
//
// NEUTRAL NO-OP (preserved for the §6h validation harness): with no non-neutral
// mutations every c_i = 0, so L = 0 and all three models return exactly 1.0 —
// no divergence and no RNG drawn. The default "additive" is byte-identical to
// the pre-toggle behavior, so only a non-neutral, non-additive run changes.
const (
	FitAdditive       = "additive"
	FitMultiplicative = "multiplicative"
	FitSynergistic    = "synergistic"
)

// FitnessModel returns the normalized fitness model for this model, falling back
// to additive for an unset or unrecognized value (so a typo can never silently
// disable load).
func FitnessModel(model *core.Model) string {
	switch m := model.StringParam("fitness_model", FitAdditive); m {
	case FitAdditive, FitMultiplicative, FitSynergistic:
		return m
	default:
		return FitAdditive
	}
}

// CombineFitness combines per-locus fitness contributions into a single relative
// fitness according to the model's fitness_model param. contribs is expected in
// a deterministic (sorted) order so the floating-point result is reproducible
// under a fixed seed. An empty slice yields exactly 1.0 in every model.
func CombineFitness(contribs []float64, model *core.Model) float64 {
	switch FitnessModel(model) {
	case FitMultiplicative:
		w := 1.0
		for _, c := range contribs {
			w *= 1.0 + c
		}
		return w
	case FitSynergistic:
		load := 0.0
		for _, c := range contribs {
			load += c
		}
		beta := model.Parameters["epistasis_coefficient"]
		return 1.0 + load + beta*load*math.Abs(load)
	default: // FitAdditive: byte-identical to the historical fitness = 1 + load.
		load := 0.0
		for _, c := range contribs {
			load += c
		}
		return 1.0 + load
	}
}
