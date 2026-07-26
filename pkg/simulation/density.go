package simulation

import (
	"drift/pkg/core"
	"math"
)

// Density dependence / logistic carrying capacity (roadmap §2).
//
// Without it, DRIFT's population is regulated only by hard caps in death.go: a
// flat max_pop_size ceiling and a constant max_growth_rate. So a population grows
// at a fixed rate every year until it slams into max_pop_size — pure
// exponential-growth-or-extinction with a sharp elbow, not a realistic S-curve.
//
// The logistic term makes the per-year growth *ceiling* density-dependent: the
// realized growth rate declines linearly to 1 as the population approaches the
// carrying capacity K, so (e.g.) a post-bottleneck population recovers along a
// smooth logistic curve toward K instead of exploding to the hard cap.
//
// It is opt-in via the "carrying_capacity" param (K). K <= 0 disables it entirely,
// leaving the original constant-max_growth_rate cull exactly in place — a strict
// no-op that keeps default/neutral runs byte-identical (§6h -validate still PASS).

// densityRegulated reports whether a logistic carrying-capacity term is active.
func densityRegulated(model *core.Model) bool {
	return model.Parameters["carrying_capacity"] > 0
}

// logisticGrowthCeiling returns the maximum allowed population size for the current
// year under a logistic carrying-capacity term. It generalizes death.go Step 3's
// constant max_growth_rate ceiling into a density-dependent one:
//
//	r     = maxGrowthRate - 1     // intrinsic per-year growth rate (the N->0 limit)
//	gEff  = 1 + r*(1 - N/K)       // realized growth factor: declines linearly to 1 at K
//	ceil  = ceil(N * gEff)
//
// where N is last year's population size. So a sparse population (N << K) still grows
// at ~maxGrowthRate, one at carrying capacity (N == K) holds steady (gEff == 1), and
// one that overshoots (N > K) declines gently back toward K. This is the standard
// discrete logistic map N_{t+1} = N_t * (1 + r(1 - N_t/K)).
//
// Two guards keep it well-behaved: the density factor (1 - N/K) is clamped at -1 so
// an overshoot can never shrink the population faster than r per year (symmetric with
// the max growth rate — no annihilation crash if K is set below the current size),
// and the ceiling is floored at 1 so the logistic term alone never empties the
// population. The absolute max_pop_size hard cap is applied by the caller on top of
// this (it remains a safety ceiling above K).
func logisticGrowthCeiling(n, k, maxGrowthRate float64) int {
	if k <= 0 {
		// Disabled: fall back to the constant-rate ceiling (caller normally guards
		// on densityRegulated, but keep this self-consistent).
		return int(math.Ceil(n * maxGrowthRate))
	}
	r := maxGrowthRate - 1.0
	density := 1.0 - n/k
	if density < -1.0 {
		density = -1.0
	}
	gEff := 1.0 + r*density
	allowed := int(math.Ceil(n * gEff))
	if allowed < 1 {
		allowed = 1
	}
	return allowed
}
