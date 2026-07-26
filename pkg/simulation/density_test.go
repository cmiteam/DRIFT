package simulation

import (
	"math"
	"testing"

	"drift/pkg/core"
)

func modelWithK(k float64) *core.Model {
	return &core.Model{
		Parameters:   map[string]float64{"carrying_capacity": k},
		StringParams: map[string]string{},
	}
}

func TestDensityRegulated(t *testing.T) {
	cases := []struct {
		name string
		k    float64
		want bool
	}{
		{"absent (zero-value) disables", 0, false},
		{"explicit zero disables", 0, false},
		{"negative disables", -1, false},
		{"positive enables", 5000, true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := modelWithK(tc.k)
			if got := densityRegulated(m); got != tc.want {
				t.Errorf("densityRegulated(k=%.0f) = %v, want %v", tc.k, got, tc.want)
			}
		})
	}
	// The zero-value case: a model whose params map has no carrying_capacity entry
	// at all must read as disabled (Go returns 0 for a missing key) — this is what
	// keeps existing models a strict no-op.
	empty := &core.Model{Parameters: map[string]float64{}, StringParams: map[string]string{}}
	if densityRegulated(empty) {
		t.Errorf("densityRegulated with no carrying_capacity key should be false")
	}
}

func TestLogisticGrowthCeiling(t *testing.T) {
	// g = 1.5 (exactly representable) => intrinsic rate r = 0.5. K = 2000 with N
	// chosen so every N/K and N*gEff product is dyadic (exact in float64), keeping
	// these expectations free of ceil() rounding noise.
	const g = 1.5
	const k = 2000.0
	cases := []struct {
		name string
		n    float64
		k    float64
		want int
	}{
		// Sparse population (N = K/2): density 0.5, gEff = 1 + 0.5*0.5 = 1.25 =>
		// 1000*1.25 = 1250. Below the constant ceiling ceil(1000*1.5)=1500.
		{"sparse grows below max rate", 1000, k, 1250},
		// At carrying capacity: gEff = 1 => holds steady at N.
		{"at K holds steady", 2000, k, 2000},
		// Above K: density -0.5, gEff = 0.75 => 3000*0.75 = 2250, a decline toward K.
		{"above K declines", 3000, k, 2250},
		// Large overshoot (N = 3K): density -2 clamped to -1 => gEff = 0.5 (decline
		// capped at r per year), 6000*0.5 = 3000, NOT 6000*(1+0.5*-2)=0 -> annihilation.
		{"overshoot decline clamped at r", 6000, k, 3000},
		// Floor: N=0 with births pending would compute ceil(0)=0; floored to 1 so the
		// logistic term alone never empties the population.
		{"floored at one", 0, k, 1},
		// Disabled (K <= 0): exactly the constant-rate ceiling ceil(N*g).
		{"disabled falls back to constant", 1000, 0, int(math.Ceil(1000 * g))},
		{"negative K falls back to constant", 250, -1, int(math.Ceil(250 * g))},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := logisticGrowthCeiling(tc.n, tc.k, g); got != tc.want {
				t.Errorf("logisticGrowthCeiling(n=%.0f, k=%.0f, g=%.2f) = %d, want %d",
					tc.n, tc.k, g, got, tc.want)
			}
		})
	}
}

// TestLogisticCeilingByteIdenticalWhenDisabled asserts the helper reproduces the
// exact constant-rate formula used by the original death.go Step 3 for K <= 0, so a
// run with carrying_capacity unset takes a byte-identical code path.
func TestLogisticCeilingByteIdenticalWhenDisabled(t *testing.T) {
	// Use DRIFT's real default growth rate.
	const g = 1.00462
	for _, n := range []float64{1, 10, 100, 999, 1000, 5000, 12345} {
		want := int(math.Ceil(n * g))
		if got := logisticGrowthCeiling(n, 0, g); got != want {
			t.Errorf("disabled ceiling(n=%.0f) = %d, want ceil(n*g)=%d", n, got, want)
		}
	}
}

// TestLogisticCeilingBitesBelowK is the distinguishable-outcome guard: with K
// active, the growth ceiling never exceeds the constant-rate ceiling (it throttles,
// never accelerates), and the *per-capita* growth allowance falls monotonically as
// the population fills toward K — the defining logistic property (the absolute yearly
// increment is humped, peaking near K/2, so per-capita, not headroom, is the monotone
// quantity). It confirms the term actually slows growth, it is not a silent no-op.
func TestLogisticCeilingBitesBelowK(t *testing.T) {
	const g = 1.00462
	const k = 10000.0
	constantCeil := func(n float64) int { return int(math.Ceil(n * g)) }

	prevPerCapita := math.Inf(1)
	sawThrottle := false
	for _, n := range []float64{1000, 3000, 6000, 9000, 9900} {
		logistic := logisticGrowthCeiling(n, k, g)
		constant := constantCeil(n)
		if logistic > constant {
			t.Errorf("at n=%.0f logistic ceiling %d exceeds constant %d (should throttle, not accelerate)", n, logistic, constant)
		}
		if logistic < constant {
			sawThrottle = true
		}
		// Per-capita growth allowance (fraction the population may grow this year)
		// must fall monotonically toward K. A tiny epsilon absorbs ceil() rounding.
		perCapita := (float64(logistic) - n) / n
		if perCapita > prevPerCapita+1e-9 {
			t.Errorf("at n=%.0f per-capita growth %.6f exceeds that at the previous (sparser) N (%.6f); per-capita growth must decelerate toward K", n, perCapita, prevPerCapita)
		}
		prevPerCapita = perCapita
	}
	if !sawThrottle {
		t.Errorf("logistic ceiling never fell below the constant ceiling — the term is a no-op below K")
	}

	// And at exactly K the allowed headroom is zero (equilibrium).
	if got := logisticGrowthCeiling(k, k, g); got != int(k) {
		t.Errorf("at N=K=%.0f ceiling = %d, want %d (no growth)", k, got, int(k))
	}
}
