package utils

import (
	"drift/pkg/core"
	"testing"
)

// Same seed must reproduce the same sequence across all primitives.
func TestSeedRNG_Reproducible(t *testing.T) {
	seq := func(seed int64) []float64 {
		SeedRNG(seed)
		out := make([]float64, 0, 30)
		for i := 0; i < 10; i++ {
			out = append(out, float64(RandIntn(1000)))
			out = append(out, RandFloat64())
			out = append(out, float64(RandPoisson(5)))
		}
		return out
	}

	a := seq(42)
	b := seq(42)
	for i := range a {
		if a[i] != b[i] {
			t.Fatalf("same seed diverged at index %d: %v vs %v", i, a[i], b[i])
		}
	}
}

func TestSeedRNG_DifferentSeedsDiffer(t *testing.T) {
	SeedRNG(1)
	a := RandFloat64()
	SeedRNG(2)
	b := RandFloat64()
	if a == b {
		t.Errorf("different seeds produced identical first draw: %v", a)
	}
}

func TestRandPoisson_Basics(t *testing.T) {
	SeedRNG(7)
	if got := RandPoisson(0); got != 0 {
		t.Errorf("RandPoisson(0) = %d, want 0", got)
	}
	if got := RandPoisson(-3); got != 0 {
		t.Errorf("RandPoisson(negative) = %d, want 0", got)
	}
	// Mean of many draws should be near lambda for both the Knuth and
	// normal-approximation branches.
	for _, lambda := range []float64{5, 50} {
		SeedRNG(11)
		sum := 0
		const n = 20000
		for i := 0; i < n; i++ {
			sum += RandPoisson(lambda)
		}
		mean := float64(sum) / n
		if mean < lambda*0.9 || mean > lambda*1.1 {
			t.Errorf("RandPoisson(%.0f) empirical mean %.2f outside [%.1f, %.1f]",
				lambda, mean, lambda*0.9, lambda*1.1)
		}
	}
}

func TestRNGStateRoundTrip(t *testing.T) {
	SeedRNG(12345)
	// Advance the stream a bit, then snapshot.
	for i := 0; i < 17; i++ {
		RandFloat64()
	}
	state := RNGState()
	if len(state) != 32 {
		t.Fatalf("RNGState length = %d, want 32", len(state))
	}

	// Record the next values from the snapshot point.
	want := []float64{RandFloat64(), RandFloat64(), RandFloat64()}

	// Perturb the generator, then restore and confirm we get the same values.
	SeedRNG(999)
	RandFloat64()
	if err := RestoreRNG(state); err != nil {
		t.Fatalf("RestoreRNG: %v", err)
	}
	got := []float64{RandFloat64(), RandFloat64(), RandFloat64()}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("post-restore value %d = %v, want %v", i, got[i], want[i])
		}
	}

	if err := RestoreRNG([]byte{1, 2, 3}); err == nil {
		t.Errorf("RestoreRNG with bad length should error")
	}
}

func TestSetGenomeBits(t *testing.T) {
	model := &core.Model{
		FreeParameters: map[string]int{},
		ChromosomeArms: map[int]map[int][]int{
			0: {0: {0, 100}, 1: {100, 50}}, // highest end = 150
			1: {0: {150, 200}, 1: {350, 25}}, // highest end = 375
		},
	}
	if got := SetGenomeBits(model); got != 375 {
		t.Errorf("SetGenomeBits returned %d, want 375", got)
	}
	if model.FreeParameters["genome_bits"] != 375 {
		t.Errorf("genome_bits = %d, want 375", model.FreeParameters["genome_bits"])
	}
}

// Regression: SetGenomeBits must initialize a nil FreeParameters map rather than panic.
func TestSetGenomeBits_NilFreeParameters(t *testing.T) {
	model := &core.Model{
		ChromosomeArms: map[int]map[int][]int{0: {0: {0, 64}, 1: {64, 64}}},
	}
	if got := SetGenomeBits(model); got != 128 {
		t.Errorf("SetGenomeBits returned %d, want 128", got)
	}
}
