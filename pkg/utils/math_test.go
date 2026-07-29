package utils

import (
	"math"
	"testing"
)

// TestGammaRandomMoments checks the Marsaglia-Tsang gamma sampler recovers the
// analytic mean (shape*scale) and variance (shape*scale^2) for both shape>1 and
// the shape<1 boost path (the human-DFE regime, alpha~0.2). Powers the opt-in
// gamma DFE (§6f).
func TestGammaRandomMoments(t *testing.T) {
	cases := []struct{ shape, scale float64 }{
		{0.2, 0.05}, // shape<1 boost path; mean 0.01 (human deleterious DFE)
		{2.0, 1.5},  // shape>1 direct path
	}
	const nDraws = 200000
	for _, c := range cases {
		SeedRNG(4242)
		var sum, sumSq float64
		var negatives int
		for i := 0; i < nDraws; i++ {
			x := GammaRandom(c.shape, c.scale)
			if x < 0 {
				negatives++
			}
			sum += x
			sumSq += x * x
		}
		mean := sum / nDraws
		variance := sumSq/nDraws - mean*mean
		wantMean := c.shape * c.scale
		wantVar := c.shape * c.scale * c.scale
		if negatives > 0 {
			t.Errorf("shape=%.2g: %d negative draws (gamma is nonnegative)", c.shape, negatives)
		}
		if math.Abs(mean-wantMean) > 0.03*wantMean {
			t.Errorf("shape=%.2g scale=%.2g: mean=%.5g want ~%.5g", c.shape, c.scale, mean, wantMean)
		}
		if math.Abs(variance-wantVar) > 0.08*wantVar {
			t.Errorf("shape=%.2g scale=%.2g: var=%.5g want ~%.5g", c.shape, c.scale, variance, wantVar)
		}
	}
}

// TestGammaRandomGuards: nonpositive shape/scale return 0 (never NaN/panic).
func TestGammaRandomGuards(t *testing.T) {
	SeedRNG(1)
	for _, tc := range [][2]float64{{0, 1}, {1, 0}, {-1, 1}, {1, -1}} {
		if g := GammaRandom(tc[0], tc[1]); g != 0 {
			t.Errorf("GammaRandom(%.1f,%.1f)=%v, want 0", tc[0], tc[1], g)
		}
	}
}

// TestGammaRandomReproducible: same seed => identical stream.
func TestGammaRandomReproducible(t *testing.T) {
	SeedRNG(777)
	a := make([]float64, 50)
	for i := range a {
		a[i] = GammaRandom(0.2, 0.05)
	}
	SeedRNG(777)
	for i := range a {
		if b := GammaRandom(0.2, 0.05); b != a[i] {
			t.Fatalf("draw %d not reproducible: %v vs %v", i, a[i], b)
		}
	}
}

// TestNormalRandomMoments: standard normal has mean ~0, sd ~1.
func TestNormalRandomMoments(t *testing.T) {
	SeedRNG(99)
	const n = 200000
	var sum, sumSq float64
	for i := 0; i < n; i++ {
		x := NormalRandom()
		sum += x
		sumSq += x * x
	}
	mean := sum / n
	sd := math.Sqrt(sumSq/n - mean*mean)
	if math.Abs(mean) > 0.02 {
		t.Errorf("normal mean=%.4f want ~0", mean)
	}
	if math.Abs(sd-1) > 0.02 {
		t.Errorf("normal sd=%.4f want ~1", sd)
	}
}
