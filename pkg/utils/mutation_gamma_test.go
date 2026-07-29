package utils

import (
	"math"
	"testing"

	"drift/pkg/core"
)

// gammaModel returns a model whose default class draws effects from the gamma DFE
// (dfe_model=gamma), with every mutation non-neutral (f_neutral=0) and deleterious
// (f_beneficial=0) so the realized effects are exactly -Gamma(shape, mean/shape).
func gammaModel(mu, gshape, gmean float64) *core.Model {
	return &core.Model{
		Parameters: map[string]float64{
			"mu":                mu,
			"f_neutral":         0,
			"f_beneficial":      0,
			"shape":             1,
			"scale":             0.05,
			"Weibull_adj":       1,
			"mu_scale_factor":   1000000,
			"fitness_dominance": 0.5,
			"dfe_gamma_shape":   gshape,
			"dfe_gamma_mean":    gmean,
		},
		StringParams:   map[string]string{"dfe_model": "gamma"},
		FreeParameters: map[string]int{"mutID": 0, "genome_bits": 100000},
	}
}

// TestGammaDFERealizedMean: with dfe_model=gamma the realized |effect| of the
// segregating pool has mean ~ dfe_gamma_mean and shape ~ dfe_gamma_shape (checked
// via the coefficient of variation, CV = 1/sqrt(shape) for a gamma).
func TestGammaDFERealizedMean(t *testing.T) {
	const gshape, gmean = 0.2, 0.01
	SeedRNG(31337)
	m := gammaModel(3, gshape, gmean)
	p := newMutPop()
	for ind := 1; ind <= 4000; ind++ {
		GenerateNewMutations(m, p, ind)
	}
	var sum, sumSq float64
	n := 0
	for _, mut := range p.MutationPool {
		if mut.Effect == 0 {
			t.Fatal("f_neutral=0 should leave no neutral mutations")
		}
		a := -mut.Effect // all deleterious (f_beneficial=0)
		if a <= 0 {
			t.Fatalf("expected deleterious effect, got %v", mut.Effect)
		}
		sum += a
		sumSq += a * a
		n++
	}
	if n < 1000 {
		t.Fatalf("too few mutations to characterize: %d", n)
	}
	mean := sum / float64(n)
	cv := math.Sqrt(sumSq/float64(n)-mean*mean) / mean
	if math.Abs(mean-gmean) > 0.06*gmean {
		t.Errorf("realized mean |effect|=%.5g want ~%.5g", mean, gmean)
	}
	// Gamma CV = 1/sqrt(shape); shape 0.2 => CV ~2.24. Wide band (heavy tail noise).
	wantCV := 1.0 / math.Sqrt(gshape)
	if math.Abs(cv-wantCV) > 0.25*wantCV {
		t.Errorf("realized CV=%.3f want ~%.3f (shape %.2g)", cv, wantCV, gshape)
	}
}

// TestGammaDFEDistinguishableFromWeibull: gamma and weibull, given the same
// nominal mean-effect target, produce detectably different DFE shapes — the gamma
// (shape 0.2) has a far heavier tail (larger CV) than the exponential Weibull
// (shape 1, CV=1). Confirms the new distribution is not a relabelled Weibull.
func TestGammaDFEDistinguishableFromWeibull(t *testing.T) {
	meanEffect := 0.01
	cvOf := func(useGamma bool) float64 {
		SeedRNG(2024)
		var m *core.Model
		if useGamma {
			m = gammaModel(3, 0.2, meanEffect)
		} else {
			// Weibull shape 1 (exponential), scale = meanEffect, adj 1 => mean meanEffect.
			m = &core.Model{
				Parameters: map[string]float64{
					"mu": 3, "f_neutral": 0, "f_beneficial": 0,
					"shape": 1, "scale": meanEffect, "Weibull_adj": 1,
					"mu_scale_factor": 1000000, "fitness_dominance": 0.5,
				},
				StringParams:   map[string]string{},
				FreeParameters: map[string]int{"mutID": 0, "genome_bits": 100000},
			}
		}
		p := newMutPop()
		for ind := 1; ind <= 3000; ind++ {
			GenerateNewMutations(m, p, ind)
		}
		var sum, sumSq float64
		n := 0
		for _, mut := range p.MutationPool {
			a := math.Abs(mut.Effect)
			sum += a
			sumSq += a * a
			n++
		}
		mean := sum / float64(n)
		return math.Sqrt(sumSq/float64(n)-mean*mean) / mean
	}
	gammaCV := cvOf(true)
	weibullCV := cvOf(false)
	if gammaCV <= weibullCV*1.3 {
		t.Errorf("gamma CV=%.3f should be >> weibull CV=%.3f (heavier tail)", gammaCV, weibullCV)
	}
}

// TestGammaDFEParse: a per-class dfe=gamma override with gshape/gmean tokens is
// parsed into the class fields; a class without the token inherits the weibull
// default.
func TestGammaDFEParse(t *testing.T) {
	base := MutationClass{DFEModel: "weibull", GammaShape: 0.2, GammaMean: 0.01, Shape: 1, Scale: 0.05, WeibullAdj: 1}
	c := parseMutationClass("cnv rate=0.5 dfe=gamma gshape=0.3 gmean=0.02", base, 100000)
	if c.DFEModel != "gamma" || c.GammaShape != 0.3 || c.GammaMean != 0.02 {
		t.Errorf("parsed class = %+v, want dfe=gamma gshape=0.3 gmean=0.02", c)
	}
	c2 := parseMutationClass("point rate=8", base, 100000)
	if c2.DFEModel != "weibull" {
		t.Errorf("class without dfe token = %q, want inherited weibull", c2.DFEModel)
	}
}
