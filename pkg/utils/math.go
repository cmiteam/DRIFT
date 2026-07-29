package utils

import (
	"math"
)

func Min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func Max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func Abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func WeibullRandom(shape, scale float64) float64 {
	u := RandFloat64()
	return scale * math.Pow(-math.Log(u), 1/shape)
}

// NormalRandom returns a standard-normal (mean 0, sd 1) draw from the shared
// seeded generator. Thin wrapper over the source's ziggurat sampler; used by
// GammaRandom for the opt-in gamma DFE (roadmap §6f).
func NormalRandom() float64 { return simRNG.NormFloat64() }

// GammaRandom draws from a Gamma(shape, scale) distribution (mean = shape*scale,
// variance = shape*scale^2) from the shared seeded generator, using the
// Marsaglia & Tsang (2000) method with the shape<1 boost. It powers the opt-in
// literature-anchored gamma DFE (dfe_model=gamma, roadmap §6f): the human
// deleterious distribution of fitness effects is well fit by a gamma with shape
// alpha ~ 0.2 (Eyre-Walker 2007; Kim/Huber/Lohmueller 2017), which a single
// Weibull cannot match. Reproducible under a fixed seed.
//
// This is the ONLY genetics-kernel path §6f touches; it is reached only when a
// mutation class selects DFEModel=="gamma", so the default (Weibull) draw is
// byte-identical and strictly-neutral runs are unaffected.
func GammaRandom(shape, scale float64) float64 {
	if shape <= 0 || scale <= 0 {
		return 0
	}
	if shape < 1 {
		// Boost: X ~ Gamma(a) can be sampled as Gamma(a+1) * U^(1/a).
		u := RandFloat64()
		for u <= 0 {
			u = RandFloat64()
		}
		return GammaRandom(shape+1, scale) * math.Pow(u, 1/shape)
	}
	d := shape - 1.0/3.0
	c := 1.0 / math.Sqrt(9*d)
	for {
		x := NormalRandom()
		v := 1 + c*x
		if v <= 0 {
			continue
		}
		v = v * v * v
		u := RandFloat64()
		if u < 1-0.0331*x*x*x*x {
			return d * v * scale
		}
		if math.Log(u) < 0.5*x*x+d*(1-v+math.Log(v)) {
			return d * v * scale
		}
	}
}
