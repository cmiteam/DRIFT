package analysis

import (
	"drift/pkg/core"
	"math"
	"strings"
	"testing"
)

// The scaling rationale (TMR4A.md W7, §6 risk 1). Every number here is hand-computed from
// the two per-generation anchors, because the whole claim of the module is that theta/rho
// is derivable rather than assumed — a test that only re-ran the same arithmetic would
// prove nothing.

// scalingModel: one chromosome, two 100-bit arms, 50 cM each (200 bits, 100 cM total),
// mu mutations per individual per generation, obligate crossover reserved at 50 cM.
func scalingModel(mu, armCM, obligateCM float64) *core.Model {
	return &core.Model{
		ChromosomeArms: map[int]map[int][]int{
			1: {0: {0, 100}, 1: {100, 100}},
		},
		ArmCM: map[int]map[int]float64{1: {0: armCM, 1: armCM}},
		Parameters: map[string]float64{
			"mu":                 mu,
			"recomb_obligate_cM": obligateCM,
			"recomb_cM_per_bit":  1,
		},
		StringParams: map[string]string{
			"recombination_model": "map",
			"linkage_model":       "arm",
		},
		FreeParameters: map[string]int{"genome_bits": 200, "run": 1, "year": 0},
		ModelName:      "ScalingTest",
	}
}

// The central identity, hand-computed: mu is per INDIVIDUAL but crossovers are per
// GAMETE, so theta/rho is (mu/2) / crossovers. Getting that factor of 2 wrong is the
// easiest way to be exactly twice as wrong as you think.
func TestRateAnchorsThetaOverRho(t *testing.T) {
	// mu = 6 per individual => 3 per gamete.
	// 100 cM total, 50 cM obligate reservation => 1 obligate + (100-50)/100 = 1.5
	// crossovers per meiosis. theta/rho = 3 / 1.5 = 2.
	model := scalingModel(6, 50, 50)
	a := ComputeRateAnchors(model, 100, buildContigs(model))

	if a.MutationsPerIndividual != 6 {
		t.Errorf("mutations per individual = %g, want 6", a.MutationsPerIndividual)
	}
	if a.MutationsPerGamete != 3 {
		t.Errorf("mutations per gamete = %g, want 3 (each mutation lands on ONE strand)",
			a.MutationsPerGamete)
	}
	if a.CrossoversPerMeiosis != 1.5 {
		t.Errorf("crossovers per meiosis = %g, want 1.5 (1 obligate + 0.5 Poisson extras)",
			a.CrossoversPerMeiosis)
	}
	if a.CrossoversNominal != 1.0 {
		t.Errorf("nominal crossovers = %g, want 1.0 (100 cM / 100)", a.CrossoversNominal)
	}
	if a.ThetaOverRho != 2 {
		t.Errorf("theta/rho = %g, want 2", a.ThetaOverRho)
	}
}

// The claim W7 rests on: bp-per-bit buys resolution and NOTHING else. Scaling it changes
// both per-bp rates by the same factor and leaves theta/rho untouched.
func TestRateAnchorsBpPerBitIsScaleInvariant(t *testing.T) {
	model := scalingModel(6, 50, 50)
	contigs := buildContigs(model)

	coarse := ComputeRateAnchors(model, 10, contigs)
	fine := ComputeRateAnchors(model, 1000, contigs)

	if coarse.ThetaOverRho != fine.ThetaOverRho {
		t.Errorf("theta/rho moved with bp-per-bit: %g at 10 vs %g at 1000 — it must be "+
			"scale-invariant, that is the whole reason X is free", coarse.ThetaOverRho, fine.ThetaOverRho)
	}
	// Both per-bp rates fall by exactly the 100x change in B.
	if math.Abs(coarse.MuPerBp/fine.MuPerBp-100) > 1e-9 {
		t.Errorf("mu per bp ratio = %g, want 100", coarse.MuPerBp/fine.MuPerBp)
	}
	if len(coarse.Regions) != 1 || len(fine.Regions) != 1 {
		t.Fatalf("expected one region each")
	}
	if math.Abs(coarse.Regions[0].RhoPerBp/fine.Regions[0].RhoPerBp-100) > 1e-9 {
		t.Errorf("rho per bp ratio = %g, want 100",
			coarse.Regions[0].RhoPerBp/fine.Regions[0].RhoPerBp)
	}
	// And the declared region length is genome_bits x B.
	if fine.RegionBpTotal != 200*1000 {
		t.Errorf("region bp = %d, want 200000", fine.RegionBpTotal)
	}
}

// The finding this report exists to surface: on a short cM map the OBLIGATE crossover
// dominates and the realized recombination rate is nothing like the map's. Telling
// ARGweaver a rho derived from the cM column would be wrong by that whole factor.
func TestRateAnchorsObligateCrossoverInflationIsReported(t *testing.T) {
	// 5 cM total against a 50 cM obligate reservation: the map implies 0.05 crossovers,
	// the engine performs 1. A 20x inflation.
	model := scalingModel(6, 2.5, 50)
	a := ComputeRateAnchors(model, 100, buildContigs(model))

	if a.CrossoversNominal != 0.05 {
		t.Fatalf("nominal crossovers = %g, want 0.05", a.CrossoversNominal)
	}
	if a.CrossoversPerMeiosis != 1 {
		t.Fatalf("realized crossovers = %g, want 1 (the obligate crossover alone)",
			a.CrossoversPerMeiosis)
	}
	if !hasWarning(a, "obligate crossover dominates") {
		t.Errorf("a 20x gap between realized and nominal recombination must be warned about; "+
			"warnings were %v", a.Warnings)
	}
	// The exported rho uses the REALIZED count — the rate the engine actually produced.
	if a.Regions[0].RhoPerBp != 1/float64(200*100) {
		t.Errorf("rho per bp = %g, want the realized 1/(200 bits x 100 bp)", a.Regions[0].RhoPerBp)
	}
}

// Being far off the human ratio is legitimate (W9 sweeps exactly this kind of thing) but
// must never happen silently: the mainstream comparison arm depends on it.
func TestRateAnchorsWarnsWhenFarFromHumanRatio(t *testing.T) {
	// mu = 0.02 per individual => 0.01 per gamete against 1.5 crossovers: theta/rho
	// 0.00667, about 0.0064x human.
	off := ComputeRateAnchors(scalingModel(0.02, 50, 50), 100, buildContigs(scalingModel(0.02, 50, 50)))
	if !hasWarning(off, "theta/rho") {
		t.Errorf("theta/rho 150x below human must be warned about; got %v", off.Warnings)
	}

	// On target: 1.5 crossovers per meiosis needs mu/2 = 1.5 x 1.0417, so mu = 3.125.
	onTarget := scalingModel(3.125, 50, 50)
	a := ComputeRateAnchors(onTarget, 100, buildContigs(onTarget))
	if math.Abs(a.RatioVsHuman-1) > 1e-6 {
		t.Fatalf("fixture is not actually on target: ratio vs human = %g", a.RatioVsHuman)
	}
	if hasWarning(a, "theta/rho") {
		t.Errorf("an on-target ratio must not warn; got %v", a.Warnings)
	}
}

// The two engine settings an export is invalid without, each warned about by name.
func TestRateAnchorsWarnsOnWrongEngineSettings(t *testing.T) {
	model := scalingModel(6, 50, 50)
	model.StringParams["recombination_model"] = "legacy"
	model.StringParams["linkage_model"] = "legacy"
	a := ComputeRateAnchors(model, 100, buildContigs(model))

	if !hasWarning(a, "recombination_model") {
		t.Errorf("legacy recombination must be warned about (TMR4A.md §8a); got %v", a.Warnings)
	}
	if !hasWarning(a, "linkage_model") {
		t.Errorf("legacy linkage must be warned about (TMR4A.md §8c); got %v", a.Warnings)
	}
}

// Coarse resolution defeats the purpose: a local tree spans hundreds of bp to ~1 kb, so
// breakpoints inside a bit are unresolvable once the bit is wider than that.
func TestRateAnchorsWarnsOnCoarseResolution(t *testing.T) {
	model := scalingModel(6, 50, 50)
	if a := ComputeRateAnchors(model, 5000, buildContigs(model)); !hasWarning(a, "bp per bit") {
		t.Errorf("5000 bp per bit must be warned about; got %v", a.Warnings)
	}
	if a := ComputeRateAnchors(model, 100, buildContigs(model)); hasWarning(a, "bp per bit") {
		t.Errorf("100 bp per bit is the target and must not warn; got %v", a.Warnings)
	}
}

// Arms that recombine at different rates are invisible to ARGweaver, which is told one
// rho for the whole region.
func TestRateAnchorsWarnsOnWithinRegionRateHeterogeneity(t *testing.T) {
	model := scalingModel(6, 50, 50)
	model.ArmCM = map[int]map[int]float64{1: {0: 10, 1: 90}} // 9x apart
	a := ComputeRateAnchors(model, 100, buildContigs(model))
	if !hasWarning(a, "differ in recombination rate") {
		t.Errorf("a 9x arm-rate spread must be warned about; got %v", a.Warnings)
	}
}

// The rationale renders the three things §6 risk 1 says must be defended in writing.
func TestRationaleStatesTheThreeRequiredNumbers(t *testing.T) {
	model := scalingModel(6, 50, 50)
	a := ComputeRateAnchors(model, 100, buildContigs(model))

	var sb strings.Builder
	a.WriteRationale(&sb, model)
	out := sb.String()

	for _, want := range []string{
		"1 genome bit = 100 bp",                   // the bp-per-bit convention
		"mutations per individual per generation", // anchor 1
		"crossovers per meiosis, realized",        // anchor 2
		"theta/rho",                               // the scale-invariant quantity
		"--mutrate",                               // what arg-sample must be told
		"--recombrate",
	} {
		if !strings.Contains(out, want) {
			t.Errorf("rationale is missing %q:\n%s", want, out)
		}
	}
}

func hasWarning(a *RateAnchors, substr string) bool {
	for _, w := range a.Warnings {
		if strings.Contains(w, substr) {
			return true
		}
	}
	return false
}
