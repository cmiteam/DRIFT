package utils

import (
	"math"
	"testing"

	"drift/pkg/core"
)

// TestParseRateMap checks the rate-map spec parser: empty/degenerate specs return
// nil (uniform no-op), regions build segments, boundaries clamp to the genome, and
// both ';' and ',' work as region separators.
func TestParseRateMap(t *testing.T) {
	// Empty / whitespace / no-region specs -> nil (fall back to uniform).
	for _, s := range []string{"", "   ", "none", "uniform", ";;", "1000-2000"} {
		if rm := parseRateMap(s, 100000); rm != nil {
			t.Errorf("parseRateMap(%q) = %+v, want nil", s, rm)
		}
	}
	// A map that is baseline everywhere (mult 1) is indistinguishable from uniform.
	if rm := parseRateMap("1000-2000:1", 100000); rm != nil {
		t.Errorf("all-baseline map should be nil, got %+v", rm)
	}
	// Degenerate genome.
	if rm := parseRateMap("0-10:5", 0); rm != nil {
		t.Errorf("genomeBits<=0 should be nil, got %+v", rm)
	}

	// A single hotspot: two baseline flanks + the hotspot = three segments.
	rm := parseRateMap("1000-2000:10", 100000)
	if rm == nil {
		t.Fatal("hotspot map should be non-nil")
	}
	if len(rm.segStart) != 3 {
		t.Fatalf("segments = %d, want 3 (flank/hotspot/flank): %+v", len(rm.segStart), rm.segStart)
	}
	// total weight = 99000 baseline (1000 left flank + 98000 right) + 1000*10
	// hotspot = 109000.
	if rm.total != 109000 {
		t.Errorf("total weight = %v, want 109000", rm.total)
	}

	// Comma separator (per-class inline form) and clamping past genome end.
	rm = parseRateMap("0-100:2,90000-999999:3", 100000)
	if rm == nil {
		t.Fatal("comma-separated map should parse")
	}
	// Last segment must be clamped to the genome length (end 100000, not 999999).
	last := len(rm.segStart) - 1
	if rm.segStart[last]+rm.segLen[last] != 100000 {
		t.Errorf("last segment end = %d, want 100000 (clamped)", rm.segStart[last]+rm.segLen[last])
	}

	// A zero-multiplier "dead" region contributes no weight and no drawable segment,
	// but the surrounding baseline still does.
	rm = parseRateMap("1000-2000:0", 10000)
	if rm == nil {
		t.Fatal("coldspot(0) with baseline flanks should be non-nil")
	}
	if rm.total != 9000 { // 10000 - 1000 dead bits (baseline flanks 1000 + 8000)
		t.Errorf("total weight with dead region = %v, want 9000", rm.total)
	}
}

// TestRateMapDrawConcentrates confirms the weighted draw actually concentrates
// positions in the high-multiplier region roughly in proportion to its weight.
func TestRateMapDrawConcentrates(t *testing.T) {
	SeedRNG(31337)
	const genomeBits = 100000
	// 1% of the genome (bits 1000..2000) at 100x rate. Its share of total weight is
	// (1000*100) / (99000*1 + 1000*100) = 100000/199000 ≈ 0.5025.
	rm := parseRateMap("1000-2000:100", genomeBits)
	if rm == nil {
		t.Fatal("map should be non-nil")
	}
	const draws = 200000
	inHot := 0
	for i := 0; i < draws; i++ {
		p := rm.Draw()
		if p < 0 || p >= genomeBits {
			t.Fatalf("draw %d out of range: %d", i, p)
		}
		if p >= 1000 && p < 2000 {
			inHot++
		}
	}
	frac := float64(inHot) / draws
	if frac < 0.47 || frac > 0.53 {
		t.Errorf("hotspot fraction = %.3f, want ~0.502 (10x concentration in 1%% of genome)", frac)
	}
}

// rateMapModel is a minimal model for exercising the spectrum + rate map together.
func rateMapModel() *core.Model {
	return &core.Model{
		Parameters: map[string]float64{
			"mu": 20, "f_neutral": 1, "f_beneficial": 0,
			"shape": 1, "scale": 0.05, "Weibull_adj": 1,
			"mu_scale_factor": 1000000, "fitness_dominance": 0.5,
		},
		FreeParameters: map[string]int{"mutID": 0, "genome_bits": 100000},
		StringParams:   map[string]string{},
	}
}

// TestRateMapByteIdenticalWhenEmpty verifies the strict no-op: with no
// mutation_rate_map set, GenerateNewMutations draws the identical position stream
// (and leaves the RNG in the identical state) as when the rate-map code is bypassed
// entirely — i.e. the default class carries a nil RateMap and uses RandIntn.
func TestRateMapByteIdenticalWhenEmpty(t *testing.T) {
	classes := MutationClasses(rateMapModel())
	if len(classes) != 1 || classes[0].RateMap != nil {
		t.Fatalf("empty mutation_rate_map must give a single class with nil RateMap, got %+v", classes)
	}

	// Full run comparison: default (nil map) vs an explicit uniform RandIntn oracle.
	SeedRNG(555)
	m := rateMapModel()
	p := newMutPop()
	for ind := 1; ind <= 100; ind++ {
		GenerateNewMutations(m, p, ind)
	}
	newTail := RandFloat64()

	SeedRNG(555)
	oldM := rateMapModel()
	oldP := newMutPop()
	for ind := 1; ind <= 100; ind++ {
		legacyGenerateNewMutations(oldM, oldP, ind)
	}
	oldTail := RandFloat64()

	if len(p.MutationPool) != len(oldP.MutationPool) {
		t.Fatalf("pool size differs: new=%d old=%d", len(p.MutationPool), len(oldP.MutationPool))
	}
	for id, om := range oldP.MutationPool {
		nm := p.MutationPool[id]
		if nm.Position != om.Position {
			t.Fatalf("position differs for mutation %d: new=%d old=%d", id, nm.Position, om.Position)
		}
	}
	if newTail != oldTail {
		t.Fatalf("RNG stream diverged: new=%v old=%v", newTail, oldTail)
	}
}

// TestRateMapEndToEndConcentration drives the real kernel with a hotspot map and
// confirms de-novo mutations concentrate in the hotspot region.
func TestRateMapEndToEndConcentration(t *testing.T) {
	SeedRNG(24680)
	m := rateMapModel()
	// 2% of the genome at 50x. Expected hotspot share ≈ (2000*50)/(98000+2000*50)
	// = 100000/198000 ≈ 0.505.
	m.StringParams["mutation_rate_map"] = "10000-12000:50"

	classes := MutationClasses(m)
	if len(classes) != 1 || classes[0].RateMap == nil {
		t.Fatalf("expected single class with a non-nil RateMap, got %+v", classes)
	}

	p := newMutPop()
	for ind := 1; ind <= 4000; ind++ {
		GenerateNewMutations(m, p, ind)
	}
	if len(p.MutationPool) == 0 {
		t.Fatal("expected mutations generated")
	}
	inHot := 0
	for _, mut := range p.MutationPool {
		if mut.Position >= 10000 && mut.Position < 12000 {
			inHot++
		}
	}
	frac := float64(inHot) / float64(len(p.MutationPool))
	// Loose bounds for Monte-Carlo noise around the ~0.505 expectation; the point is
	// that ~half the mutations land in 2% of the genome (vs ~0.02 under uniform).
	if frac < 0.42 || frac > 0.58 {
		t.Errorf("hotspot fraction = %.3f (%d/%d), want ~0.5", frac, inHot, len(p.MutationPool))
	}
}

// TestPerClassRateMapOverride confirms a class can override the global map (and
// another class can opt back out to uniform) via the ratemap= token.
func TestPerClassRateMapOverride(t *testing.T) {
	m := rateMapModel()
	m.StringParams["mutation_rate_map"] = "0-1000:5" // global map
	// classA: inherits the global map; classB: own hotspot; classC: opts out.
	m.StringParams["mutation_classes"] = "classA rate=1; classB rate=1 ratemap=5000-6000:9,7000-8000:3; classC rate=1 ratemap=none"

	classes := MutationClasses(m)
	if len(classes) != 3 {
		t.Fatalf("expected 3 classes, got %d", len(classes))
	}
	// classA inherits the global map (non-nil, 3 segments over 0-1000 hotspot).
	if classes[0].RateMap == nil {
		t.Errorf("classA should inherit the global rate map")
	}
	// classB has its own two-region map (5 segments: flank/hot/gap/hot/flank).
	if classes[1].RateMap == nil {
		t.Fatalf("classB should have its own rate map")
	}
	if classes[1].RateMap == classes[0].RateMap {
		t.Errorf("classB map should differ from the global (classA) map")
	}
	if len(classes[1].RateMap.segStart) != 5 {
		t.Errorf("classB segments = %d, want 5", len(classes[1].RateMap.segStart))
	}
	// classC opted out -> uniform (nil).
	if classes[2].RateMap != nil {
		t.Errorf("classC ratemap=none should be nil (uniform), got %+v", classes[2].RateMap)
	}
}

// sanity: the weighted expectation math used in the concentration tests.
func TestRateMapWeightSanity(t *testing.T) {
	rm := parseRateMap("1000-2000:100", 100000)
	share := (1000.0 * 100.0) / rm.total
	if math.Abs(share-100000.0/199000.0) > 1e-9 {
		t.Errorf("hotspot weight share = %v, want %v", share, 100000.0/199000.0)
	}
}
