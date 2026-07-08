package utils

import (
	"math"
	"testing"

	"drift/pkg/core"
)

// legacyModel returns a model with the pre-spectrum global mutation params only
// (no mutation_classes), used both as the default-spectrum input and to drive the
// legacy-kernel oracle in the byte-identity test.
func legacyModel(mu float64) *core.Model {
	return &core.Model{
		Parameters: map[string]float64{
			"mu":                mu,
			"f_neutral":         0.3,
			"f_beneficial":      0.1,
			"shape":             1,
			"scale":             0.05,
			"Weibull_adj":       1,
			"mu_scale_factor":   1000000,
			"fitness_dominance": 0.5,
		},
		FreeParameters: map[string]int{"mutID": 0, "genome_bits": 100000},
	}
}

func newMutPop() *core.Pop {
	return &core.Pop{
		IndMutations: map[int]map[int][]int{},
		MutationPool: map[int]core.Mutation{},
		MutationHist: map[int]int{},
	}
}

// legacyGenerateNewMutations is a verbatim copy of the pre-spectrum single-class
// kernel. It is the oracle for TestDefaultSpectrumByteIdentical: the refactored
// GenerateNewMutations, given no mutation_classes, must consume the identical RNG
// stream and produce the identical mutation pool.
func legacyGenerateNewMutations(model *core.Model, pop *core.Pop, ind int) {
	numNewMutations := RandPoisson(model.Parameters["mu"])
	dominancePct := 50
	if h, ok := model.Parameters["fitness_dominance"]; ok {
		dominancePct = int(math.Round(h * 100))
	}
	for i := 0; i < numNewMutations; i++ {
		model.FreeParameters["mutID"]++
		mutationID := model.FreeParameters["mutID"]
		position := RandIntn(int(model.FreeParameters["genome_bits"]))
		mutationEffect := 0.0
		isMutationNonNeutral := RandFloat64()
		if isMutationNonNeutral >= model.Parameters["f_neutral"] {
			mutationEffect = WeibullRandom(model.Parameters["shape"], model.Parameters["scale"]) / model.Parameters["Weibull_adj"]
			isMutationDeleterious := RandFloat64()
			if isMutationDeleterious > model.Parameters["f_beneficial"] {
				mutationEffect = -mutationEffect
			}
		}
		pop.MutationHist[int(mutationEffect*model.Parameters["mu_scale_factor"])]++
		strand := RandIntn(2)
		if pop.IndMutations[ind] == nil {
			pop.IndMutations[ind] = make(map[int][]int)
		}
		if pop.IndMutations[ind][strand] == nil {
			pop.IndMutations[ind][strand] = []int{}
		}
		pop.IndMutations[ind][strand] = append(pop.IndMutations[ind][strand], mutationID)
		pop.MutationPool[mutationID] = core.Mutation{
			Id: mutationID, Position: position, Effect: mutationEffect,
			Origin: ind, Count: 1, Dominance: dominancePct,
		}
	}
}

// TestDefaultSpectrumByteIdentical is the compatibility guarantee for the §6h
// neutral harness: with mutation_classes unset, GenerateNewMutations must draw
// the exact same RNG sequence and build the exact same pool as the legacy kernel,
// leaving the shared generator in the identical state afterward.
func TestDefaultSpectrumByteIdentical(t *testing.T) {
	const mu = 8

	// New code path.
	SeedRNG(9271)
	newM := legacyModel(mu)
	newP := newMutPop()
	for ind := 1; ind <= 200; ind++ {
		GenerateNewMutations(newM, newP, ind)
	}
	newTail := []float64{RandFloat64(), RandFloat64(), RandFloat64()}

	// Legacy oracle.
	SeedRNG(9271)
	oldM := legacyModel(mu)
	oldP := newMutPop()
	for ind := 1; ind <= 200; ind++ {
		legacyGenerateNewMutations(oldM, oldP, ind)
	}
	oldTail := []float64{RandFloat64(), RandFloat64(), RandFloat64()}

	if len(newP.MutationPool) != len(oldP.MutationPool) {
		t.Fatalf("pool size: new=%d old=%d", len(newP.MutationPool), len(oldP.MutationPool))
	}
	if len(newP.MutationPool) == 0 {
		t.Fatal("expected some mutations generated")
	}
	for id, om := range oldP.MutationPool {
		nm, ok := newP.MutationPool[id]
		if !ok {
			t.Fatalf("mutation %d missing from new pool", id)
		}
		// Every legacy field must match exactly (Class/Size are new bookkeeping).
		if nm.Id != om.Id || nm.Position != om.Position || nm.Effect != om.Effect ||
			nm.Origin != om.Origin || nm.Count != om.Count || nm.Dominance != om.Dominance {
			t.Fatalf("mutation %d differs:\n new=%+v\n old=%+v", id, nm, om)
		}
		// The default class is point (index 0, size 1).
		if nm.Class != 0 || nm.Size != 1 {
			t.Fatalf("mutation %d default class/size = %d/%d, want 0/1", id, nm.Class, nm.Size)
		}
	}
	// The RNG streams must be aligned after identical consumption.
	for i := range newTail {
		if newTail[i] != oldTail[i] {
			t.Fatalf("RNG stream diverged at tail draw %d: new=%v old=%v", i, newTail[i], oldTail[i])
		}
	}
}

// TestParseMutationClasses checks the spec parser: default single class, field
// inheritance from the model globals, per-field overrides, and that malformed
// tokens/entries are skipped rather than fatal.
func TestParseMutationClasses(t *testing.T) {
	model := legacyModel(10)

	// Empty spec -> single point class carrying the model globals.
	got := MutationClasses(model)
	if len(got) != 1 || got[0].Name != "point" {
		t.Fatalf("empty spec: got %d classes %+v", len(got), got)
	}
	if got[0].Rate != 10 || got[0].Scale != 0.05 || got[0].Dominance != 50 || got[0].Size != 1 {
		t.Fatalf("default class did not inherit globals: %+v", got[0])
	}

	model.StringParams = map[string]string{
		"mutation_classes": "point rate=8; indel rate=1 scale=0.2 size=10 dom=0.2; large_del rate=0.1 fneutral=0 fbeneficial=0 scale=0.5 adj=2 size=1000 bogus=1",
	}
	got = MutationClasses(model)
	if len(got) != 3 {
		t.Fatalf("expected 3 classes, got %d: %+v", len(got), got)
	}

	// point: only rate overridden; everything else inherits globals.
	if got[0].Name != "point" || got[0].Rate != 8 || got[0].Scale != 0.05 ||
		got[0].FNeutral != 0.3 || got[0].FBeneficial != 0.1 || got[0].Size != 1 {
		t.Errorf("point class wrong: %+v", got[0])
	}
	// indel: rate/scale/size/dom overridden, shape/fneutral inherited.
	if got[1].Name != "indel" || got[1].Rate != 1 || got[1].Scale != 0.2 ||
		got[1].Size != 10 || got[1].Dominance != 20 || got[1].Shape != 1 || got[1].FNeutral != 0.3 {
		t.Errorf("indel class wrong: %+v", got[1])
	}
	// large_del: all effect fields overridden; unknown key 'bogus' ignored.
	if got[2].Name != "large_del" || got[2].Rate != 0.1 || got[2].FNeutral != 0 ||
		got[2].FBeneficial != 0 || got[2].Scale != 0.5 || got[2].WeibullAdj != 2 || got[2].Size != 1000 {
		t.Errorf("large_del class wrong: %+v", got[2])
	}

	// All-separator spec falls back to the single default class (no crash).
	model.StringParams["mutation_classes"] = " ; ;; "
	got = MutationClasses(model)
	if len(got) != 1 || got[0].Name != "point" {
		t.Fatalf("all-separator spec should fall back to default: %+v", got)
	}
}

// TestMutationSpectrumDistinguishable drives the real kernel across many
// individuals with a two-class spec and confirms the classes are distinguishable
// end-to-end: their realized counts track their configured rates, their mean
// effect magnitudes track their configured Weibull scales, and each mutation is
// stamped with its class index and size for downstream filtering.
func TestMutationSpectrumDistinguishable(t *testing.T) {
	SeedRNG(4242)
	model := legacyModel(0) // mu unused; rates come from the spec
	// small: rate 5, tiny effects (scale 0.01, size 1)
	// large: rate 1, big effects (scale 0.5, size 1000); all non-neutral+deleterious
	model.StringParams = map[string]string{
		"mutation_classes": "small rate=5 fneutral=0 fbeneficial=0 shape=1 scale=0.01 size=1; large rate=1 fneutral=0 fbeneficial=0 shape=1 scale=0.5 size=1000",
	}

	pop := newMutPop()
	const N = 3000
	for ind := 1; ind <= N; ind++ {
		GenerateNewMutations(model, pop, ind)
	}

	counts := map[int]int{}
	sumAbs := map[int]float64{}
	sizes := map[int]int{}
	for _, m := range pop.MutationPool {
		counts[m.Class]++
		sumAbs[m.Class] += math.Abs(m.Effect)
		sizes[m.Class] = m.Size
	}
	if counts[0] == 0 || counts[1] == 0 {
		t.Fatalf("both classes should be present: %v", counts)
	}
	// Rates 5:1 -> realized counts should be near 5x (loose bounds for MC noise).
	ratio := float64(counts[0]) / float64(counts[1])
	if ratio < 4 || ratio > 6 {
		t.Errorf("class count ratio = %.2f (counts %v), want ~5", ratio, counts)
	}
	// Sizes recorded per class.
	if sizes[0] != 1 || sizes[1] != 1000 {
		t.Errorf("recorded sizes = %v, want {0:1, 1:1000}", sizes)
	}
	// Mean |effect|: class 1 (scale 0.5) should dwarf class 0 (scale 0.01).
	mean0 := sumAbs[0] / float64(counts[0])
	mean1 := sumAbs[1] / float64(counts[1])
	if !(mean1 > 10*mean0) {
		t.Errorf("mean |effect| class0=%.4g class1=%.4g; class1 should be >>class0", mean0, mean1)
	}
}
