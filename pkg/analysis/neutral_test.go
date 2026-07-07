package analysis

import (
	"drift/pkg/core"
	"math"
	"os"
	"testing"
)

// buildMutPop constructs a pop whose mutation pool gives KNOWN, hand-computed
// neutral statistics. Four diploid individuals (ids 1..4, numChrom = 8) carry:
//
//	mut 101: 1 derived copy  (ind1 strand0)                      -> singleton, het
//	mut 102: 2 derived copies(ind1 strand1 + ind2 strand0)       -> doubleton, 2 het
//	mut 103: 4 derived copies(ind3 both strands + ind4 both)     -> 0 het (homozygous)
//	mut 104: 8 derived copies(everyone, both strands)            -> FIXED (excluded)
//	mut 105: 0 derived copies                                    -> absent (excluded)
func buildMutPop() (*core.Pop, []int) {
	// mut 104 is placed on all 8 strands (fixed); mut 105 on none (absent).
	im := map[int]map[int][]int{
		1: {0: {101, 104}, 1: {102, 104}},
		2: {0: {102, 104}, 1: {104}},
		3: {0: {103, 104}, 1: {103, 104}},
		4: {0: {103, 104}, 1: {103, 104}},
	}
	pop := &core.Pop{
		IndData:      map[int][]int{1: {}, 2: {}, 3: {}, 4: {}},
		IndMutations: im,
		MutationPool: map[int]core.Mutation{
			101: {Id: 101, Position: 10},
			102: {Id: 102, Position: 20},
			103: {Id: 103, Position: 30},
			104: {Id: 104, Position: 40},
			105: {Id: 105, Position: 50},
		},
	}
	return pop, []int{1, 2, 3, 4}
}

func approx(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s = %.8f, want %.8f (tol %.1g)", name, got, want, tol)
	}
}

func TestComputeNeutralStats_HandComputed(t *testing.T) {
	pop, ids := buildMutPop()
	s := ComputeNeutralStats(pop, ids)

	if s.NumSamples != 4 || s.NumChrom != 8 {
		t.Fatalf("NumSamples/NumChrom = %d/%d, want 4/8", s.NumSamples, s.NumChrom)
	}
	if s.SegSites != 3 {
		t.Fatalf("SegSites = %d, want 3 (mut104 fixed and mut105 absent are excluded)", s.SegSites)
	}
	if s.Singletons != 1 {
		t.Errorf("Singletons = %d, want 1", s.Singletons)
	}
	if s.Doubletons != 1 {
		t.Errorf("Doubletons = %d, want 1", s.Doubletons)
	}

	// Unfolded SFS: derived counts 1, 2, 4.
	wantUnfolded := map[int]int{1: 1, 2: 1, 4: 1}
	for i, c := range s.UnfoldedSFS {
		if wantUnfolded[i] != c {
			t.Errorf("UnfoldedSFS[%d] = %d, want %d", i, c, wantUnfolded[i])
		}
	}

	// theta_pi = 0.25 + 0.428571 + 0.571429 = 1.25 exactly.
	approx(t, "ThetaPi", s.ThetaPi, 1.25, 1e-9)

	// theta_W = S / a_n, a_n = sum_{i=1}^{7} 1/i.
	a1 := 0.0
	for i := 1; i < 8; i++ {
		a1 += 1.0 / float64(i)
	}
	approx(t, "ThetaW", s.ThetaW, 3.0/a1, 1e-9)

	// Tajima's D, hand-computed above ~ 0.33040.
	approx(t, "TajimasD", s.TajimasD, 0.33040, 1e-3)

	// HWE: meanHe = 0.364583..., meanHo = 0.25, F_IS = 0.314286.
	approx(t, "MeanHe", s.MeanHe, (0.21875+0.375+0.5)/3, 1e-9)
	approx(t, "MeanHo", s.MeanHo, 0.25, 1e-9)
	approx(t, "FIS", s.FIS, 0.3142857, 1e-6)
}

func TestExpectedUnfoldedSFS_Shape(t *testing.T) {
	n := 20
	exp := ExpectedUnfoldedSFS(n)
	// Sums to 1 over classes 1..n-1.
	var sum float64
	for i := 1; i < n; i++ {
		sum += exp[i]
	}
	approx(t, "SFS proportion sum", sum, 1.0, 1e-12)
	// Proportional to 1/i: exp[i] * i should be constant.
	c1 := exp[1] * 1
	for i := 2; i < n; i++ {
		approx(t, "1/i proportionality", exp[i]*float64(i), c1, 1e-12)
	}
	// Boundary classes are zero.
	if exp[0] != 0 || exp[n] != 0 {
		t.Errorf("boundary classes nonzero: exp[0]=%g exp[n]=%g", exp[0], exp[n])
	}
}

func TestSFSShapeChiSquare_PerfectFit(t *testing.T) {
	// An SFS drawn EXACTLY in 1/i proportions should give a chi2 of ~0.
	numChrom := 40
	segSites := 100000
	exp := ExpectedUnfoldedSFS(numChrom)
	unfolded := make([]int, numChrom+1)
	for i := 1; i < numChrom; i++ {
		unfolded[i] = int(math.Round(exp[i] * float64(segSites)))
	}
	chi2, dof := SFSShapeChiSquare(unfolded, segSites, numChrom)
	if dof <= 0 {
		t.Fatalf("dof = %d, want > 0", dof)
	}
	if chi2/float64(dof) > 0.05 {
		t.Errorf("chi2/dof = %.4f for a perfect 1/i fit, want ~0", chi2/float64(dof))
	}
}

func TestSFSShapeChiSquare_Distorted(t *testing.T) {
	// A star-like (all-singleton) SFS is a gross departure from 1/i and must
	// produce a large chi2/dof.
	numChrom := 40
	segSites := 10000
	unfolded := make([]int, numChrom+1)
	unfolded[1] = segSites // everything is a singleton
	chi2, dof := SFSShapeChiSquare(unfolded, segSites, numChrom)
	if dof <= 0 || chi2/float64(dof) < 5 {
		t.Errorf("chi2/dof = %.4f for an all-singleton SFS, want large (>5)", chi2/float64(dof))
	}
}

func TestBuildReport_PassFail(t *testing.T) {
	pop, ids := buildMutPop()
	s := ComputeNeutralStats(pop, ids)

	// Wide tolerances centered on the hand pop's own values so scoring logic (not
	// the characterized defaults) is what's exercised here.
	wide := NeutralTolerances{
		TajimasDCenter: 0, TajimasD: 1.0,
		ThetaRatioCenter: 1.0, ThetaRatio: 0.5,
		SFSShape: 100, FIS: 0.5, FuLiD: 5,
	}
	r := BuildReport(s, 1.0, 4, wide)
	if !r.AllPass() {
		for _, c := range r.Checks {
			if !c.Pass {
				t.Errorf("unexpected FAIL: %s obs=%.4f exp=%.4f tol=%.4f", c.Name, c.Observed, c.Expected, c.Tolerance)
			}
		}
	}
	// ImpliedNe = thetaW / (2*mu).
	approx(t, "ImpliedNe", r.ImpliedNe, s.ThetaW/2.0, 1e-9)

	// Zero tolerances: everything that is not exactly on expectation must FAIL,
	// so the report is not vacuously passing.
	strict := NeutralTolerances{}
	r2 := BuildReport(s, 1.0, 4, strict)
	if r2.AllPass() {
		t.Errorf("report passed under zero tolerances; scoring is not discriminating")
	}
}

func TestComputeNeutralStats_TooFewSamples(t *testing.T) {
	pop, _ := buildMutPop()
	s := ComputeNeutralStats(pop, []int{1})
	if s.SegSites != 0 {
		t.Errorf("SegSites = %d for a single-individual sample, want 0", s.SegSites)
	}
}

func TestSaveNeutralValidation(t *testing.T) {
	pop, _ := buildMutPop()
	model := &core.Model{
		Parameters: map[string]float64{
			"track_mutations":        1,
			"validation_sample_size": 0,
			"mu":                     1.0,
		},
		FreeParameters: map[string]int{"run": 1, "year": 100},
		ModelName:      "Neutral",
		ResultsDir:     t.TempDir(),
	}
	if err := SaveNeutralValidation(model, pop); err != nil {
		t.Fatalf("SaveNeutralValidation: %v", err)
	}
	// A report file must have been written for this run/year.
	path := model.ResultsDir + "/Neutral_validation_run1_year100.csv"
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("expected report at %s: %v", path, err)
	}
	if !containsSub(string(data), "TajimasD") || !containsSub(string(data), "ImpliedNe") {
		t.Errorf("report missing expected rows:\n%s", data)
	}

	// track_mutations off => no-op (nil, no file).
	model.Parameters["track_mutations"] = 0
	model.FreeParameters["year"] = 200
	if err := SaveNeutralValidation(model, pop); err != nil {
		t.Errorf("no-op path returned error: %v", err)
	}
	if _, err := os.Stat(model.ResultsDir + "/Neutral_validation_run1_year200.csv"); err == nil {
		t.Errorf("report written despite track_mutations=0")
	}
}

func containsSub(h, n string) bool {
	for i := 0; i+len(n) <= len(h); i++ {
		if h[i:i+len(n)] == n {
			return true
		}
	}
	return false
}

func TestSampleLiving(t *testing.T) {
	pop := &core.Pop{IndData: map[int][]int{5: {}, 2: {}, 9: {}, 1: {}}}
	all := SampleLiving(pop, 0)
	if len(all) != 4 || all[0] != 1 || all[3] != 9 {
		t.Errorf("SampleLiving(0) = %v, want sorted all ids", all)
	}
	sub := SampleLiving(pop, 2)
	if len(sub) != 2 {
		t.Errorf("SampleLiving(2) returned %d ids, want 2", len(sub))
	}
	// Output is sorted even when subsampling.
	if sub[0] > sub[1] {
		t.Errorf("SampleLiving subsample not sorted: %v", sub)
	}
}
