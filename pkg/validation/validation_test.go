package validation

import (
	"os"
	"path/filepath"
	"testing"
)

// TestNeutralValidation_EndToEnd drives the REAL engine through a strictly neutral
// scenario and asserts the §6h credibility-backbone checks against DRIFT's
// characterized neutral baseline. This is the end-to-end reproducibility run: if a
// future change perturbs the neutral dynamics (selection leaking in, a broken
// mutation kernel, a demography no-op regressing), these bands catch it.
//
// The two textbook expectations DRIFT meets exactly are asserted tightly (HWE
// F_IS ~ 0; theta_pi/theta_W consistency; stable, reproducible implied Ne). Tajima's
// D is asserted against DRIFT's reproducible rare-variant-excess baseline (~ -0.7,
// NOT the Wright-Fisher 0) — see CharacterizedTolerances for why.
func TestNeutralValidation_EndToEnd(t *testing.T) {
	cfg := DefaultNeutralConfig()
	if testing.Short() {
		cfg.Replicates = 6
		cfg.BurnInYears = 800
	}
	res := RunNeutral(cfg)
	tol := CharacterizedTolerances()

	if res.MeanSegSites < 20 {
		t.Fatalf("only %.0f segregating sites on average; the neutral run produced too little "+
			"variation to validate (raise mu / burn-in)", res.MeanSegSites)
	}

	for _, c := range res.Checks(tol) {
		if !c.Pass {
			t.Errorf("FAIL %s: observed %.4f, expected %.4f +/- %.4f  (%s)",
				c.Name, c.Observed, c.Expected, c.Tolerance, c.Note)
		}
	}

	// Sanity: mean final census should stay near the target N (no drift to
	// extinction or runaway growth over the burn-in).
	if res.MeanFinalN < 0.5*float64(cfg.N) {
		t.Errorf("mean final N = %.0f collapsed well below target N = %d", res.MeanFinalN, cfg.N)
	}

	t.Logf("neutral: R=%d N=%d mu=%.2g burn=%dy | S=%.0f thetaW=%.2f thetaPi=%.2f pi/W=%.3f "+
		"D=%.3f(SEM%.3f) F_IS=%.3f Ne=%.1f(Ne/N %.2f, CV %.2f) SFSchi2/dof=%.2f",
		cfg.Replicates, cfg.N, cfg.Mu, cfg.BurnInYears, res.MeanSegSites, res.MeanThetaW,
		res.MeanThetaPi, res.MeanThetaR, res.MeanTajimasD, res.SEMTajimasD, res.MeanFIS,
		res.MeanImpliedNe, res.MeanImpliedNe/float64(cfg.N), res.CVImpliedNe, res.SFSChi2PerDof)
}

// TestNeutralValidation_Reproducible asserts that the whole end-to-end run is
// bit-reproducible under a fixed base seed — the determinism guarantee the rest of
// DRIFT depends on. Uses a small, fast config.
func TestNeutralValidation_Reproducible(t *testing.T) {
	cfg := DefaultNeutralConfig()
	cfg.Replicates = 3
	cfg.BurnInYears = 400

	a := RunNeutral(cfg)
	b := RunNeutral(cfg)

	if a.MeanTajimasD != b.MeanTajimasD ||
		a.MeanThetaW != b.MeanThetaW ||
		a.MeanThetaPi != b.MeanThetaPi ||
		a.MeanFIS != b.MeanFIS ||
		a.MeanSegSites != b.MeanSegSites {
		t.Errorf("neutral run not bit-reproducible:\n  A: D=%v thetaW=%v thetaPi=%v FIS=%v S=%v\n  B: D=%v thetaW=%v thetaPi=%v FIS=%v S=%v",
			a.MeanTajimasD, a.MeanThetaW, a.MeanThetaPi, a.MeanFIS, a.MeanSegSites,
			b.MeanTajimasD, b.MeanThetaW, b.MeanThetaPi, b.MeanFIS, b.MeanSegSites)
	}

	// Per-replicate stats must match exactly too.
	for i := range a.Replicates {
		ra, rb := a.Replicates[i], b.Replicates[i]
		if ra.FinalN != rb.FinalN || ra.Stats.SegSites != rb.Stats.SegSites ||
			ra.Stats.TajimasD != rb.Stats.TajimasD {
			t.Errorf("replicate %d (seed %d) not reproducible: A(finalN=%d,S=%d,D=%v) B(finalN=%d,S=%d,D=%v)",
				i, ra.Seed, ra.FinalN, ra.Stats.SegSites, ra.Stats.TajimasD,
				rb.FinalN, rb.Stats.SegSites, rb.Stats.TajimasD)
		}
	}
}

// TestWriteAveragedReport checks the report CSV is written with the expected
// structure (context block, pooled SFS, scored checks).
func TestWriteAveragedReport(t *testing.T) {
	cfg := DefaultNeutralConfig()
	cfg.Replicates = 3
	cfg.BurnInYears = 300
	res := RunNeutral(cfg)

	path := filepath.Join(t.TempDir(), "neutral_validation.csv")
	if err := WriteAveragedReport(path, res, CharacterizedTolerances()); err != nil {
		t.Fatalf("WriteAveragedReport: %v", err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read report: %v", err)
	}
	s := string(data)
	for _, want := range []string{"MeanTajimasD", "UnfoldedFreqClass", "TajimasD_characterized", "HWE_FIS"} {
		if !containsLine(s, want) {
			t.Errorf("report missing expected content %q", want)
		}
	}
}

func containsLine(haystack, needle string) bool {
	return len(haystack) > 0 && len(needle) > 0 && (indexOf(haystack, needle) >= 0)
}

func indexOf(h, n string) int {
	for i := 0; i+len(n) <= len(h); i++ {
		if h[i:i+len(n)] == n {
			return i
		}
	}
	return -1
}
