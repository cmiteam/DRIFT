package validation

import (
	"math"
	"testing"
)

// The TMR-K-A analytic anchor — TMR4A.md W8, fourth bullet. See tmrka_anchor.go for what
// is being anchored against what and why the absolute band is characterized rather than
// textbook.

// AnchorTolerances are the anchor's pass/fail bands.
//
// The two exact ones are exact because they hold for ANY genealogy: TMR-K-A cannot rise
// with K, and a locus the walker has no provenance for must not appear at all.
//
// The two banded ones are characterized to what DRIFT actually does, measured across four
// independent base seeds at the default anchor settings:
//
//	Ne_pool/Ne_geneal   0.951, 0.959, 0.971, 1.047     band [0.75, 1.30]
//	CV of Ne across K   0.040, 0.017, 0.031, 0.051     band <= 0.20
//
// The bands are wide enough to pass reliably on a different seed or a shortened run, and
// narrow enough to catch what they exist to catch: a walk whose absolute timescale is
// wrong by a factor (the ratio moves), or a trajectory whose shape is wrong (the spread
// across K opens up). A systematic half-depth bug moves the ratio to ~2.
func AnchorTolerances() (ratioLo, ratioHi, maxCV, maxCensored float64) {
	return 0.75, 1.30, 0.20, 0.10
}

// TestTMRKAAnchorMatchesCoalescentExpectation is the absolute-timescale check the rest of
// W8 does not make. Every other walker test either fixes the genealogy by hand or checks
// the walk against itself; all of them would still pass if the walk systematically
// reported the wrong depth. This one measures Ne twice from one neutral run — once from
// the mutation pool's segregating sites (theta_W/2mu, the §6h route) and once from the
// walker's own TMR-K-A depths inverted through 4Ne(1/K - 1/n) — and requires the two to
// agree. The two share no code and no substrate.
func TestTMRKAAnchorMatchesCoalescentExpectation(t *testing.T) {
	cfg := DefaultTMRKAAnchorConfig()
	if testing.Short() {
		cfg.Neutral.Replicates = 3
	}
	res := RunTMRKAAnchor(cfg)

	t.Logf("n=%d lineages, realised generation time %.2f y, census N=%d",
		res.Lineages, res.GenerationTime, res.CensusN)
	t.Logf("Ne from mutation pool (theta_W/2mu) = %.2f  (Ne/N = %.3f)",
		res.NeFromPool, res.NeFromPool/float64(res.CensusN))
	t.Logf("Ne from true genealogy            = %.2f  (CV across K = %.3f)",
		res.MeanNeFromGen, res.CVNeFromGen)
	t.Logf("ratio pool/genealogy = %.3f over %d locus-genealogies", res.NeRatio, res.TotalLoci)
	for _, row := range res.ByK {
		t.Logf("  K=%2d: reached %3d, censored %3d, mean depth %7.1f gens (sd %.1f) => Ne %6.2f",
			row.K, row.Reached, row.Censored, row.MeanGens, row.SDGens, row.ImpliedNe)
	}

	if res.TotalLoci == 0 {
		t.Fatal("no locus genealogies were walked; the anchor measured nothing")
	}
	if res.GenerationTime <= 0 {
		t.Fatal("no realised generation time; every depth below would be in arbitrary units")
	}
	if res.MissingARG != 0 {
		t.Errorf("%d loci had no strand provenance; the walk had no substrate at them", res.MissingARG)
	}

	// EXACT. Raising K can only bring the target nearer the present, at every locus,
	// under every model. A rise is a broken trajectory, not a tolerance question.
	if res.OrderViolations != 0 {
		t.Errorf("%d loci reported a deeper TMR-K-A at a LARGER K; the trajectory is not monotone",
			res.OrderViolations)
	}

	ratioLo, ratioHi, maxCV, maxCensored := AnchorTolerances()

	// The anchor only means anything on loci that actually coalesced to K. A run too
	// short for that reports `censored_at_founding`, which is a true answer to a
	// different question (TMR4A.md §6 risk 5) and must not be quietly averaged in.
	if cf := res.CensoredFraction(); cf > maxCensored {
		t.Errorf("%.1f%% of loci were censored at founding before reaching the smallest K; "+
			"the run is too short for a coalescence depth to exist (raise BurnInYears)", 100*cf)
	}

	// Ne-FREE SHAPE CHECK. 4Ne(1/K - 1/n) is the only K-dependence a constant-Ne
	// coalescent admits, so the Ne implied at different K must agree whatever Ne is.
	// Spread here is a direct measurement of departure from Kingman's shape and needs no
	// second estimate to compare against.
	if res.CVNeFromGen > maxCV {
		t.Errorf("implied Ne varies by CV %.3f across K (max %.2f); the lineage-count trajectory "+
			"does not follow 4Ne(1/K - 1/n)", res.CVNeFromGen, maxCV)
	}

	// ABSOLUTE CHECK, in the characterized band. See AnchorTolerances.
	if res.NeRatio < ratioLo || res.NeRatio > ratioHi {
		t.Errorf("Ne from the mutation pool (%.2f) and Ne from the true genealogy (%.2f) disagree: "+
			"ratio %.3f outside the characterized band [%.2f, %.2f]",
			res.NeFromPool, res.MeanNeFromGen, res.NeRatio, ratioLo, ratioHi)
	}

	// Sanity on the bridge between the two: DRIFT's compressed neutral life history
	// (maturity 5, menopause at 0.9*24) cannot produce a generation time outside this.
	if res.GenerationTime < 5 || res.GenerationTime > 22 {
		t.Errorf("realised generation time %.2f y is outside the reproductive window the life "+
			"history allows; the years-to-generations conversion is wrong", res.GenerationTime)
	}
}

// TestTMRKAAnchorDepthDependsOnLocusPlacement pins the reason the anchor spreads its
// focal loci uniformly over genome positions instead of putting one per chromosome.
//
// Under DRIFT's default `recombination_model="legacy"` a meiosis mask is ONE contiguous
// interior segment per chromosome, with endpoints drawn uniformly in the two arms. So the
// probability that a position is inherited from parental copy 0 falls from ~1/2 at the
// centromere to ~0 at a chromosome end: a locus at the start of an arm is taken from the
// same parental copy at nearly every meiosis, its lineage almost never switches strand,
// and its marginal genealogy is correspondingly shallow.
//
// The effect is large — the implied Ne spans more than a factor of two across placements —
// and it is invisible unless someone measures it. Without this test the anchor's locus
// spec looks like an arbitrary detail and the obvious "simplification" to one locus per
// chromosome would break it while looking tidier.
//
// This is a property of the engine's recombination model, not of the walker: the walker
// reports each of these genealogies correctly. It matters for W7, where DRIFT output is
// handed to ARGweaver, whose model assumes recombination is position-homogeneous.
func TestTMRKAAnchorDepthDependsOnLocusPlacement(t *testing.T) {
	if testing.Short() {
		t.Skip("locus-placement contrast runs three neutral scenarios")
	}
	placement := func(loci string) *TMRKAAnchorResult {
		cfg := DefaultTMRKAAnchorConfig()
		cfg.Neutral.ARGLoci = loci
		cfg.Neutral.BurnInYears = 2500 // the contrast is structural; it shows up early
		cfg.Neutral.Replicates = 3
		cfg.Ks = []int{4}
		return RunTMRKAAnchor(cfg)
	}
	// BitsPerChrom is 2048, so chromosome c spans [2048c, 2048c+2048) with its
	// centromere at 2048c+1024.
	starts := placement("0-14336:2048")
	centro := placement("1024-15360:2048")
	uniform := placement("0-16383:997")

	t.Logf("implied Ne by focal-locus placement, same scenario, %d loci each:",
		starts.TotalLoci)
	t.Logf("  chromosome starts %6.2f   uniform %6.2f   centromeres %6.2f   (pool says %6.2f)",
		starts.MeanNeFromGen, uniform.MeanNeFromGen, centro.MeanNeFromGen, uniform.NeFromPool)

	for _, r := range []*TMRKAAnchorResult{starts, centro, uniform} {
		if r.TotalLoci == 0 || r.MeanNeFromGen <= 0 {
			t.Fatal("a placement produced no usable genealogy")
		}
		if r.OrderViolations != 0 {
			t.Errorf("%d order violations; the walk is wrong at this placement regardless of depth",
				r.OrderViolations)
		}
	}

	if !(starts.MeanNeFromGen < uniform.MeanNeFromGen && uniform.MeanNeFromGen < centro.MeanNeFromGen) {
		t.Errorf("expected genealogical depth to increase from chromosome starts to uniform to "+
			"centromeres; got %.2f / %.2f / %.2f",
			starts.MeanNeFromGen, uniform.MeanNeFromGen, centro.MeanNeFromGen)
	}
	if spread := centro.MeanNeFromGen / starts.MeanNeFromGen; spread < 2 {
		t.Errorf("centromere-to-chromosome-start depth ratio is only %.2f; the legacy mask's "+
			"position dependence has changed and the anchor's locus spec needs revisiting", spread)
	}
	// Only the uniform placement is comparable with the mutation pool, which drops its
	// mutations at uniformly-drawn positions and therefore averages over exactly this
	// heterogeneity. That is the whole justification for the default spec.
	ratioLo, ratioHi, _, _ := AnchorTolerances()
	if uniform.NeRatio < ratioLo || uniform.NeRatio > ratioHi {
		t.Errorf("uniform placement should be the one that matches the pool; ratio %.3f is outside [%.2f, %.2f]",
			uniform.NeRatio, ratioLo, ratioHi)
	}
	for _, bad := range []*TMRKAAnchorResult{starts, centro} {
		if bad.NeRatio >= ratioLo && bad.NeRatio <= ratioHi {
			t.Errorf("a single-locus-per-chromosome placement landed inside the anchor band "+
				"(ratio %.3f); the contrast this test documents has gone away", bad.NeRatio)
		}
	}
}

// THE GATE (TMR4A.md §2): every item must leave a default run byte-identical. The
// genealogy capture the anchor switches on draws no RNG of its own — the pedigree path
// never did, and the neutral scenario already tracks mutations so the meiosis masks
// track_arg reads are drawn either way — so a neutral run with capture on must produce
// EXACTLY the statistics it produces with capture off, to the last bit of every float.
func TestTMRKAAnchorIsByteIdentical(t *testing.T) {
	base := DefaultNeutralConfig()
	base.Replicates = 3
	base.BurnInYears = 400

	off := RunNeutral(base)

	on := base
	on.TrackPedigree = true
	on.ARGLoci = "0-16383:997"
	got := RunNeutral(on)

	if off.MeanTajimasD != got.MeanTajimasD || off.MeanThetaW != got.MeanThetaW ||
		off.MeanThetaPi != got.MeanThetaPi || off.MeanFIS != got.MeanFIS ||
		off.MeanSegSites != got.MeanSegSites || off.MeanFinalN != got.MeanFinalN {
		t.Errorf("genealogy capture perturbed the neutral run:\n  off: D=%v thetaW=%v thetaPi=%v FIS=%v S=%v N=%v\n  on:  D=%v thetaW=%v thetaPi=%v FIS=%v S=%v N=%v",
			off.MeanTajimasD, off.MeanThetaW, off.MeanThetaPi, off.MeanFIS, off.MeanSegSites, off.MeanFinalN,
			got.MeanTajimasD, got.MeanThetaW, got.MeanThetaPi, got.MeanFIS, got.MeanSegSites, got.MeanFinalN)
	}
	for i := range off.Replicates {
		a, b := off.Replicates[i], got.Replicates[i]
		if a.FinalN != b.FinalN || a.Stats.SegSites != b.Stats.SegSites ||
			a.Stats.TajimasD != b.Stats.TajimasD || a.Stats.ThetaW != b.Stats.ThetaW {
			t.Errorf("replicate %d (seed %d) differs: off(N=%d,S=%d,D=%v) on(N=%d,S=%d,D=%v)",
				i, a.Seed, a.FinalN, a.Stats.SegSites, a.Stats.TajimasD,
				b.FinalN, b.Stats.SegSites, b.Stats.TajimasD)
		}
	}

	// With capture off nothing is retained at all — no silent allocation on the §6h path.
	for i, rr := range off.Replicates {
		if rr.Pop != nil || rr.Model != nil || rr.SampleIDs != nil {
			t.Errorf("replicate %d retained its population with capture off", i)
		}
	}
	// With it on, the substrate the walker needs is actually there.
	for i, rr := range got.Replicates {
		if rr.Pop == nil || !rr.Pop.PedigreeOn() || len(rr.SampleIDs) == 0 {
			t.Fatalf("replicate %d has no pedigree to walk despite capture being on", i)
		}
	}
}

// The anchor's own aggregation must not depend on Go's randomized map iteration order:
// the same neutral run walked twice has to produce the identical numbers, or every band
// above is measuring noise as well as signal.
func TestTMRKAAnchorIsReproducible(t *testing.T) {
	cfg := DefaultTMRKAAnchorConfig()
	cfg.Neutral.BurnInYears = 800
	cfg.Neutral.Replicates = 2

	a := RunTMRKAAnchor(cfg)
	b := RunTMRKAAnchor(cfg)

	if a.GenerationTime != b.GenerationTime || a.NeFromPool != b.NeFromPool ||
		a.MeanNeFromGen != b.MeanNeFromGen || a.TotalLoci != b.TotalLoci ||
		a.OrderViolations != b.OrderViolations {
		t.Fatalf("anchor not reproducible:\n  A: genTime=%v Ne_pool=%v Ne_gen=%v loci=%d\n  B: genTime=%v Ne_pool=%v Ne_gen=%v loci=%d",
			a.GenerationTime, a.NeFromPool, a.MeanNeFromGen, a.TotalLoci,
			b.GenerationTime, b.NeFromPool, b.MeanNeFromGen, b.TotalLoci)
	}
	for i := range a.ByK {
		if a.ByK[i] != b.ByK[i] {
			t.Errorf("K=%d row differs between identical anchor runs: %+v vs %+v",
				a.ByK[i].K, a.ByK[i], b.ByK[i])
		}
	}
	if math.IsNaN(a.NeRatio) {
		t.Error("Ne ratio is NaN")
	}
}
