package analysis

import (
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// W7 parser tests. The hand-built trees below have TMR-K-A values that can be read
// off by eye, which is the point: the walker on the inferred side has to agree with
// the walker on the truth side (pkg/analysis/tmrka.go) about what "K lineages" means,
// and the only way to know that is to check both against arithmetic done by hand.

func TestParseNewickHeightsAndTMRKA(t *testing.T) {
	// A ladder: (A,B) meet at 1, C joins at 2, D joins at 3.
	const nwk = "(((A:1,B:1):1,C:2):1,D:3);"
	root, leaves, err := parseNewick(nwk)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if leaves != 4 {
		t.Fatalf("leaves = %d, want 4", leaves)
	}
	if math.Abs(root.height-3) > 1e-9 {
		t.Errorf("root height = %v, want 3", root.height)
	}

	for _, tc := range []struct {
		k     int
		want  float64
		exact bool
	}{
		{k: 3, want: 1, exact: true},
		{k: 2, want: 2, exact: true},
		{k: 1, want: 3, exact: true},
	} {
		got, exact, ok := tmrKA(root, leaves, tc.k)
		if !ok {
			t.Errorf("K=%d: no result", tc.k)
			continue
		}
		if math.Abs(got-tc.want) > 1e-9 || exact != tc.exact {
			t.Errorf("K=%d: got %v exact=%v, want %v exact=%v", tc.k, got, exact, tc.want, tc.exact)
		}
	}
}

// The case that matters most for this study, and the one the DRIFT-side walker also
// has: ARGweaver's discretised time grid makes simultaneous coalescences common, so
// the lineage count can fall PAST K in a single step. TMR-K-A is still the first time
// the count is K or below, but "K lineages existed" is then FALSE and must be
// reported as such rather than averaged into the coalesced pile.
func TestTMRKASkipsPastK(t *testing.T) {
	// Three lineages meet simultaneously at height 1: the count goes 4 -> 2,
	// never resting at 3.
	root, leaves, err := parseNewick("((A:1,B:1,C:1):1,D:2);")
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	got, exact, ok := tmrKA(root, leaves, 3)
	if !ok {
		t.Fatal("K=3: no result")
	}
	if math.Abs(got-1) > 1e-9 {
		t.Errorf("K=3: got %v, want 1", got)
	}
	if exact {
		t.Error("K=3: reported an exact hit, but the count skipped 3 entirely")
	}

	// Same tree, K=2: the count lands on 2 exactly at the same height.
	if got, exact, _ := tmrKA(root, leaves, 2); math.Abs(got-1) > 1e-9 || !exact {
		t.Errorf("K=2: got %v exact=%v, want 1 exact=true", got, exact)
	}
}

// Real arg-summarize output carries NHX annotations describing the transition to the
// NEXT local tree. They are not part of this tree's topology and must not perturb it.
func TestParseNewickStripsNHX(t *testing.T) {
	plain := "((A:1,B:1):1,C:2);"
	nhx := "((A:1,B:1):1[&&NHX:coal_time=7574.2],C:2)[&&NHX:recomb_time=4975.8];"
	a, la, err := parseNewick(plain)
	if err != nil {
		t.Fatalf("plain: %v", err)
	}
	b, lb, err := parseNewick(nhx)
	if err != nil {
		t.Fatalf("nhx: %v", err)
	}
	if la != lb || math.Abs(a.height-b.height) > 1e-9 {
		t.Errorf("NHX changed the tree: leaves %d vs %d, height %v vs %v", la, lb, a.height, b.height)
	}
}

// Zero-length branches are pervasive in real output (many coalescences land on the
// same grid point). They must not produce a parent shallower than its child.
func TestParseNewickZeroBranches(t *testing.T) {
	root, leaves, err := parseNewick("((A:0.0,(B:0.0,C:0.0):0.0):85.9,D:85.9);")
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if leaves != 4 {
		t.Fatalf("leaves = %d, want 4", leaves)
	}
	if math.Abs(root.height-85.9) > 1e-9 {
		t.Errorf("root height = %v, want 85.9", root.height)
	}
	// Three lineages collapse at height 0, so at K=3 the count has already fallen
	// to 2 before any positive time: TMR-3-A is 0 and is not an exact hit.
	got, exact, ok := tmrKA(root, leaves, 3)
	if !ok || math.Abs(got) > 1e-9 || exact {
		t.Errorf("K=3: got %v exact=%v ok=%v, want 0 exact=false ok=true", got, exact, ok)
	}
}

// The ultrametry check done by hand against real ARGweaver output: the two root
// paths in the sample tree agree at 26422.7 generations.
func TestParseNewickRealShapeUltrametric(t *testing.T) {
	const nwk = "((((A:23.0,B:23.0):2607.4,(C:23.0,D:23.0):2607.4):1397.7," +
		"((E:0.0,F:0.0):1705.9,(G:0.0,H:0.0):1705.9):2322.2):22394.6[&&NHX:coal_time=7574.2]," +
		"(I:51.2,(J:23.0,K:23.0):28.2):26371.5[&&NHX:recomb_time=4975.8]);"
	root, leaves, err := parseNewick(nwk)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if leaves != 11 {
		t.Fatalf("leaves = %d, want 11", leaves)
	}
	if math.Abs(root.height-26422.7) > 1e-6 {
		t.Errorf("root height = %v, want 26422.7", root.height)
	}
}

// End to end: a tiny --tree file plus a DRIFT truth CSV, checking the bp mapping,
// the locus-centric join, and that a locus outside the sampled region is reported as
// uncovered rather than silently dropped.
func TestParseARGweaverTreesJoinsTruth(t *testing.T) {
	dir := t.TempDir()

	// Two MCMC samples over one interval covering bp 1..1000, i.e. bits 0..9 at
	// 100 bp/bit. Locus bit 0 -> bp 1 is inside; locus bit 500 -> bp 50001 is not.
	trees := strings.Join([]string{
		"## arg-summarize 0.3",
		"##chrom\tchromStart\tchromEnd\tMCMC_sample\ttree",
		"1\t0\t1000\t100\t(((A:10,B:10):10,C:20):10,D:30);",
		"1\t0\t1000\t110\t(((A:20,B:20):20,C:40):20,D:60);",
	}, "\n") + "\n"
	treesPath := filepath.Join(dir, "t.tree.txt")
	if err := os.WriteFile(treesPath, []byte(trees), 0644); err != nil {
		t.Fatal(err)
	}

	truth := "locus,k,outcome,tmrka_year,tmrka_gens_ago,tmrca_year,sample_year," +
		"initial_lineages,lineages_at_founding,coalescences,founder_strands_surviving," +
		"created_alleles_surviving\n" +
		"0,3,coalesced,-100,10.00,,500,4000,4,3996,4,4\n" +
		"500,3,censored_at_founding,,,,500,4000,12,3988,12,10\n"
	truthPath := filepath.Join(dir, "truth.csv")
	if err := os.WriteFile(truthPath, []byte(truth), 0644); err != nil {
		t.Fatal(err)
	}

	res, err := ParseARGweaverTrees(ARGweaverParseOptions{
		TreesPath: treesPath, TruthPath: truthPath, K: 3, BpPerBit: 100, GenTime: 25,
	})
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if res.TreesParsed != 2 {
		t.Errorf("TreesParsed = %d, want 2", res.TreesParsed)
	}
	if res.NumLeaves != 4 {
		t.Errorf("NumLeaves = %d, want 4", res.NumLeaves)
	}

	byBit := map[int]ARGweaverTreeStats{}
	for _, l := range res.Loci {
		byBit[l.LocusBit] = l
	}

	in, ok := byBit[0]
	if !ok {
		t.Fatal("locus bit 0 missing")
	}
	if in.Samples != 2 {
		t.Errorf("locus 0: Samples = %d, want 2", in.Samples)
	}
	// TMR-3-A is the first coalescence: 10 in one sample, 20 in the other -> mean 15.
	if math.Abs(in.TMRKAMean-15) > 1e-9 {
		t.Errorf("locus 0: TMRKAMean = %v, want 15", in.TMRKAMean)
	}
	// Root heights 30 and 60 -> mean 45.
	if math.Abs(in.TMRCAMean-45) > 1e-9 {
		t.Errorf("locus 0: TMRCAMean = %v, want 45", in.TMRCAMean)
	}
	if !in.HasTruth || in.TruthOutcome != OutcomeCoalesced {
		t.Errorf("locus 0: truth not joined: %+v", in)
	}
	if in.TruthGens != 10 {
		t.Errorf("locus 0: TruthGens = %v, want 10", in.TruthGens)
	}

	out, ok := byBit[500]
	if !ok {
		t.Fatal("locus bit 500 missing")
	}
	if out.Samples != 0 {
		t.Errorf("locus 500 is outside the sampled region; Samples = %d, want 0", out.Samples)
	}
	if out.TruthOutcome != OutcomeCensored {
		t.Errorf("locus 500: TruthOutcome = %q, want %q", out.TruthOutcome, OutcomeCensored)
	}

	// Span-weighted mean over both trees, both spanning 1000 bp: (10+20)/2 = 15.
	if math.Abs(res.SpanWeighted-15) > 1e-9 {
		t.Errorf("SpanWeighted = %v, want 15", res.SpanWeighted)
	}

	// The comparison CSV must be writable and carry a header.
	outPath := filepath.Join(dir, "cmp.csv")
	if err := WriteARGweaverComparison(outPath, res); err != nil {
		t.Fatalf("write: %v", err)
	}
	b, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(b), "inferred_tmr3a_mean_gen") {
		t.Errorf("comparison CSV header missing the K-specific column:\n%s", string(b))
	}
}

// burnInFixture writes four MCMC samples over one interval. Samples 100 and 110
// are "unconverged" and shallow (TMR-3-A 10 and 20); samples 900 and 910 are
// "burned in" and deep (100 and 200). Pooling all four gives a mean of 82.5,
// which is the error the filter exists to prevent.
func burnInFixture(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	trees := strings.Join([]string{
		"## arg-summarize 0.3",
		"##chrom\tchromStart\tchromEnd\tMCMC_sample\ttree",
		"1\t0\t1000\t100\t(((A:10,B:10):10,C:20):10,D:30);",
		"1\t0\t1000\t110\t(((A:20,B:20):20,C:40):20,D:60);",
		"1\t0\t1000\t900\t(((A:100,B:100):100,C:200):100,D:300);",
		"1\t0\t1000\t910\t(((A:200,B:200):200,C:400):200,D:600);",
	}, "\n") + "\n"
	p := filepath.Join(dir, "t.tree.txt")
	if err := os.WriteFile(p, []byte(trees), 0644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestParseARGweaverMinSampleDropsBurnIn(t *testing.T) {
	p := burnInFixture(t)
	base := ARGweaverParseOptions{TreesPath: p, K: 3, BpPerBit: 100, GenTime: 25, LocusBits: []int{0}}

	// Unfiltered: all four samples pooled.
	all, err := ParseARGweaverTrees(base)
	if err != nil {
		t.Fatal(err)
	}
	if all.TreesParsed != 4 {
		t.Fatalf("unfiltered TreesParsed = %d, want 4", all.TreesParsed)
	}
	if got := all.Loci[0].TMRKAMean; math.Abs(got-82.5) > 1e-9 {
		t.Errorf("unfiltered mean = %v, want 82.5 (10+20+100+200)/4", got)
	}
	if all.SampleLo != 100 || all.SampleHi != 910 {
		t.Errorf("unfiltered sample range = %d-%d, want 100-910", all.SampleLo, all.SampleHi)
	}
	if all.RowsFiltered != 0 {
		t.Errorf("unfiltered RowsFiltered = %d, want 0", all.RowsFiltered)
	}

	// Burn-in filter: keep only the converged tail.
	opts := base
	opts.MinSample = 900
	burned, err := ParseARGweaverTrees(opts)
	if err != nil {
		t.Fatal(err)
	}
	if burned.TreesParsed != 2 {
		t.Fatalf("filtered TreesParsed = %d, want 2", burned.TreesParsed)
	}
	if burned.RowsFiltered != 2 {
		t.Errorf("RowsFiltered = %d, want 2", burned.RowsFiltered)
	}
	if got := burned.Loci[0].TMRKAMean; math.Abs(got-150) > 1e-9 {
		t.Errorf("filtered mean = %v, want 150 (100+200)/2", got)
	}
	if burned.SampleLo != 900 || burned.SampleHi != 910 {
		t.Errorf("filtered sample range = %d-%d, want 900-910", burned.SampleLo, burned.SampleHi)
	}
	// The span-weighted genome-wide figure must respect the filter too, otherwise
	// the cross-check against arg-summarize --tmrca would silently use burn-in.
	if math.Abs(burned.SpanWeighted-150) > 1e-9 {
		t.Errorf("filtered SpanWeighted = %v, want 150", burned.SpanWeighted)
	}
}

// The early/late split: capping and flooring the same chain must give the two
// halves separately, which is how you demonstrate the answer has stopped moving.
func TestParseARGweaverMaxSampleEarlyLateSplit(t *testing.T) {
	p := burnInFixture(t)
	base := ARGweaverParseOptions{TreesPath: p, K: 3, BpPerBit: 100, GenTime: 25, LocusBits: []int{0}}

	early := base
	early.MaxSample = 110
	e, err := ParseARGweaverTrees(early)
	if err != nil {
		t.Fatal(err)
	}
	if got := e.Loci[0].TMRKAMean; math.Abs(got-15) > 1e-9 {
		t.Errorf("early mean = %v, want 15", got)
	}

	late := base
	late.MinSample = 900
	l, err := ParseARGweaverTrees(late)
	if err != nil {
		t.Fatal(err)
	}
	if got := l.Loci[0].TMRKAMean; math.Abs(got-150) > 1e-9 {
		t.Errorf("late mean = %v, want 150", got)
	}

	// A window that selects the middle of the chain.
	mid := base
	mid.MinSample, mid.MaxSample = 110, 900
	m, err := ParseARGweaverTrees(mid)
	if err != nil {
		t.Fatal(err)
	}
	if m.TreesParsed != 2 {
		t.Errorf("windowed TreesParsed = %d, want 2 (samples 110 and 900)", m.TreesParsed)
	}
	if got := m.Loci[0].TMRKAMean; math.Abs(got-60) > 1e-9 {
		t.Errorf("windowed mean = %v, want 60 (20+100)/2", got)
	}
}

// A filter that matches nothing must produce an empty result that says so, not a
// clean-looking zero. This is the failure mode that would quietly turn a typo in
// the cut-off into a published number.
func TestParseARGweaverFilterMatchingNothingIsVisible(t *testing.T) {
	p := burnInFixture(t)
	res, err := ParseARGweaverTrees(ARGweaverParseOptions{
		TreesPath: p, K: 3, BpPerBit: 100, GenTime: 25, LocusBits: []int{0},
		MinSample: 99999,
	})
	if err != nil {
		t.Fatal(err)
	}
	if res.TreesParsed != 0 {
		t.Errorf("TreesParsed = %d, want 0", res.TreesParsed)
	}
	if res.RowsFiltered != 4 {
		t.Errorf("RowsFiltered = %d, want 4", res.RowsFiltered)
	}
	if res.RowsRead-res.RowsFiltered != 0 {
		t.Errorf("retained = %d, want 0", res.RowsRead-res.RowsFiltered)
	}
	if len(res.Loci) > 0 && res.Loci[0].Samples != 0 {
		t.Errorf("locus Samples = %d, want 0", res.Loci[0].Samples)
	}
}
