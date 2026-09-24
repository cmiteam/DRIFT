package analysis

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// writeSites builds a minimal `.sites` file from per-haplotype sequences, so the tests
// exercise the real parser rather than a hand-built HaploSample.
func writeSites(t *testing.T, names []string, seqs []string) string {
	t.Helper()
	if len(names) != len(seqs) {
		t.Fatalf("test setup: %d names for %d sequences", len(names), len(seqs))
	}
	var b strings.Builder
	b.WriteString("NAMES")
	for _, n := range names {
		b.WriteString("\t" + n)
	}
	b.WriteString("\nREGION\t1\t1\t100000\n")
	for s := 0; s < len(seqs[0]); s++ {
		// One site per 100 bp, matching the export's bp-per-bit convention.
		fmt.Fprintf(&b, "%d\t", s*100+1)
		for _, seq := range seqs {
			b.WriteByte(seq[s])
		}
		b.WriteByte('\n')
	}
	path := filepath.Join(t.TempDir(), "test.sites")
	if err := os.WriteFile(path, []byte(b.String()), 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

// Two designed classes separated by a wide gulf, with one de-novo difference inside a
// class: exactly the geometry the created-allele model builds.
func twoClassSample(t *testing.T) *HaploSample {
	t.Helper()
	names := []string{"IND1_1", "IND1_2", "IND2_1", "IND2_2"}
	seqs := []string{
		"AAAAAAAAAA",
		"AAAAAAAAAT", // one site from its class-mate
		"TTTTTTTTTT",
		"TTTTTTTTTT",
	}
	h, err := LoadSitesWindow(writeSites(t, names, seqs), 0)
	if err != nil {
		t.Fatal(err)
	}
	return h
}

func TestDistanceIsFractionOfDifferingSites(t *testing.T) {
	h := twoClassSample(t)
	if h.NSites != 10 {
		t.Fatalf("expected 10 sites, got %d", h.NSites)
	}
	if got := h.Distance(0, 1); math.Abs(got-0.1) > 1e-9 {
		t.Errorf("one differing site in ten: want 0.1, got %v", got)
	}
	if got := h.Distance(0, 2); math.Abs(got-1.0) > 1e-9 {
		t.Errorf("every site differs: want 1.0, got %v", got)
	}
	if got := h.Distance(2, 3); got != 0 {
		t.Errorf("identical haplotypes: want 0, got %v", got)
	}
}

// The window is not cosmetic: arg-sample saw only --region, so clustering the whole file
// would measure a different object from the one the MCMC was given.
func TestWindowExcludesSitesBeyondMaxPos(t *testing.T) {
	names := []string{"h1", "h2"}
	seqs := []string{"AAAAAAAAAA", "AAAAATTTTT"}
	path := writeSites(t, names, seqs)

	full, err := LoadSitesWindow(path, 0)
	if err != nil {
		t.Fatal(err)
	}
	if full.NSites != 10 || full.Distance(0, 1) != 0.5 {
		t.Fatalf("whole file: got %d sites, distance %v", full.NSites, full.Distance(0, 1))
	}

	// Sites sit at 1, 101, 201 ... so a 500 bp window keeps the first five — all of
	// which are identical between the two haplotypes.
	win, err := LoadSitesWindow(path, 500)
	if err != nil {
		t.Fatal(err)
	}
	if win.NSites != 5 {
		t.Fatalf("500 bp window: want 5 sites, got %d", win.NSites)
	}
	if d := win.Distance(0, 1); d != 0 {
		t.Errorf("the differing sites are outside the window: want 0, got %v", d)
	}
}

func TestIndividualIndexPairsStrands(t *testing.T) {
	h := twoClassSample(t)
	names, cols := h.IndividualIndex()
	if len(names) != 2 {
		t.Fatalf("want 2 individuals, got %d: %v", len(names), names)
	}
	// Both strands of a diploid must travel together — SampleIDs resamples
	// individuals, not haplotypes, and splitting them would understate the loss of
	// rare classes.
	for _, n := range names {
		if len(cols[n]) != 2 {
			t.Errorf("individual %s has %d columns, want 2", n, len(cols[n]))
		}
	}
}

func TestSingleLinkageRecoversDesignedClasses(t *testing.T) {
	h := twoClassSample(t)
	all := []int{0, 1, 2, 3}

	// A cut above the within-class difference (0.1) and below the between-class
	// gulf (1.0) must find exactly the two designed classes.
	c := h.Cluster(all, 0.5)
	if len(c.Classes) != 2 {
		t.Fatalf("want 2 classes at t=0.5, got %d (%v)", len(c.Classes), c.Sizes())
	}
	if got := c.Sizes(); got[0] != 2 || got[1] != 2 {
		t.Errorf("want sizes [2 2], got %v", got)
	}

	// Below the within-class difference the de-novo mutation splits its own class off.
	if c := h.Cluster(all, 0.05); len(c.Classes) != 3 {
		t.Errorf("want 3 classes at t=0.05, got %d (%v)", len(c.Classes), c.Sizes())
	}
	// Above the gulf everything merges.
	if c := h.Cluster(all, 1.0); len(c.Classes) != 1 {
		t.Errorf("want 1 class at t=1.0, got %d", len(c.Classes))
	}
}

// The headline rule commits to no constant, so it has to find the designed gulf on its
// own. Single linkage joins the two classes at their CLOSEST pair, so the merge heights
// are 0 (the identical pair), 0.1 (the de-novo difference) and 0.9 (the nearest
// cross-class pair, not the 1.0 farthest one) — a gulf of 0.8 to be found.
func TestLargestGapCutFindsTheDesignedGulf(t *testing.T) {
	h := twoClassSample(t)
	edges := h.mstEdges([]int{0, 1, 2, 3})
	heights := make([]float64, len(edges))
	for i, e := range edges {
		heights[i] = e.dist
	}
	cut, gap := LargestGapCut(heights, 0.01)
	if cut <= 0.1 || cut >= 0.9 {
		t.Errorf("cut %v is not inside the gulf (0.1, 0.9)", cut)
	}
	if math.Abs(gap-0.8) > 1e-9 {
		t.Errorf("want gap 0.8, got %v", gap)
	}
	if c := h.Cluster([]int{0, 1, 2, 3}, cut); len(c.Classes) != 2 {
		t.Errorf("the parameter-free cut should give 2 classes, got %d", len(c.Classes))
	}
}

// The whole point of the differing-site-count metric is that it is a rescaling, so it
// must not change which haplotypes group together — only the units a threshold is quoted
// in. If this ever fails, the two metrics are answering different questions and the
// n-scaling series is not comparable to the calibration.
func TestScaleChangesUnitsNotClasses(t *testing.T) {
	h := twoClassSample(t)
	all := []int{0, 1, 2, 3}

	fractionSizes := h.Cluster(all, 0.5).Sizes()

	h.SetScale(1) // distances become differing-site counts
	if d := h.Distance(0, 1); math.Abs(d-1) > 1e-9 {
		t.Errorf("one differing site: want distance 1, got %v", d)
	}
	// 0.5 of ten sites is five sites, so the equivalent cut is 5.
	countSizes := h.Cluster(all, 5).Sizes()

	if len(fractionSizes) != len(countSizes) {
		t.Fatalf("metric changed the class count: %v vs %v", fractionSizes, countSizes)
	}
	for i := range fractionSizes {
		if fractionSizes[i] != countSizes[i] {
			t.Fatalf("metric changed the classes: %v vs %v", fractionSizes, countSizes)
		}
	}
	if got := h.ThresholdGrid(); &got[0] != &SiteCountThresholds[0] {
		t.Error("at scale 1 the sweep grid should be the site-count grid")
	}
}

// Without the guard the widest gap can be the one between "identical" and "the first real
// merge", which would report every distinct sequence as its own class.
func TestLargestGapCutIgnoresGapsBelowMinHeight(t *testing.T) {
	// Merges at 0, 0, 0.30, 0.34: the 0->0.30 jump is wider than 0.30->0.34, so
	// without a guard the cut lands below 0.30 and splits the real class.
	heights := []float64{0, 0, 0.30, 0.34}
	if cut, _ := LargestGapCut(heights, 0); cut >= 0.30 {
		t.Fatalf("test premise wrong: unguarded cut was %v", cut)
	}
	cut, _ := LargestGapCut(heights, 0.32)
	if cut <= 0.30 {
		t.Errorf("with minHeight 0.32 the cut should sit above 0.30, got %v", cut)
	}
}

func TestPlateauReturnsWidestAgreeingRun(t *testing.T) {
	sweep := []SweepPoint{
		{Threshold: 0.01, Classes: 5},
		{Threshold: 0.02, Classes: 3},
		{Threshold: 0.10, Classes: 2},
		{Threshold: 0.15, Classes: 2},
		{Threshold: 0.20, Classes: 2},
		{Threshold: 0.25, Classes: 1},
	}
	classes, lo, hi := Plateau(sweep)
	if classes != 2 || lo != 0.10 || hi != 0.20 {
		t.Errorf("want 2 classes over 0.10-0.20, got %d over %v-%v", classes, lo, hi)
	}
}

func TestClusteringSizesAreLargestFirstAndStable(t *testing.T) {
	names := []string{"a_1", "a_2", "b_1", "b_2", "c_1", "c_2"}
	seqs := []string{
		"AAAA", "AAAA", "AAAA", // one class of 3
		"TTTT", "TTTT", // one class of 2
		"ACAC", // a singleton
	}
	h, err := LoadSitesWindow(writeSites(t, names, seqs), 0)
	if err != nil {
		t.Fatal(err)
	}
	c := h.Cluster([]int{0, 1, 2, 3, 4, 5}, 0.1)
	want := []int{3, 2, 1}
	got := c.Sizes()
	if len(got) != len(want) {
		t.Fatalf("want sizes %v, got %v", want, got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("want sizes %v, got %v", want, got)
		}
	}
	// Effective classes: 1/((1/2)^2+(1/3)^2+(1/6)^2) = 2.57...
	if e := c.Simpson(); math.Abs(e-2.5714) > 1e-3 {
		t.Errorf("want effective classes 2.5714, got %v", e)
	}
}

// The point of stage 1: a class carried by one individual is missed by a draw that does
// not include that individual, and the bootstrap has to see that.
func TestBootstrapLosesARareClassAtTheExpectedRate(t *testing.T) {
	// Ten individuals: nine in one class, one carrying a divergent class alone.
	var names []string
	var seqs []string
	for i := 0; i < 10; i++ {
		names = append(names, fmt.Sprintf("IND%d_1", i), fmt.Sprintf("IND%d_2", i))
		if i == 9 {
			seqs = append(seqs, "TTTTTTTTTT", "TTTTTTTTTT")
			continue
		}
		seqs = append(seqs, "AAAAAAAAAA", "AAAAAAAAAA")
	}
	h, err := LoadSitesWindow(writeSites(t, names, seqs), 0)
	if err != nil {
		t.Fatal(err)
	}

	res, err := RunBootstrap(h, "test", 0, BootstrapOptions{
		NIndividuals: 5, Draws: 2000, Seed: 1, Threshold: 0.5, MinHeight: 0.01,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Full.Classes) != 2 {
		t.Fatalf("full population should hold 2 classes, got %d", len(res.Full.Classes))
	}

	// Drawing 5 of 10 individuals includes the rare one exactly half the time, so the
	// class count is 2 about half the draws and 1 the rest.
	hit := res.ClassHitRate[1]
	if math.Abs(hit-0.5) > 0.05 {
		t.Errorf("rare class capture rate %v, want ~0.5", hit)
	}
	two := float64(res.CountHistogram[2]) / float64(len(res.Draws))
	if math.Abs(two-0.5) > 0.05 {
		t.Errorf("draws returning 2 classes: %v, want ~0.5", two)
	}
	if res.CountHistogram[1]+res.CountHistogram[2] != len(res.Draws) {
		t.Errorf("every draw should return 1 or 2 classes, got %v", res.CountHistogram)
	}
}

// A fixed seed must give a fixed answer — the resampling is an input to a study result
// and has to be reproducible from the recorded seed alone.
func TestBootstrapIsDeterministicForASeed(t *testing.T) {
	h := twoClassSample(t)
	opts := BootstrapOptions{NIndividuals: 1, Draws: 50, Seed: 7, Threshold: 0.5}
	a, err := RunBootstrap(h, "test", 0, opts)
	if err != nil {
		t.Fatal(err)
	}
	b, err := RunBootstrap(h, "test", 0, opts)
	if err != nil {
		t.Fatal(err)
	}
	for i := range a.Draws {
		if a.Draws[i].ClassesRecluster != b.Draws[i].ClassesRecluster ||
			a.Draws[i].ClassesRepresented != b.Draws[i].ClassesRepresented {
			t.Fatalf("draw %d differs between two runs at the same seed", i)
		}
	}
}

func TestBootstrapRejectsMoreIndividualsThanExist(t *testing.T) {
	h := twoClassSample(t)
	if _, err := RunBootstrap(h, "test", 0, BootstrapOptions{
		NIndividuals: 99, Draws: 1, Seed: 1, Threshold: 0.5,
	}); err == nil {
		t.Error("asking for more individuals than the file holds should be an error")
	}
}

func TestWriteBootstrapProducesTheFourCSVs(t *testing.T) {
	h := twoClassSample(t)
	res, err := RunBootstrap(h, "test", 0, BootstrapOptions{
		NIndividuals: 2, Draws: 5, Seed: 1, Threshold: 0.5,
	})
	if err != nil {
		t.Fatal(err)
	}
	prefix := filepath.Join(t.TempDir(), "out")
	written, err := WriteBootstrap(prefix, res)
	if err != nil {
		t.Fatal(err)
	}
	if len(written) != 4 {
		t.Fatalf("want 4 CSVs, got %d: %v", len(written), written)
	}
	for _, p := range written {
		info, err := os.Stat(p)
		if err != nil {
			t.Errorf("%s: %v", p, err)
			continue
		}
		if info.Size() == 0 {
			t.Errorf("%s is empty", p)
		}
	}
}
