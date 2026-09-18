package analysis

import (
	"bytes"
	"compress/gzip"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"testing"
)

// A hand-checkable ultrametric tree on 6 leaves. Coalescences sit at heights
// 10, 20, 30, 40, 50, so the lineage count runs 6→5 at 10, 5→4 at 20, 4→3 at 30,
// 3→2 at 40, 2→1 at 50. Therefore TMR-4-A = 20 and TMRCA = 50.
const smcTree1 = `(((((a:10,b:10):10,c:20):10,d:30):10,e:40):10,f:50);`

// The same shape scaled by two: TMR-4-A = 40, TMRCA = 100.
const smcTree2 = `(((((a:20,b:20):20,c:40):20,d:60):20,e:80):20,f:100);`

// smcFile writes a minimal but format-faithful .smc: NAMES, REGION, then TREE
// lines with 1-based inclusive coordinates, plus an SPR line that must be ignored.
func smcFile(t *testing.T, path string, spans [][2]int, trees []string) {
	t.Helper()
	var b bytes.Buffer
	b.WriteString("NAMES\ta\tb\tc\td\te\tf\n")
	b.WriteString("REGION\tchr1\t1\t1000\n")
	for i, tr := range trees {
		fmtLine := "TREE\t" + strconv.Itoa(spans[i][0]) + "\t" + strconv.Itoa(spans[i][1]) + "\t" + tr + "\n"
		b.WriteString(fmtLine)
		b.WriteString("SPR\t" + strconv.Itoa(spans[i][1]) + "\t3\t12.5\t4\t19.0\n")
	}
	if filepath.Ext(path) == ".gz" {
		var gz bytes.Buffer
		w := gzip.NewWriter(&gz)
		if _, err := w.Write(b.Bytes()); err != nil {
			t.Fatal(err)
		}
		w.Close()
		if err := os.WriteFile(path, gz.Bytes(), 0o644); err != nil {
			t.Fatal(err)
		}
		return
	}
	if err := os.WriteFile(path, b.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
}

func smcCloseTo(t *testing.T, what string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s = %.4f, want %.4f", what, got, want)
	}
}

// The core acceptance test, and the one that matters for the published claim: bp
// weighting must be able to move the median, because the whole 431 kya → 495 kya
// step in Swamidass's analysis is exactly that.
func TestScanSMCWeightingMovesTheMedian(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "chr1.2400.smc")
	// Tree 1 spans 100 bp (TMR-4-A 20); tree 2 spans 900 bp (TMR-4-A 40).
	smcFile(t, p, [][2]int{{1, 100}, {101, 1000}}, []string{smcTree1, smcTree2})

	res, err := ScanSMC(SMCScanOptions{Paths: []string{p}, K: 4, GenTime: 25})
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Samples) != 1 {
		t.Fatalf("samples = %d, want 1", len(res.Samples))
	}
	s := res.Samples[0]
	if s.Iteration != "2400" {
		t.Errorf("iteration = %q, want %q", s.Iteration, "2400")
	}
	if s.Skipped != 0 {
		t.Errorf("skipped = %d, want 0", s.Skipped)
	}
	d := s.TMRKA
	if d.Trees != 2 {
		t.Fatalf("trees = %d, want 2", d.Trees)
	}
	if d.SpanBp != 1000 {
		t.Errorf("span = %d bp, want 1000 (1-based inclusive)", d.SpanBp)
	}
	if d.DistinctVals != 2 {
		t.Errorf("distinct grid values = %d, want 2", d.DistinctVals)
	}

	smcCloseTo(t, "unweighted median", d.Median, 20, 1e-9)
	smcCloseTo(t, "bp-weighted median", d.WMedian, 40, 1e-9)
	smcCloseTo(t, "unweighted mean", d.Mean, 30, 1e-9)
	smcCloseTo(t, "bp-weighted mean", d.WMean, 38, 1e-9)
	smcCloseTo(t, "unweighted stdev", d.Stdev, math.Sqrt(200), 1e-9)
	smcCloseTo(t, "min", d.Min, 20, 1e-9)
	smcCloseTo(t, "max", d.Max, 40, 1e-9)

	// TMRCA is the root height, so the ratio is computed on the same trees.
	smcCloseTo(t, "TMRCA bp-weighted median", s.TMRCA.WMedian, 100, 1e-9)
	smcCloseTo(t, "TMR4A/TMRCA", s.RatioWeighted(), 0.4, 1e-9)
}

// Percentiles must be read off occupied grid points, never interpolated between
// them — the sampler never produced an intermediate value, so reporting one would
// invent precision. With weights 100/900 on values 20/40, everything at or above
// the 25th percentile is 40 and only the 5th percentile is still 20.
func TestScanSMCQuantilesSnapToGrid(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "chr1.2400.smc")
	smcFile(t, p, [][2]int{{1, 100}, {101, 1000}}, []string{smcTree1, smcTree2})

	res, err := ScanSMC(SMCScanOptions{Paths: []string{p}, K: 4})
	if err != nil {
		t.Fatal(err)
	}
	d := res.Samples[0].TMRKA
	for _, c := range []struct {
		name string
		got  float64
		want float64
	}{
		{"WP5", d.WP5, 20},
		{"WQ25", d.WQ25, 40},
		{"WQ75", d.WQ75, 40},
		{"WP95", d.WP95, 40},
	} {
		smcCloseTo(t, c.name, c.got, c.want, 1e-9)
	}
}

// The point of the whole exercise: several MCMC samples must produce a spread of
// genome-wide medians rather than a single number.
func TestScanSMCAcrossMCMCSamples(t *testing.T) {
	dir := t.TempDir()
	// Iteration 100: both trees shallow → weighted median 20.
	smcFile(t, filepath.Join(dir, "chr1.100.smc"),
		[][2]int{{1, 500}, {501, 1000}}, []string{smcTree1, smcTree1})
	// Iteration 200: both trees deep → weighted median 40.
	smcFile(t, filepath.Join(dir, "chr1.200.smc"),
		[][2]int{{1, 500}, {501, 1000}}, []string{smcTree2, smcTree2})
	// Iteration 300, gzipped, mixed → weighted median 40 (900 of 1000 bp deep).
	smcFile(t, filepath.Join(dir, "chr1.300.smc.gz"),
		[][2]int{{1, 100}, {101, 1000}}, []string{smcTree1, smcTree2})

	res, err := ScanSMC(SMCScanOptions{Paths: []string{dir}, K: 4, GenTime: 25})
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Samples) != 3 {
		t.Fatalf("samples = %d, want 3", len(res.Samples))
	}
	// Samples must be ordered numerically by iteration, not lexically.
	for i, want := range []string{"100", "200", "300"} {
		if res.Samples[i].Iteration != want {
			t.Errorf("sample[%d] iteration = %q, want %q", i, res.Samples[i].Iteration, want)
		}
	}
	smcCloseTo(t, "across-sample min median", res.MedianMin, 20, 1e-9)
	smcCloseTo(t, "across-sample max median", res.MedianMax, 40, 1e-9)
	smcCloseTo(t, "across-sample mean median", res.MedianMean, (20+40+40)/3.0, 1e-9)
	if res.MedianStdev <= 0 {
		t.Errorf("across-sample stdev = %.4f, want > 0", res.MedianStdev)
	}
}

// One genome-wide sample is assembled from many per-block files. They must pool
// into a single distribution, not be mistaken for a posterior.
func TestScanSMCPoolsBlocksOfOneSample(t *testing.T) {
	dir := t.TempDir()
	smcFile(t, filepath.Join(dir, "chr1_0-1000.2400.smc"),
		[][2]int{{1, 1000}}, []string{smcTree1})
	smcFile(t, filepath.Join(dir, "chr2_0-1000.2400.smc"),
		[][2]int{{1, 1000}}, []string{smcTree2})

	res, err := ScanSMC(SMCScanOptions{Paths: []string{dir}, K: 4})
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Samples) != 1 {
		t.Fatalf("samples = %d, want 1 (both blocks are iteration 2400)", len(res.Samples))
	}
	s := res.Samples[0]
	if s.Files != 2 {
		t.Errorf("files = %d, want 2", s.Files)
	}
	if s.TMRKA.Trees != 2 || s.TMRKA.SpanBp != 2000 {
		t.Errorf("trees/span = %d/%d, want 2/2000", s.TMRKA.Trees, s.TMRKA.SpanBp)
	}
	// Equal spans, so the weighted median takes the lower of the two grid points.
	smcCloseTo(t, "bp-weighted median", s.TMRKA.WMedian, 20, 1e-9)
}

// An iteration held in a directory rather than the filename must still group.
func TestScanSMCIterationFromPath(t *testing.T) {
	dir := t.TempDir()
	for _, it := range []string{"100", "200"} {
		sub := filepath.Join(dir, it)
		if err := os.MkdirAll(sub, 0o755); err != nil {
			t.Fatal(err)
		}
		tree := smcTree1
		if it == "200" {
			tree = smcTree2
		}
		smcFile(t, filepath.Join(sub, "chr1.smc"), [][2]int{{1, 1000}}, []string{tree})
	}

	// Without the option the files carry no iteration and pool into one sample —
	// which would silently pass off a pooled posterior as a single draw.
	pooled, err := ScanSMC(SMCScanOptions{Paths: []string{dir}, K: 4})
	if err != nil {
		t.Fatal(err)
	}
	if len(pooled.Samples) != 1 {
		t.Fatalf("pooled samples = %d, want 1", len(pooled.Samples))
	}

	split, err := ScanSMC(SMCScanOptions{Paths: []string{dir}, K: 4, IterFromPath: true})
	if err != nil {
		t.Fatal(err)
	}
	if len(split.Samples) != 2 {
		t.Fatalf("split samples = %d, want 2", len(split.Samples))
	}
	smcCloseTo(t, "across-sample min", split.MedianMin, 20, 1e-9)
	smcCloseTo(t, "across-sample max", split.MedianMax, 40, 1e-9)
}

// TMR-K-A must reduce to the root height at K=1, which is what makes the
// TMR4A/TMRCA ratio an internal consistency check rather than two pipelines.
func TestScanSMCKEqualsOneIsTMRCA(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "chr1.2400.smc")
	smcFile(t, p, [][2]int{{1, 1000}}, []string{smcTree1})

	res, err := ScanSMC(SMCScanOptions{Paths: []string{p}, K: 1})
	if err != nil {
		t.Fatal(err)
	}
	s := res.Samples[0]
	smcCloseTo(t, "TMR-1-A", s.TMRKA.WMedian, 50, 1e-9)
	smcCloseTo(t, "TMRCA", s.TMRCA.WMedian, 50, 1e-9)
	smcCloseTo(t, "ratio", s.RatioWeighted(), 1.0, 1e-9)
}

func TestWriteSMCScan(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "chr1.2400.smc")
	smcFile(t, p, [][2]int{{1, 100}, {101, 1000}}, []string{smcTree1, smcTree2})

	res, err := ScanSMC(SMCScanOptions{Paths: []string{p}, K: 4, GenTime: 25})
	if err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(dir, "scan.csv")
	if err := WriteSMCScan(out, res); err != nil {
		t.Fatal(err)
	}
	b, err := os.ReadFile(out)
	if err != nil {
		t.Fatal(err)
	}
	got := string(b)
	for _, want := range []string{"tmr4a_wmedian_gen", "tmr4a_wmedian_yr", "2400,1,2,1000,2,40.0000"} {
		if !bytes.Contains([]byte(got), []byte(want)) {
			t.Errorf("CSV missing %q\ngot:\n%s", want, got)
		}
	}
	// 40 generations at 25 yr/gen must surface as 1000 years.
	if !bytes.Contains(b, []byte("1000.0\n")) {
		t.Errorf("CSV missing year conversion 1000.0\ngot:\n%s", got)
	}
}
