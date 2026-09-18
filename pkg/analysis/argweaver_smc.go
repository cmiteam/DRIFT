package analysis

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

// Raw ARGweaver .smc scanner — the genome-wide TMR-K-A distribution, and its
// spread ACROSS MCMC samples.
//
// WHY THIS EXISTS ALONGSIDE argweaver_parse.go. That file answers "what does
// ARGweaver infer at the twenty focal loci where DRIFT knows the truth" — a
// locus-centric join against ground truth, which is the calibration study. This
// file answers a different question, and it is the one the PUBLISHED claim rests
// on: given a genome-wide set of local trees, what is the distribution of TMR-K-A
// over the whole autosome, and how much does that distribution move between one
// MCMC sample and another? There is no truth to join against here — the input is
// real human data (Rasmussen et al. 2014, 54 Complete Genomics diploids = 108
// haploids), where nobody knows the answer.
//
// WHAT THE PUBLISHED ANALYSIS DID, AND DID NOT, DO. Swamidass's TMR4A ≈ 495 kya
// took ONE MCMC sample — iteration 2400 of one chain — parsed ~12.5M local trees
// covering ~90% of the autosome, took the 4th coalescent in each, and reported the
// bp-span-weighted MEDIAN (495 kya; 431 kya unweighted). Two spreads were never
// reported and neither is in print anywhere:
//
//	(1) across loci, WITHIN that one sample — the CDF was plotted but only the
//	    median was ever quoted; no mean, no quartiles, no range;
//	(2) across MCMC samples — never examined at all. Rasmussen ran ~2,000
//	    iterations per ~2 Mb block retaining every tenth, and the posted archive is
//	    424 GB at ~1.5 GB per genome-wide sample, i.e. of order 280 samples.
//
// The published "± 100 kya" is not either of these: it is a stated judgement about
// grid quantization and mutation-rate constancy, not a computed interval.
//
// So this scanner reports (1) per sample and (2) across samples, and reproducing
// the published 431/495 kya from a single sample is the acceptance test that the
// pipeline is reading the data the same way before any new number is quoted.
//
// THE INPUT is ARGweaver's native .smc(.gz), confirmed against argweaver/smc.py in
// CshlSiepelLab/argweaver: tab-delimited, one record per line,
//
//	NAMES   <name> <name> ...
//	REGION  <chrom> <start> <end>
//	TREE    <start> <end> <newick>
//	SPR     <pos> <recomb_node> <recomb_time> <coal_node> <coal_time>
//
// TREE start/end are 1-based inclusive, so a tree's bp span is end-start+1. Branch
// lengths are in GENERATIONS. SPR lines describe the transition to the next local
// tree and are ignored: each TREE line already carries the full topology for its
// own interval. Reading .smc directly means no smc2bed/tabix/bedops dependency and
// no arg-summarize round trip.
//
// THE STATISTIC is computed by tmrKA() in argweaver_parse.go — the same walker,
// with the same "first time the count is K or below" rule as the DRIFT-side
// pkg/analysis/tmrka.go:177 — so the inferred and true sides cannot drift apart.
//
// QUANTILES ARE EXACT, NOT SAMPLED, and the reason is worth stating. ARGweaver
// discretises time onto a grid of --ntimes points, so TMR-K-A over 12.5M trees
// takes only a few dozen DISTINCT values. Accumulating (value → count, bp) in a
// map therefore gives exact weighted quantiles in kilobytes rather than holding
// 12.5M float64s, and it makes the grid coarseness directly visible: DistinctVals
// is the number of grid points the whole genome-wide answer is being read off.
// That is not incidental — grid quantization is one of the two things the
// published ±100 kya was justified by, so it should be measured, not assumed.

// smcBucket accumulates the trees that landed on one discrete time-grid value.
type smcBucket struct {
	trees int64
	spanB int64
}

// SMCDist is one distribution of local-tree depths, summarised both unweighted
// (one vote per local tree) and weighted by bp span. Swamidass reported the
// weighted median; the unweighted median is the 431 kya figure. Keeping both is
// what makes the 431→495 step auditable rather than a black box.
type SMCDist struct {
	Trees        int64
	SpanBp       int64
	DistinctVals int // grid points actually occupied — the quantization footprint

	// Unweighted, one vote per local tree.
	Mean   float64
	Stdev  float64
	Min    float64
	Max    float64
	P5     float64
	Q25    float64
	Median float64
	Q75    float64
	P95    float64

	// Weighted by bp span. WMedian is the published statistic.
	WMean   float64
	WStdev  float64
	WP5     float64
	WQ25    float64
	WMedian float64
	WQ75    float64
	WP95    float64
}

// SMCSampleStats is one MCMC sample — every .smc block file sharing an iteration.
type SMCSampleStats struct {
	Iteration string // parsed from the filename; "" when the files carry no iteration
	Files     int
	Skipped   int64 // TREE lines whose tree had <= K leaves, or that failed to parse

	TMRKA SMCDist // the statistic under test, K as configured (4 = published)
	TMRCA SMCDist // K=1, for the TMR4A/TMRCA ratio cross-check
}

// RatioWeighted is TMR-K-A / TMRCA on the bp-weighted medians. Swamidass reports
// 0.38 against a Kingman expectation of ~0.24 at n=108 and never explains the gap.
// Computing it here per MCMC sample shows whether the gap is stable across the
// posterior or is a property of the one sample that was examined.
func (s *SMCSampleStats) RatioWeighted() float64 {
	if s.TMRCA.WMedian == 0 {
		return math.NaN()
	}
	return s.TMRKA.WMedian / s.TMRCA.WMedian
}

// SMCScanResult is the whole scan: one entry per MCMC sample, plus the
// across-sample spread that the published analysis never computed.
type SMCScanResult struct {
	K       int
	GenTime float64
	Samples []SMCSampleStats

	// Across-sample spread of the genome-wide bp-weighted median TMR-K-A, in
	// generations. This is axis (2) — the posterior. With one sample it is a
	// point and MedianRange is zero, which is exactly the published situation.
	MedianMean   float64
	MedianStdev  float64
	MedianMin    float64
	MedianMax    float64
	MedianAcross float64 // median across samples of the per-sample weighted median
}

// SMCScanOptions configures ScanSMC.
type SMCScanOptions struct {
	Paths   []string // files and/or directories; directories are walked for *.smc(.gz)
	K       int      // 4 reproduces the published statistic
	GenTime float64  // years per generation; 25 is the published constant
	Verbose bool

	// IterFromPath matches the iteration against the whole path rather than the
	// basename. Some ARGweaver layouts put the iteration in a directory
	// (.../2400/chr1.smc.gz) rather than the filename, and grouping every block
	// under one nameless sample would silently pool the whole posterior into a
	// single "distribution" — which would look like an answer and be a mistake.
	IterFromPath bool
}

// smcIterRe pulls the MCMC iteration out of an ARGweaver output filename, e.g.
// `chr1.2400.smc.gz` or `out.iter_2400.smc`. Files that match nothing are grouped
// together under "" so a single hand-picked file still works.
var smcIterRe = regexp.MustCompile(`(?:^|[._])(?:iter[_-]?)?(\d+)\.smc(?:\.gz)?$`)

// smcDirIterRe is the -argweaver-smc-iter-from-path fallback: an iteration held in
// a directory component, e.g. .../out/2400/chr1.smc.gz or .../iter_2400/....
var smcDirIterRe = regexp.MustCompile(`/(?:iter[_-]?)?(\d+)/[^/]*\.smc(?:\.gz)?$`)

// openMaybeGz opens a path, transparently decompressing .gz.
func openMaybeGz(path string) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if !strings.HasSuffix(strings.ToLower(path), ".gz") {
		return f, nil
	}
	zr, err := gzip.NewReader(f)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return &gzReadCloser{zr: zr, f: f}, nil
}

type gzReadCloser struct {
	zr *gzip.Reader
	f  *os.File
}

func (g *gzReadCloser) Read(p []byte) (int, error) { return g.zr.Read(p) }
func (g *gzReadCloser) Close() error {
	err := g.zr.Close()
	if cerr := g.f.Close(); err == nil {
		err = cerr
	}
	return err
}

// collectSMCFiles expands the option paths into concrete .smc/.smc.gz files.
func collectSMCFiles(paths []string) ([]string, error) {
	var out []string
	for _, p := range paths {
		info, err := os.Stat(p)
		if err != nil {
			return nil, err
		}
		if !info.IsDir() {
			out = append(out, p)
			continue
		}
		err = filepath.Walk(p, func(path string, fi os.FileInfo, err error) error {
			if err != nil {
				return err
			}
			if fi.IsDir() {
				return nil
			}
			l := strings.ToLower(path)
			if strings.HasSuffix(l, ".smc") || strings.HasSuffix(l, ".smc.gz") {
				out = append(out, path)
			}
			return nil
		})
		if err != nil {
			return nil, err
		}
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("no .smc or .smc.gz files found under %s", strings.Join(paths, ", "))
	}
	sort.Strings(out)
	return out, nil
}

// smcAccum is the per-sample working state: one bucket map per statistic.
type smcAccum struct {
	kBuckets map[float64]*smcBucket
	rBuckets map[float64]*smcBucket
	files    int
	skipped  int64
}

func newSMCAccum() *smcAccum {
	return &smcAccum{
		kBuckets: map[float64]*smcBucket{},
		rBuckets: map[float64]*smcBucket{},
	}
}

func (a *smcAccum) add(m map[float64]*smcBucket, v float64, span int64) {
	b := m[v]
	if b == nil {
		b = &smcBucket{}
		m[v] = b
	}
	b.trees++
	b.spanB += span
}

// scanSMCFile streams one .smc file into the accumulator.
func scanSMCFile(path string, k int, a *smcAccum) error {
	rc, err := openMaybeGz(path)
	if err != nil {
		return err
	}
	defer rc.Close()

	sc := bufio.NewScanner(rc)
	// Local trees for 108 haploids run to tens of kB; the default 64 kB token
	// limit is not enough and the failure mode is a silent short read.
	sc.Buffer(make([]byte, 0, 1<<20), 64<<20)

	for sc.Scan() {
		line := sc.Text()
		if !strings.HasPrefix(line, "TREE\t") {
			continue
		}
		tok := strings.SplitN(line, "\t", 4)
		if len(tok) < 4 {
			a.skipped++
			continue
		}
		start, err1 := strconv.ParseInt(strings.TrimSpace(tok[1]), 10, 64)
		end, err2 := strconv.ParseInt(strings.TrimSpace(tok[2]), 10, 64)
		if err1 != nil || err2 != nil || end < start {
			a.skipped++
			continue
		}
		// .smc TREE coordinates are 1-based inclusive.
		span := end - start + 1

		root, leaves, err := parseNewick(strings.TrimSpace(tok[3]))
		if err != nil {
			a.skipped++
			continue
		}
		if tk, _, ok := tmrKA(root, leaves, k); ok {
			a.add(a.kBuckets, tk, span)
		} else {
			a.skipped++
		}
		if tr, _, ok := tmrKA(root, leaves, 1); ok {
			a.add(a.rBuckets, tr, span)
		}
	}
	if err := sc.Err(); err != nil {
		return fmt.Errorf("%s: %w", path, err)
	}
	a.files++
	return nil
}

// summariseBuckets turns a bucket map into an SMCDist, computing both the
// unweighted and the bp-weighted summary exactly.
func summariseBuckets(m map[float64]*smcBucket) SMCDist {
	var d SMCDist
	if len(m) == 0 {
		return d
	}
	vals := make([]float64, 0, len(m))
	for v := range m {
		vals = append(vals, v)
	}
	sort.Float64s(vals)
	d.DistinctVals = len(vals)
	d.Min = vals[0]
	d.Max = vals[len(vals)-1]

	for _, v := range vals {
		b := m[v]
		d.Trees += b.trees
		d.SpanBp += b.spanB
		d.Mean += v * float64(b.trees)
		d.WMean += v * float64(b.spanB)
	}
	if d.Trees > 0 {
		d.Mean /= float64(d.Trees)
	}
	if d.SpanBp > 0 {
		d.WMean /= float64(d.SpanBp)
	}
	for _, v := range vals {
		b := m[v]
		d.Stdev += float64(b.trees) * (v - d.Mean) * (v - d.Mean)
		d.WStdev += float64(b.spanB) * (v - d.WMean) * (v - d.WMean)
	}
	if d.Trees > 1 {
		d.Stdev = math.Sqrt(d.Stdev / float64(d.Trees-1))
	}
	if d.SpanBp > 1 {
		d.WStdev = math.Sqrt(d.WStdev / float64(d.SpanBp-1))
	}

	cw := func(pick func(*smcBucket) int64) func(float64) float64 {
		var total int64
		for _, v := range vals {
			total += pick(m[v])
		}
		return func(q float64) float64 {
			if total == 0 {
				return 0
			}
			// Lower-weighted-quantile: the smallest grid value whose cumulative
			// weight reaches q. On a discretised grid this is the honest answer —
			// interpolating between grid points would invent precision the
			// sampler never had.
			target := q * float64(total)
			var run int64
			for _, v := range vals {
				run += pick(m[v])
				if float64(run) >= target {
					return v
				}
			}
			return vals[len(vals)-1]
		}
	}
	uq := cw(func(b *smcBucket) int64 { return b.trees })
	wq := cw(func(b *smcBucket) int64 { return b.spanB })

	d.P5, d.Q25, d.Median, d.Q75, d.P95 = uq(0.05), uq(0.25), uq(0.50), uq(0.75), uq(0.95)
	d.WP5, d.WQ25, d.WMedian, d.WQ75, d.WP95 = wq(0.05), wq(0.25), wq(0.50), wq(0.75), wq(0.95)
	return d
}

// ScanSMC reads every .smc file, groups them by MCMC iteration, and reports the
// per-sample genome-wide distribution plus the across-sample spread.
func ScanSMC(opts SMCScanOptions) (*SMCScanResult, error) {
	if opts.K < 1 {
		opts.K = 4
	}
	if opts.GenTime <= 0 {
		opts.GenTime = 25
	}
	files, err := collectSMCFiles(opts.Paths)
	if err != nil {
		return nil, err
	}

	// Group by iteration so that one "sample" is a whole genome-wide draw
	// assembled from its per-block files, which is the unit the published
	// statistic is computed over.
	groups := map[string]*smcAccum{}
	var order []string
	for _, f := range files {
		iter := ""
		subject := filepath.Base(f)
		if opts.IterFromPath {
			subject = filepath.ToSlash(f)
		}
		if m := smcIterRe.FindStringSubmatch(subject); m != nil {
			iter = m[1]
		} else if opts.IterFromPath {
			if m := smcDirIterRe.FindStringSubmatch(filepath.ToSlash(f)); m != nil {
				iter = m[1]
			}
		}
		a := groups[iter]
		if a == nil {
			a = newSMCAccum()
			groups[iter] = a
			order = append(order, iter)
		}
		if opts.Verbose {
			fmt.Printf("  scanning %s (iteration %q)\n", f, iter)
		}
		if err := scanSMCFile(f, opts.K, a); err != nil {
			return nil, err
		}
	}
	sort.Slice(order, func(i, j int) bool {
		a, ea := strconv.Atoi(order[i])
		b, eb := strconv.Atoi(order[j])
		if ea == nil && eb == nil {
			return a < b
		}
		return order[i] < order[j]
	})

	res := &SMCScanResult{K: opts.K, GenTime: opts.GenTime}
	var medians []float64
	for _, it := range order {
		a := groups[it]
		s := SMCSampleStats{
			Iteration: it,
			Files:     a.files,
			Skipped:   a.skipped,
			TMRKA:     summariseBuckets(a.kBuckets),
			TMRCA:     summariseBuckets(a.rBuckets),
		}
		res.Samples = append(res.Samples, s)
		if s.TMRKA.Trees > 0 {
			medians = append(medians, s.TMRKA.WMedian)
		}
	}
	if len(medians) > 0 {
		res.MedianMean, res.MedianStdev = meanStdev(medians)
		sorted := append([]float64(nil), medians...)
		sort.Float64s(sorted)
		res.MedianMin = sorted[0]
		res.MedianMax = sorted[len(sorted)-1]
		res.MedianAcross = quantile(sorted, 0.50)
	}
	return res, nil
}

// PrintSMCScan renders the scan. Everything is shown in generations AND years,
// because every year-denominated claim inherits the 25 yr/generation constant and
// the reader should be able to see where it entered.
func PrintSMCScan(res *SMCScanResult) {
	y := func(g float64) float64 { return g * res.GenTime }
	fmt.Printf("\nARGweaver .smc scan — TMR-%dA over local trees (%.0f yr/generation)\n",
		res.K, res.GenTime)
	fmt.Printf("%s\n", strings.Repeat("=", 78))

	for _, s := range res.Samples {
		label := s.Iteration
		if label == "" {
			label = "(no iteration in filename)"
		}
		fmt.Printf("\nMCMC sample %s — %d file(s), %d local trees, %.1f Mb spanned",
			label, s.Files, s.TMRKA.Trees, float64(s.TMRKA.SpanBp)/1e6)
		if s.Skipped > 0 {
			fmt.Printf(", %d skipped", s.Skipped)
		}
		fmt.Printf("\n  time grid: %d distinct TMR-%dA values occupied\n", s.TMRKA.DistinctVals, res.K)

		d := s.TMRKA
		fmt.Printf("  %-26s %12s %12s\n", "TMR-"+strconv.Itoa(res.K)+"A", "generations", "years")
		fmt.Printf("  %-26s %12.1f %12.0f\n", "median (bp-weighted)", d.WMedian, y(d.WMedian))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "median (unweighted)", d.Median, y(d.Median))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "mean (bp-weighted)", d.WMean, y(d.WMean))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "mean (unweighted)", d.Mean, y(d.Mean))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "stdev (bp-weighted)", d.WStdev, y(d.WStdev))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "5th pct (bp-weighted)", d.WP5, y(d.WP5))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "25th pct (bp-weighted)", d.WQ25, y(d.WQ25))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "75th pct (bp-weighted)", d.WQ75, y(d.WQ75))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "95th pct (bp-weighted)", d.WP95, y(d.WP95))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "min", d.Min, y(d.Min))
		fmt.Printf("  %-26s %12.1f %12.0f\n", "max", d.Max, y(d.Max))

		r := s.TMRCA
		fmt.Printf("  %-26s %12.1f %12.0f\n", "TMRCA median (bp-wtd)", r.WMedian, y(r.WMedian))
		fmt.Printf("  %-26s %12.3f\n", "TMR-"+strconv.Itoa(res.K)+"A / TMRCA", s.RatioWeighted())
	}

	fmt.Printf("\n%s\n", strings.Repeat("-", 78))
	if len(res.Samples) < 2 {
		fmt.Printf("ACROSS MCMC SAMPLES: only %d sample supplied, so there is no posterior\n"+
			"spread to report. This is exactly the published situation — the 495 kya\n"+
			"figure rests on iteration 2400 alone. Supply more .smc files to fill it in.\n",
			len(res.Samples))
		return
	}
	fmt.Printf("ACROSS %d MCMC SAMPLES — genome-wide bp-weighted median TMR-%dA\n",
		len(res.Samples), res.K)
	fmt.Printf("  %-26s %12.1f %12.0f\n", "mean", res.MedianMean, y(res.MedianMean))
	fmt.Printf("  %-26s %12.1f %12.0f\n", "median", res.MedianAcross, y(res.MedianAcross))
	fmt.Printf("  %-26s %12.1f %12.0f\n", "stdev", res.MedianStdev, y(res.MedianStdev))
	fmt.Printf("  %-26s %12.1f %12.0f\n", "min", res.MedianMin, y(res.MedianMin))
	fmt.Printf("  %-26s %12.1f %12.0f\n", "max", res.MedianMax, y(res.MedianMax))
	fmt.Printf("  %-26s %12.1f %12.0f\n", "range", res.MedianMax-res.MedianMin,
		y(res.MedianMax-res.MedianMin))
}

// WriteSMCScan writes one row per MCMC sample, so the across-sample spread can be
// replotted without re-reading hundreds of GB.
func WriteSMCScan(path string, res *SMCScanResult) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()

	k := strconv.Itoa(res.K)
	fmt.Fprintf(w, "iteration,files,trees,span_bp,distinct_grid_values,"+
		"tmr%sa_wmedian_gen,tmr%sa_median_gen,tmr%sa_wmean_gen,tmr%sa_mean_gen,tmr%sa_wstdev_gen,"+
		"tmr%sa_wp5_gen,tmr%sa_wq25_gen,tmr%sa_wq75_gen,tmr%sa_wp95_gen,tmr%sa_min_gen,tmr%sa_max_gen,"+
		"tmrca_wmedian_gen,ratio_tmr%sa_tmrca,gen_time_yr,tmr%sa_wmedian_yr\n",
		k, k, k, k, k, k, k, k, k, k, k, k, k)
	for _, s := range res.Samples {
		d := s.TMRKA
		fmt.Fprintf(w, "%s,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.6f,%.2f,%.1f\n",
			s.Iteration, s.Files, d.Trees, d.SpanBp, d.DistinctVals,
			d.WMedian, d.Median, d.WMean, d.Mean, d.WStdev,
			d.WP5, d.WQ25, d.WQ75, d.WP95, d.Min, d.Max,
			s.TMRCA.WMedian, s.RatioWeighted(), res.GenTime, d.WMedian*res.GenTime)
	}
	return nil
}
