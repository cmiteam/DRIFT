package analysis

// The sampling bootstrap — TMR4A.md W9c stage 1, the fourth variance component.
//
// THE QUESTION. `arg-sample` has no sample-size setting; n is whatever the `.sites` file
// holds, frozen upstream by DRIFT's `vcf_sample_size`. So the 20 individuals exported for
// `createdvar` are one draw from 2,000, and every structure statistic the study leads on
// — the cliff position, the zero-skip K, K4/K10 — is a count of divergent classes *in
// that draw*. Miss a class and the cliff moves by one, indistinguishably from the
// recombinant-mosaic explanation §9.12 gave for `createdvar_excl`.
//
// WHY SEED REPLICATES CANNOT ANSWER IT. `SampleIDs` draws from the same seeded stream as
// the simulation, so a new `rng_seed` changes the realization and the draw together. The
// two are confounded by construction. Isolating the sampling component requires the
// realization held fixed and only the draw varied — which is what exporting at
// `vcf_sample_size = 0` and resampling here does. That export draws no RNG at all
// (pkg/analysis/ibd.go:202) and every export runs after the simulation completes, so it
// is the same realization, more fully written down.
//
// WHY IT COSTS NO ARGweaver TIME. §9.12 showed the class count predicts the cliff
// position exactly on `createdvar` (7 = 7). So the distribution of class counts over
// redraws is the distribution of cliff positions under resampling, obtained without
// running the MCMC once. Stage 2 then spends ~27 h confirming three points of that
// distribution rather than discovering its shape at ~9 h a draw.
//
// TWO STATISTICS PER DRAW, AND THEY ARE NOT THE SAME NUMBER. `ClassesRecluster` clusters
// the draw's own haplotypes from scratch, which is what an analyst holding only that
// export would do and therefore what the study actually did. `ClassesRepresented` counts
// how many of the full population's classes the draw touches. They diverge exactly when
// clustering 40 haplotypes gives a different answer from clustering 4,000 — chaining
// through a mosaic that the smaller draw happens to include or exclude — and that
// divergence is a measurement of how fragile the cliff is, not an inconsistency.

import (
	"encoding/csv"
	"fmt"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

// BootstrapOptions configures one stage-1 run.
type BootstrapOptions struct {
	// NIndividuals per draw. 20 is the study's `vcf_sample_size`, giving the 40
	// haplotypes every archived export carries.
	NIndividuals int
	Draws        int

	// Seed for the resampling only. This deliberately uses its own generator rather
	// than pkg/utils' simulation stream: the whole point of stage 1 is that the
	// realization is fixed and only the draw varies, so the draw must NOT come from
	// the stream that made the realization. Phase 0's ban on ad-hoc `math/rand` is
	// about reproducing simulations; this is a post-hoc analysis whose seed is an
	// explicit input and is recorded in the output.
	Seed int64

	// Threshold for single linkage. NaN selects the parameter-free per-draw
	// largest-gap cut (LargestGapCut), which is the headline rule.
	Threshold float64

	// MinHeight guards LargestGapCut against cutting at the gap between identical
	// haplotypes and the first real merge. Ignored when Threshold is set.
	MinHeight float64
}

// DrawResult is one resampled draw.
type DrawResult struct {
	Draw               int
	ClassesRecluster   int
	ClassesRepresented int
	Cut                float64
	GapBelow           float64
	Simpson            float64
	Sizes              []int
	MissedClasses      []int // full-population class ranks absent from this draw

	// Individuals drawn, by name. Recorded because the draws are not merely a
	// distribution to look at: stage 2 has to EXPORT three specific draws and run them
	// through arg-sample, so a draw that cannot be named cannot be re-run. The seed and
	// the draw index reproduce it too, but only if every option is identical, and that
	// is a fragile thing to ask of a hand-off months later.
	Individuals []string
}

// BootstrapResult is the whole of stage 1 for one case.
type BootstrapResult struct {
	Path    string
	Window  int
	Opts    BootstrapOptions
	NHaps   int
	NInds   int
	NSites  int
	MultiAl int

	// Full is the clustering of every exported haplotype: the population's true class
	// structure, and the per-class frequency vector that §9.15 records as computed
	// inside summariseSources and then discarded.
	Full      *Clustering
	FullCut   float64
	FullGap   float64
	FullSweep []SweepPoint

	Draws []DrawResult

	// CountHistogram maps a class count to the number of draws returning it — the
	// distribution of cliff positions under resampling.
	CountHistogram map[int]int
	// ClassHitRate maps a full-population class rank (0 = largest) to the fraction of
	// draws that sampled it at least once.
	ClassHitRate map[int]float64
}

// RunBootstrap performs stage 1 on an already-loaded sample.
func RunBootstrap(h *HaploSample, path string, window int, opts BootstrapOptions) (*BootstrapResult, error) {
	if opts.NIndividuals <= 0 {
		opts.NIndividuals = 20
	}
	if opts.Draws <= 0 {
		opts.Draws = 1000
	}

	all := make([]int, len(h.Names))
	for i := range all {
		all[i] = i
	}
	indNames, indCols := h.IndividualIndex()
	if opts.NIndividuals > len(indNames) {
		return nil, fmt.Errorf("asked for %d individuals but the file holds %d",
			opts.NIndividuals, len(indNames))
	}

	res := &BootstrapResult{
		Path: path, Window: window, Opts: opts,
		NHaps: len(h.Names), NInds: len(indNames), NSites: h.NSites, MultiAl: h.MultiAllelic,
		CountHistogram: map[int]int{}, ClassHitRate: map[int]float64{},
	}

	// The full population first: it supplies both the frequency vector and the
	// reference class membership that ClassesRepresented is counted against.
	fullEdges := h.mstEdges(all)
	heights := make([]float64, len(fullEdges))
	for i, e := range fullEdges {
		heights[i] = e.dist
	}
	res.FullCut, res.FullGap = LargestGapCut(heights, opts.MinHeight)
	if !math.IsNaN(opts.Threshold) {
		res.FullCut = opts.Threshold
	}
	res.Full = clusterFromEdges(all, fullEdges, res.FullCut)
	res.FullSweep = h.SweepThresholds(all, h.ThresholdGrid())

	// Rank each haplotype's full-population class, so a draw can report which classes
	// it missed by rank rather than by an arbitrary member index.
	classOf := make([]int, len(h.Names))
	for rank, members := range res.Full.Classes {
		for _, m := range members {
			classOf[m] = rank
		}
	}

	rng := rand.New(rand.NewSource(opts.Seed))
	hits := make([]int, len(res.Full.Classes))

	for d := 0; d < opts.Draws; d++ {
		// Draw individuals without replacement, matching SampleIDs' uniform draw.
		perm := rng.Perm(len(indNames))[:opts.NIndividuals]
		var subset []int
		for _, p := range perm {
			subset = append(subset, indCols[indNames[p]]...)
		}
		sort.Ints(subset)

		edges := h.mstEdges(subset)
		hh := make([]float64, len(edges))
		for i, e := range edges {
			hh[i] = e.dist
		}
		cut := opts.Threshold
		gap := math.NaN()
		if math.IsNaN(cut) {
			cut, gap = LargestGapCut(hh, opts.MinHeight)
		}
		c := clusterFromEdges(subset, edges, cut)
		if math.IsNaN(gap) {
			gap = c.GapBelow
		}

		seen := map[int]bool{}
		for _, s := range subset {
			seen[classOf[s]] = true
		}
		var missed []int
		for rank := range res.Full.Classes {
			if seen[rank] {
				hits[rank]++
			} else {
				missed = append(missed, rank)
			}
		}

		drawn := make([]string, 0, opts.NIndividuals)
		for _, p := range perm {
			drawn = append(drawn, indNames[p])
		}
		sort.Strings(drawn)

		dr := DrawResult{
			Draw: d, ClassesRecluster: len(c.Classes), ClassesRepresented: len(seen),
			Cut: cut, GapBelow: gap, Simpson: c.Simpson(), Sizes: c.Sizes(),
			MissedClasses: missed, Individuals: drawn,
		}
		res.Draws = append(res.Draws, dr)
		res.CountHistogram[dr.ClassesRecluster]++
	}

	for rank, n := range hits {
		res.ClassHitRate[rank] = float64(n) / float64(opts.Draws)
	}
	return res, nil
}

// --- Reporting -------------------------------------------------------------------------

// quantile returns the q-th quantile of an already-sorted int slice, nearest-rank.
func quantileInt(sorted []int, q float64) int {
	if len(sorted) == 0 {
		return 0
	}
	idx := int(q * float64(len(sorted)-1))
	return sorted[idx]
}

// PrintBootstrap writes the end-of-run summary, in the shape the other analyses use.
func PrintBootstrap(r *BootstrapResult) {
	if r == nil {
		return
	}
	fmt.Printf("\nHaplotype clustering — %s\n", filepath.Base(r.Path))
	if r.Window > 0 {
		fmt.Printf("  window            1-%d bp (%d sites retained)\n", r.Window, r.NSites)
	} else {
		fmt.Printf("  window            whole file (%d sites)\n", r.NSites)
	}
	fmt.Printf("  sample            %d haplotypes / %d individuals\n", r.NHaps, r.NInds)
	if r.MultiAl > 0 {
		fmt.Printf("  WARNING: %d sites have more than two alleles; packing is biallelic\n", r.MultiAl)
	}

	fmt.Printf("\n  FULL POPULATION, cut %.4f (largest merge gap %.4f)\n", r.FullCut, r.FullGap)
	fmt.Printf("    classes         %d\n", len(r.Full.Classes))
	fmt.Printf("    sizes           %v\n", r.Full.Sizes())
	fmt.Printf("    effective (1/sum p^2) %.2f\n", r.Full.Simpson())

	pc, lo, hi := Plateau(r.FullSweep)
	fmt.Printf("    widest plateau  %d classes over thresholds %.4f-%.4f\n", pc, lo, hi)

	fmt.Printf("\n  RESAMPLING %d individuals, %d draws, seed %d\n",
		r.Opts.NIndividuals, r.Opts.Draws, r.Opts.Seed)

	counts := make([]int, 0, len(r.Draws))
	repr := make([]int, 0, len(r.Draws))
	for _, d := range r.Draws {
		counts = append(counts, d.ClassesRecluster)
		repr = append(repr, d.ClassesRepresented)
	}
	sort.Ints(counts)
	sort.Ints(repr)

	mean := 0.0
	for _, c := range counts {
		mean += float64(c)
	}
	mean /= float64(len(counts))

	fmt.Printf("    classes per draw (re-clustered): mean %.2f  median %d  p05 %d  p95 %d  min %d  max %d\n",
		mean, quantileInt(counts, 0.5), quantileInt(counts, 0.05), quantileInt(counts, 0.95),
		counts[0], counts[len(counts)-1])
	fmt.Printf("    full classes represented:        median %d  p05 %d  min %d  max %d\n",
		quantileInt(repr, 0.5), quantileInt(repr, 0.05), repr[0], repr[len(repr)-1])

	fmt.Println("    distribution of class counts (= cliff positions under resampling):")
	keys := make([]int, 0, len(r.CountHistogram))
	for k := range r.CountHistogram {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	for _, k := range keys {
		n := r.CountHistogram[k]
		frac := float64(n) / float64(len(r.Draws))
		fmt.Printf("      %2d classes  %6d draws  %5.1f%%  %s\n",
			k, n, 100*frac, strings.Repeat("#", int(frac*50+0.5)))
	}

	fmt.Println("    per-class capture rate (rank: size, frequency, P(sampled)):")
	for rank, members := range r.Full.Classes {
		p := float64(len(members)) / float64(r.NHaps)
		fmt.Printf("      %2d  size %5d  freq %.4f  captured %5.1f%%\n",
			rank, len(members), p, 100*r.ClassHitRate[rank])
	}
}

// WriteBootstrap writes the four CSVs stage 1 produces, all prefixed alike.
func WriteBootstrap(prefix string, r *BootstrapResult) ([]string, error) {
	var written []string

	// 1. The threshold sweep — so no class count in this study is ever quoted without
	//    the evidence that it is or is not threshold-robust.
	sweepPath := prefix + "_sweep.csv"
	if err := writeCSV(sweepPath,
		[]string{"threshold", "classes", "effective_classes", "sizes"},
		func(w *csv.Writer) error {
			for _, p := range r.FullSweep {
				if err := w.Write([]string{
					strconv.FormatFloat(p.Threshold, 'g', -1, 64),
					strconv.Itoa(p.Classes),
					strconv.FormatFloat(p.Simpson, 'f', 4, 64),
					joinInts(p.Sizes),
				}); err != nil {
					return err
				}
			}
			return nil
		}); err != nil {
		return written, err
	}
	written = append(written, sweepPath)

	// 2. The per-class frequency vector — the object §9.15 records as computed and
	//    discarded, now recoverable from an output file.
	classPath := prefix + "_classes.csv"
	if err := writeCSV(classPath,
		[]string{"rank", "size", "frequency", "capture_rate", "example_haplotype"},
		func(w *csv.Writer) error {
			for rank, members := range r.Full.Classes {
				example := ""
				if len(members) > 0 {
					example = strconv.Itoa(members[0])
				}
				if err := w.Write([]string{
					strconv.Itoa(rank),
					strconv.Itoa(len(members)),
					strconv.FormatFloat(float64(len(members))/float64(r.NHaps), 'f', 6, 64),
					strconv.FormatFloat(r.ClassHitRate[rank], 'f', 4, 64),
					example,
				}); err != nil {
					return err
				}
			}
			return nil
		}); err != nil {
		return written, err
	}
	written = append(written, classPath)

	// 3. Every draw, so the distribution can be re-analysed without re-running.
	drawPath := prefix + "_draws.csv"
	if err := writeCSV(drawPath,
		[]string{"draw", "classes_reclustered", "full_classes_represented", "cut",
			"gap_below_cut", "effective_classes", "sizes", "missed_class_ranks",
			"individuals"},
		func(w *csv.Writer) error {
			for _, d := range r.Draws {
				if err := w.Write([]string{
					strconv.Itoa(d.Draw),
					strconv.Itoa(d.ClassesRecluster),
					strconv.Itoa(d.ClassesRepresented),
					strconv.FormatFloat(d.Cut, 'f', 6, 64),
					strconv.FormatFloat(d.GapBelow, 'f', 6, 64),
					strconv.FormatFloat(d.Simpson, 'f', 4, 64),
					joinInts(d.Sizes),
					joinInts(d.MissedClasses),
					strings.Join(d.Individuals, " "),
				}); err != nil {
					return err
				}
			}
			return nil
		}); err != nil {
		return written, err
	}
	written = append(written, drawPath)

	// 4. The histogram on its own, because it is the deliverable: the distribution of
	//    cliff positions under resampling.
	histPath := prefix + "_cliff_distribution.csv"
	if err := writeCSV(histPath,
		[]string{"classes", "draws", "fraction"},
		func(w *csv.Writer) error {
			keys := make([]int, 0, len(r.CountHistogram))
			for k := range r.CountHistogram {
				keys = append(keys, k)
			}
			sort.Ints(keys)
			for _, k := range keys {
				if err := w.Write([]string{
					strconv.Itoa(k),
					strconv.Itoa(r.CountHistogram[k]),
					strconv.FormatFloat(float64(r.CountHistogram[k])/float64(len(r.Draws)), 'f', 6, 64),
				}); err != nil {
					return err
				}
			}
			return nil
		}); err != nil {
		return written, err
	}
	written = append(written, histPath)

	return written, nil
}

func joinInts(v []int) string {
	parts := make([]string, len(v))
	for i, x := range v {
		parts[i] = strconv.Itoa(x)
	}
	return strings.Join(parts, " ")
}

func writeCSV(path string, header []string, body func(*csv.Writer) error) error {
	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()
	w := csv.NewWriter(file)
	defer w.Flush()
	if err := w.Write(header); err != nil {
		return err
	}
	if err := body(w); err != nil {
		return err
	}
	w.Flush()
	return w.Error()
}
