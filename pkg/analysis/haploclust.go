package analysis

// Haplotype clustering by pairwise divergence — the DRIFT-side cliff predictor, and the
// sampling bootstrap built on it (TMR4A.md W9c stage 1).
//
// WHAT THIS MEASURES AND WHY IT IS NOT AN ARGweaver RESULT. Given an exported `.sites`
// file, this groups the haplotype columns into divergence classes and counts them. It
// never reads an ARG, never runs arg-sample, and consults no created-allele label: it is
// a measurement on the exported sequences alone. §9.12 established that on `createdvar`
// the class count equals the inferred cliff position exactly (7 = 7), which is what makes
// it a *predictor* — a case's cliff can be computed here, for free, before any MCMC is
// booked.
//
// WHY SINGLE LINKAGE. Created alleles are separated by a designed divergence `d` and
// their descendants differ from one another only by the handful of de-novo mutations that
// 500 years affords. So the true structure is tight balls at radius ~0 separated by gulfs
// of ~2d, and single linkage is the method that recovers exactly that geometry without
// assuming ball shape or equal class sizes. Its known weakness — chaining through
// intermediate points — is not a defect here but the *signal*: a recombinant mosaic is
// literally an intermediate point, and §9.12's `createdvar_excl` residual is one.
//
// WHY THE THRESHOLD IS SWEPT AND NEVER ASSUMED. Calibration against the two archived runs
// (2026-09-21) reproduces both published class counts, but not at a single threshold:
//
//	threshold   createdvar          createdvar_excl
//	0.12        7  (= its cliff)    10 (sizes 8 6 5 5 5 4 3 2 1 1, the §9.15 vector)
//	0.13        7                    9 (= its cliff; the 0.128 mosaic merges in)
//
// The `_excl` mosaic sits 0.128 from its parent class against a 0.2204 median
// between-class distance, so §9.12's "predicted 10, answered 9" discrepancy IS this
// threshold sensitivity, and a tool that hard-coded one cut would have hidden it. Every
// report below therefore carries the full sweep, and the headline number comes from the
// parameter-free largest-gap rule (LargestGapCut) rather than from a chosen constant.
//
// WHAT THE BOOTSTRAP ADDS. §10.8 listed three variance components and missed a fourth:
// arg-sample has no sample-size setting, so n is frozen into the `.sites` file upstream
// and seed replicates cannot vary it (`SampleIDs` draws from the same seeded stream as
// the simulation, confounding the draw with the realization). Holding one realization
// fixed and re-drawing only the individuals isolates it. Because the class count predicts
// the cliff, the distribution of class counts over redraws IS the distribution of cliff
// positions under resampling — obtained with no ARGweaver time at all.

import (
	"bufio"
	"fmt"
	"math"
	"math/bits"
	"os"
	"sort"
	"strconv"
	"strings"
)

// HaploSample is the haplotype columns of one `.sites` file over one window, packed for
// fast pairwise comparison.
//
// Sites are stored as a minor-allele indicator bitset per haplotype rather than as the
// original characters: every comparison this file performs is "do these two differ", so
// one bit per site per haplotype is sufficient and lets the distance be a popcount over
// XORed words instead of a byte loop. At the full-population scale W9c stage 1 needs
// (4,000 haplotypes, ~8M pairs) that difference is the difference between seconds and
// hours.
type HaploSample struct {
	Names  []string // haplotype column names, in file order
	Inds   []string // individual name per haplotype — the name with its _1/_2 suffix removed
	NSites int      // sites retained in the window
	Region string   // the REGION line, carried through for provenance

	// MultiAllelic counts sites with more than two observed characters. The packing is
	// biallelic, so such a site is recorded as "differs from the major allele", which
	// merges two minor alleles into one. On a DRIFT export this is zero — the genome is
	// a bit array — and it is reported so a non-zero count on foreign data is visible
	// rather than silent.
	MultiAllelic int

	// Scale is the denominator every distance is divided by. It defaults to NSites,
	// which makes a distance the fraction of window sites at which two haplotypes
	// differ — the metric §9.12 used, reproduced exactly by the calibration above.
	//
	// THAT METRIC IS NOT COMPARABLE ACROSS SAMPLE SIZES, and W9c needs one that is. A
	// site is written only when it is variable within the sample, so a bigger sample
	// emits more sites and the denominator grows with n. The NUMERATOR does not: if
	// two haplotypes differ at a site, that site is variable in any sample containing
	// them both, so their differing-site COUNT is identical in every file that holds
	// them. Setting Scale to 1 therefore puts every sample size on one ruler, and is
	// what the n-scaling series must use.
	Scale float64

	words int        // uint64 words per haplotype
	bits  [][]uint64 // minor-allele indicator, one row per haplotype
}

// SetScale changes the distance denominator. Single-linkage clustering is invariant
// under it — rescaling every distance by a constant reorders nothing — so this changes
// the NUMBERS a threshold is quoted in and never the classes a given cut produces.
func (h *HaploSample) SetScale(scale float64) {
	if scale > 0 {
		h.Scale = scale
	}
}

// IndividualIndex maps each individual name to its haplotype column indices, in file
// order. Individuals are the unit the bootstrap resamples, because that is the unit
// `SampleIDs` resamples: a diploid contributes both its strands or neither, and drawing
// haplotypes independently would understate the loss of rare classes by breaking the
// correlation between a genome's two copies.
func (h *HaploSample) IndividualIndex() (names []string, cols map[string][]int) {
	cols = map[string][]int{}
	for i, ind := range h.Inds {
		if _, seen := cols[ind]; !seen {
			names = append(names, ind)
		}
		cols[ind] = append(cols[ind], i)
	}
	return names, cols
}

// Distance is the fraction of retained sites at which two haplotypes differ.
//
// The denominator is the number of sites in the window, not the number of sites
// segregating within the pair, so distances from different subsets are on the same scale
// and a bootstrap draw's distances are comparable to the full population's.
func (h *HaploSample) Distance(i, j int) float64 {
	if h.NSites == 0 || h.Scale == 0 {
		return 0
	}
	a, b := h.bits[i], h.bits[j]
	diff := 0
	for w := 0; w < h.words; w++ {
		diff += bits.OnesCount64(a[w] ^ b[w])
	}
	return float64(diff) / h.Scale
}

// LoadSitesWindow reads an ARGweaver `.sites` file, keeping only sites whose position is
// at most maxPos (0 = keep the whole file).
//
// The window matters and is not cosmetic: arg-sample was run on `--region 1-2000000`
// while the exported file spans the whole 100 Mb genome, so clustering the whole file
// would measure a different object from the one the MCMC saw. Passing the region's upper
// bound here is what makes the class count comparable to the cliff.
func LoadSitesWindow(path string, maxPos int) (*HaploSample, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open sites file: %w", err)
	}
	defer file.Close()

	sc := bufio.NewScanner(file)
	sc.Buffer(make([]byte, 1<<20), 1<<24)

	h := &HaploSample{}
	var chars [][]byte // per-haplotype characters, accumulated then packed

	for sc.Scan() {
		line := sc.Text()
		switch {
		case strings.HasPrefix(line, "NAMES"):
			h.Names = strings.Split(line, "\t")[1:]
			chars = make([][]byte, len(h.Names))
			for _, n := range h.Names {
				// Column names are IND<id>_1 / IND<id>_2; the individual is the name
				// without that strand suffix. A name with no suffix is its own
				// individual, so foreign files still load.
				ind := n
				if k := strings.LastIndexByte(n, '_'); k > 0 {
					ind = n[:k]
				}
				h.Inds = append(h.Inds, ind)
			}
			continue
		case strings.HasPrefix(line, "REGION"):
			h.Region = line
			continue
		}

		tab := strings.IndexByte(line, '\t')
		if tab < 0 {
			continue
		}
		pos, err := strconv.Atoi(line[:tab])
		if err != nil {
			continue // not a site line
		}
		if maxPos > 0 && pos > maxPos {
			continue
		}
		bases := line[tab+1:]
		if len(bases) != len(h.Names) {
			return nil, fmt.Errorf("site at %d has %d characters for %d haplotypes",
				pos, len(bases), len(h.Names))
		}
		for i := 0; i < len(bases); i++ {
			chars[i] = append(chars[i], bases[i])
		}
		h.NSites++
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("read sites file: %w", err)
	}
	if len(h.Names) == 0 {
		return nil, fmt.Errorf("%s: no NAMES line", path)
	}

	h.pack(chars)
	h.Scale = float64(h.NSites)
	return h, nil
}

// pack turns per-haplotype characters into minor-allele indicator bitsets. The major
// allele is whichever character is commonest at that site; ties break on the character's
// byte value so the packing is deterministic and two runs of this tool agree.
func (h *HaploSample) pack(chars [][]byte) {
	n := len(h.Names)
	h.words = (h.NSites + 63) / 64
	h.bits = make([][]uint64, n)
	for i := range h.bits {
		h.bits[i] = make([]uint64, h.words)
	}

	counts := map[byte]int{}
	for s := 0; s < h.NSites; s++ {
		for k := range counts {
			delete(counts, k)
		}
		for i := 0; i < n; i++ {
			counts[chars[i][s]]++
		}
		if len(counts) > 2 {
			h.MultiAllelic++
		}
		var major byte
		best := -1
		for c, ct := range counts {
			if ct > best || (ct == best && c < major) {
				major, best = c, ct
			}
		}
		word, bit := s/64, uint(s%64)
		for i := 0; i < n; i++ {
			if chars[i][s] != major {
				h.bits[i][word] |= 1 << bit
			}
		}
	}
}

// --- Single-linkage clustering -------------------------------------------------------

// Clustering is one single-linkage cut of one haplotype subset.
type Clustering struct {
	Subset  []int   // haplotype indices clustered, in the caller's order
	Classes [][]int // class members, largest class first; indices into the ORIGINAL sample
	Cut     float64 // the threshold actually applied

	// MergeHeights is the sorted single-linkage merge distance of every merge
	// performed to reduce the subset to one class — the minimum spanning tree's edge
	// weights. The whole dendrogram is in here, so any threshold's class count is
	// recoverable from it without recomputing distances.
	MergeHeights []float64

	// GapBelow is the jump in merge height at the cut: the distance from the largest
	// merge accepted to the smallest merge refused. A wide gap means the class count is
	// robust to the threshold; a narrow one means it is not, and on this data that is
	// exactly the `createdvar_excl` mosaic.
	GapBelow float64
}

// Sizes returns the class sizes, largest first — the vector the study quotes.
func (c *Clustering) Sizes() []int {
	out := make([]int, len(c.Classes))
	for i, cl := range c.Classes {
		out[i] = len(cl)
	}
	return out
}

// Simpson returns the effective number of classes, 1/Σp², computed over the subset.
//
// This is the same statistic §9.15 applied to `created_pair_frac` at the population
// level, and it is reported beside the raw count because ten classes behaving like six
// and a half is the thing that makes a 40-haplotype draw lose one.
func (c *Clustering) Simpson() float64 {
	total := 0
	for _, cl := range c.Classes {
		total += len(cl)
	}
	if total == 0 {
		return 0
	}
	sum := 0.0
	for _, cl := range c.Classes {
		p := float64(len(cl)) / float64(total)
		sum += p * p
	}
	return 1 / sum
}

// mstEdges builds the single-linkage dendrogram of a subset as a minimum spanning tree,
// by Prim's algorithm on the implicit complete distance graph.
//
// Prim's is used rather than a distance matrix plus sort because the matrix is the memory
// bottleneck at full-population scale: 4,000 haplotypes would need 8M float64s per matrix
// where Prim's needs 4,000. Single-linkage clustering at any threshold t is exactly "cut
// every MST edge longer than t", which is why one MST pass answers every threshold in the
// sweep below.
type mstEdge struct {
	a, b int
	dist float64
}

func (h *HaploSample) mstEdges(subset []int) []mstEdge {
	n := len(subset)
	if n < 2 {
		return nil
	}
	inTree := make([]bool, n)
	best := make([]float64, n)
	from := make([]int, n)
	for i := range best {
		best[i] = math.Inf(1)
		from[i] = -1
	}
	best[0] = 0

	edges := make([]mstEdge, 0, n-1)
	for iter := 0; iter < n; iter++ {
		// Cheapest vertex not yet in the tree.
		next, nextDist := -1, math.Inf(1)
		for i := 0; i < n; i++ {
			if !inTree[i] && best[i] < nextDist {
				next, nextDist = i, best[i]
			}
		}
		if next < 0 {
			break
		}
		inTree[next] = true
		if from[next] >= 0 {
			edges = append(edges, mstEdge{a: subset[from[next]], b: subset[next], dist: nextDist})
		}
		for i := 0; i < n; i++ {
			if inTree[i] {
				continue
			}
			if d := h.Distance(subset[next], subset[i]); d < best[i] {
				best[i], from[i] = d, next
			}
		}
	}
	sort.Slice(edges, func(i, j int) bool { return edges[i].dist < edges[j].dist })
	return edges
}

// Cluster groups a haplotype subset by single linkage at threshold t.
func (h *HaploSample) Cluster(subset []int, t float64) *Clustering {
	return clusterFromEdges(subset, h.mstEdges(subset), t)
}

func clusterFromEdges(subset []int, edges []mstEdge, t float64) *Clustering {
	parent := map[int]int{}
	var find func(int) int
	find = func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	for _, s := range subset {
		parent[s] = s
	}

	c := &Clustering{Subset: subset, Cut: t, GapBelow: math.NaN()}
	lastAccepted := 0.0
	for _, e := range edges {
		c.MergeHeights = append(c.MergeHeights, e.dist)
		if e.dist > t {
			if math.IsNaN(c.GapBelow) {
				c.GapBelow = e.dist - lastAccepted
			}
			continue
		}
		lastAccepted = e.dist
		if a, b := find(e.a), find(e.b); a != b {
			parent[a] = b
		}
	}

	groups := map[int][]int{}
	for _, s := range subset {
		r := find(s)
		groups[r] = append(groups[r], s)
	}
	for _, members := range groups {
		sort.Ints(members)
		c.Classes = append(c.Classes, members)
	}
	// Largest first, then by first member, so the size vector is stable across runs.
	sort.Slice(c.Classes, func(i, j int) bool {
		if len(c.Classes[i]) != len(c.Classes[j]) {
			return len(c.Classes[i]) > len(c.Classes[j])
		}
		return c.Classes[i][0] < c.Classes[j][0]
	})
	return c
}

// LargestGapCut returns the parameter-free cut: the midpoint of the widest gap between
// consecutive single-linkage merge heights, and the class count it yields.
//
// This is the headline rule because it commits to no constant. The geometry it exploits
// is the one the created-allele model builds: within-class merges are near zero and
// between-class merges are near 2d, so the widest gap in the merge heights is the
// boundary between them. Where that geometry does not hold the gap is narrow, and a
// narrow gap is itself the finding — it says the class count is threshold-dependent and
// must be quoted with the sweep rather than alone.
//
// minHeight guards the degenerate case: with many exactly-identical haplotypes the
// largest gap can be the one between "0" and "the first real merge", which would report
// every distinct sequence as its own class. Gaps whose upper end is below minHeight are
// therefore ignored.
func LargestGapCut(edges []float64, minHeight float64) (cut float64, gap float64) {
	if len(edges) == 0 {
		return 0, 0
	}
	cut, gap = edges[len(edges)-1]*1.5, 0
	prev := 0.0
	for _, h := range edges {
		if h >= minHeight {
			if g := h - prev; g > gap {
				gap, cut = g, (prev+h)/2
			}
		}
		prev = h
	}
	return cut, gap
}

// --- Threshold sweep ------------------------------------------------------------------

// SweepPoint is one threshold and the clustering it produces.
type SweepPoint struct {
	Threshold float64
	Classes   int
	Sizes     []int
	Simpson   float64
}

// SweepThresholds clusters a subset at each threshold, reusing one MST pass.
func (h *HaploSample) SweepThresholds(subset []int, thresholds []float64) []SweepPoint {
	edges := h.mstEdges(subset)
	out := make([]SweepPoint, 0, len(thresholds))
	for _, t := range thresholds {
		c := clusterFromEdges(subset, edges, t)
		out = append(out, SweepPoint{
			Threshold: t,
			Classes:   len(c.Classes),
			Sizes:     c.Sizes(),
			Simpson:   c.Simpson(),
		})
	}
	return out
}

// DefaultThresholds is the sweep grid in fraction-of-sites units. It is deliberately
// dense between 0.10 and 0.20, where the calibration above showed both archived cases
// change their answer.
var DefaultThresholds = []float64{
	0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08,
	0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.18, 0.20, 0.22, 0.25,
}

// SiteCountThresholds is the sweep grid in differing-site counts, for Scale = 1.
//
// It is log-uniform rather than linear, and that is not cosmetic: Plateau measures how
// wide a range of thresholds agrees on a class count, so a grid whose spacing widens at
// the top would hand the prize to whatever the coarsest region says. Log spacing gives
// every decade equal weight, and the grid brackets the band the calibration identified
// on the archived n=40 exports, where 7 classes hold from 371 to 607 differing sites.
var SiteCountThresholds = []float64{
	1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 70, 100, 140, 180, 220, 260,
	300, 340, 380, 420, 460, 500, 550, 600, 700, 800, 1000, 1200, 1500, 2000,
}

// ThresholdGrid picks the grid matching the sample's current scale.
func (h *HaploSample) ThresholdGrid() []float64 {
	if h.Scale == 1 {
		return SiteCountThresholds
	}
	return DefaultThresholds
}

// Plateau returns the widest run of consecutive sweep points sharing a class count, and
// that count. The widest plateau is the most defensible single answer when a
// parameter-free cut is not wanted: it is the class count that the largest range of
// thresholds agrees on.
//
// Width is measured MULTIPLICATIVELY (hi/lo), not by subtraction. A divergence threshold
// has no natural zero and spans decades, so an additive width would make a plateau from
// 1,000 to 1,500 look wider than one from 20 to 200 — ranking the grid's coarse tail
// above a genuinely robust answer near the real class boundary.
func Plateau(sweep []SweepPoint) (classes int, lo, hi float64) {
	bestWidth := -1.0
	for i := 0; i < len(sweep); {
		j := i
		for j+1 < len(sweep) && sweep[j+1].Classes == sweep[i].Classes {
			j++
		}
		w := 0.0
		if sweep[i].Threshold > 0 {
			w = math.Log(sweep[j].Threshold / sweep[i].Threshold)
		}
		if w > bestWidth {
			bestWidth = w
			classes, lo, hi = sweep[i].Classes, sweep[i].Threshold, sweep[j].Threshold
		}
		i = j + 1
	}
	return classes, lo, hi
}
