package analysis

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

// ARGweaver output parser — TMR4A.md work item W7, the half that was deliberately
// left unwritten until a real ARGweaver output file existed to write it against.
//
// WHY THIS EXISTS AT ALL, when arg-summarize already reports a depth. Because
// `arg-summarize --tmrca` is the time to the FULL MRCA — TMR-1-A — and nothing in
// arg-summarize computes TMR-K-A for arbitrary K. The published statistic is K=4,
// and on a created-origin dataset the full MRCA frequently does not exist at all
// (DRIFT's own truth file leaves tmrca_year EMPTY for every ordinaryvar locus,
// because four founder copies seeded at t=0 never coalesce). So the only route to
// the statistic under study is to walk the per-sample local trees ourselves.
//
// THE INPUT is `arg-summarize -a <run>.bed.gz --tree`, whose columns are
//
//	chrom  chromStart  chromEnd  MCMC_sample  tree
//
// with `tree` a newick string carrying ARGweaver's NHX annotations, e.g.
// `...):22394.6[&&NHX:coal_time=7574.2],(...):26371.5[&&NHX:recomb_time=4975.8]);`
// The annotations describe the transition to the NEXT local tree and are not part
// of this tree's topology, so they are stripped. Branch lengths are in
// GENERATIONS and the trees are ultrametric (verified by hand on real output:
// 23.0+2607.4+1397.7+22394.6 = 51.2+26371.5 = 26422.7).
//
// THE JOIN IS LOCUS-CENTRIC, NOT INTERVAL-CENTRIC, and that is a deliberate
// choice. Every MCMC sample infers its OWN recombination breakpoints, so the
// intervals differ from sample to sample and cannot be lined up into a common
// grid without merging that would blur the very quantity being measured. DRIFT's
// truth is per focal locus, so for each focal locus we take, from every MCMC
// sample, whichever interval covers that locus. The result is a genuine posterior
// distribution of TMR-K-A at a point, directly comparable to the W5 walker's
// ground truth at the same point.
//
// TIES MATTER, and they are the same problem the DRIFT-side walker has.
// ARGweaver's discretised time grid puts many coalescences at identical heights —
// real output is full of `:0.0` branch lengths — so the lineage count can fall
// PAST K in a single step. TMR-K-A is therefore defined here exactly as
// pkg/analysis/tmrka.go:177 defines it: the first time the count is K OR BELOW.
// Whether it landed on K exactly or skipped past is recorded separately rather
// than silently averaged, because a locus that never had exactly K lineages is a
// different object from one that did, and §3.1's discipline about not collapsing
// outcome classes applies to the inferred side too.

// ARGweaverTreeStats is the inferred TMR-K-A at one focal locus, across MCMC samples.
type ARGweaverTreeStats struct {
	LocusBit int
	LocusBp  int

	Samples int // MCMC samples whose interval covered this locus

	TMRKAMean  float64 // generations
	TMRKAStdev float64
	TMRKAQ025  float64
	TMRKAQ50   float64
	TMRKAQ975  float64

	TMRCAMean float64 // root height, i.e. TMR-1-A, for cross-checking --tmrca

	ExactHits int // samples where the count landed on exactly K
	Skips     int // samples where it fell past K in one step

	// Truth, when a DRIFT tmrka_loci CSV was supplied.
	HasTruth       bool
	TruthOutcome   string
	TruthGens      float64 // -1 when censored (no TMR-K-A exists)
	TruthCreatedNA int     // created alleles surviving
}

// ARGweaverParseResult is the whole comparison.
type ARGweaverParseResult struct {
	K            int
	BpPerBit     int
	GenTime      float64
	Loci         []ARGweaverTreeStats
	RowsRead     int
	TreesParsed  int
	NumLeaves    int
	SpanWeighted float64 // span-weighted mean TMR-K-A over EVERY interval, generations
	SpanTotal    float64
	SpanTMRCA    float64 // same, for the root, to cross-check arg-summarize --tmrca

	// Burn-in accounting. MinSample/MaxSample echo the filter that was applied;
	// SampleLo/SampleHi are the range actually RETAINED, so a filter that silently
	// matched nothing is visible rather than being reported as a clean result.
	MinSample      int
	MaxSample      int
	RowsFiltered   int // rows dropped by the MCMC-sample filter
	SampleLo       int
	SampleHi       int
	sampleSeenOnce bool
}

// ---------------------------------------------------------------- newick

type nwNode struct {
	children []*nwNode
	blen     float64
	height   float64 // distance from the tips; leaves are 0
	isLeaf   bool
}

// parseNewick reads one newick string, discarding NHX comments and taxon names —
// only the topology and branch lengths bear on a lineage count. Returns the root.
func parseNewick(s string) (*nwNode, int, error) {
	p := &nwParser{s: s}
	root, err := p.parseNode()
	if err != nil {
		return nil, 0, err
	}
	leaves := 0
	var setHeights func(n *nwNode)
	setHeights = func(n *nwNode) {
		if n.isLeaf {
			n.height = 0
			leaves++
			return
		}
		h := 0.0
		for _, c := range n.children {
			setHeights(c)
			// Ultrametric in principle; take the max so a rounding wobble in the
			// input cannot make a parent shallower than one of its children.
			if ch := c.height + c.blen; ch > h {
				h = ch
			}
		}
		n.height = h
	}
	setHeights(root)
	return root, leaves, nil
}

type nwParser struct {
	s string
	i int
}

func (p *nwParser) parseNode() (*nwNode, error) {
	n := &nwNode{}
	p.skipSpace()
	if p.i < len(p.s) && p.s[p.i] == '(' {
		p.i++ // consume '('
		for {
			c, err := p.parseNode()
			if err != nil {
				return nil, err
			}
			n.children = append(n.children, c)
			p.skipSpace()
			if p.i >= len(p.s) {
				return nil, fmt.Errorf("newick: unexpected end inside clade")
			}
			switch p.s[p.i] {
			case ',':
				p.i++
				continue
			case ')':
				p.i++
			default:
				return nil, fmt.Errorf("newick: expected ',' or ')' at offset %d, got %q", p.i, p.s[p.i])
			}
			break
		}
	} else {
		n.isLeaf = true
	}
	// Taxon name (leaves) or internal label — read and discard.
	for p.i < len(p.s) && !strings.ContainsRune(":,()[;", rune(p.s[p.i])) {
		p.i++
	}
	p.skipComment()
	if p.i < len(p.s) && p.s[p.i] == ':' {
		p.i++
		start := p.i
		for p.i < len(p.s) && !strings.ContainsRune(",()[;", rune(p.s[p.i])) {
			p.i++
		}
		v, err := strconv.ParseFloat(strings.TrimSpace(p.s[start:p.i]), 64)
		if err != nil {
			return nil, fmt.Errorf("newick: bad branch length %q: %w", p.s[start:p.i], err)
		}
		n.blen = v
	}
	p.skipComment()
	return n, nil
}

// skipComment consumes an NHX block such as [&&NHX:coal_time=7574.2]. These
// annotate the transition to the next local tree, not this tree's topology.
func (p *nwParser) skipComment() {
	for p.i < len(p.s) && p.s[p.i] == '[' {
		depth := 0
		for p.i < len(p.s) {
			if p.s[p.i] == '[' {
				depth++
			} else if p.s[p.i] == ']' {
				depth--
				if depth == 0 {
					p.i++
					break
				}
			}
			p.i++
		}
	}
}

func (p *nwParser) skipSpace() {
	for p.i < len(p.s) && (p.s[p.i] == ' ' || p.s[p.i] == '\t') {
		p.i++
	}
}

// ---------------------------------------------------------------- the statistic

type coalEvent struct {
	height float64
	drops  int // lineages removed here: children - 1
}

// tmrKA walks the lineage count backwards and returns the first time it stands at
// K or below, mirroring pkg/analysis/tmrka.go:177 exactly. `exact` reports whether
// the count landed on K precisely; false means it fell past K in one step, which
// the discretised time grid makes common and which must not be reported as though
// K lineages had ever coexisted.
func tmrKA(root *nwNode, leaves, k int) (t float64, exact bool, ok bool) {
	if leaves <= k {
		return 0, leaves == k, true
	}
	var ev []coalEvent
	var walk func(n *nwNode)
	walk = func(n *nwNode) {
		if n.isLeaf {
			return
		}
		if d := len(n.children) - 1; d > 0 {
			ev = append(ev, coalEvent{height: n.height, drops: d})
		}
		for _, c := range n.children {
			walk(c)
		}
	}
	walk(root)
	sort.Slice(ev, func(i, j int) bool { return ev[i].height < ev[j].height })

	count := leaves
	for i := 0; i < len(ev); {
		h := ev[i].height
		// Apply every coalescence at this height together: within one discretised
		// time point they are simultaneous, exactly as one DRIFT generation is.
		for i < len(ev) && ev[i].height == h {
			count -= ev[i].drops
			i++
		}
		if count <= k {
			return h, count == k, true
		}
	}
	return 0, false, false
}

// ---------------------------------------------------------------- driving

// ARGweaverParseOptions configures ParseARGweaverTrees.
type ARGweaverParseOptions struct {
	TreesPath string // arg-summarize --tree output (plain or .gz)
	TruthPath string // DRIFT <model>_tmrka_loci_<run>_<year>.csv (optional)
	K         int
	BpPerBit  int
	GenTime   float64
	LocusBits []int // focal loci in BITS; taken from the truth file when empty

	// MinSample discards MCMC samples below this iteration — the burn-in filter.
	// It is not a refinement: pooling pre-burn-in samples is what produced the
	// bogus 1.48x on the nt60 run (Study Design §10.2), where the pooled median
	// was 35.6 gen while the post-burn-in value was >=50.1 and still climbing.
	// A chain whose recombination count is still falling has not converged, and
	// its shallow early trees drag the median down. Convergence is established
	// from the .log, NOT from this parser; this only excludes what that check
	// identified. Zero = keep everything.
	//
	// MaxSample discards samples above an iteration. Its purpose is the early/late
	// split: running twice, once capped low and once floored high, is how you show
	// the answer has stopped moving. Zero = no upper bound.
	MinSample int
	MaxSample int
}

// ParseARGweaverTrees turns arg-summarize --tree output into an inferred TMR-K-A
// per focal locus, joined against the W5 ground truth when a truth file is given.
func ParseARGweaverTrees(opts ARGweaverParseOptions) (*ARGweaverParseResult, error) {
	if opts.K < 1 {
		opts.K = 4
	}
	if opts.BpPerBit < 1 {
		opts.BpPerBit = 100
	}
	if opts.GenTime <= 0 {
		opts.GenTime = 25
	}

	res := &ARGweaverParseResult{
		K: opts.K, BpPerBit: opts.BpPerBit, GenTime: opts.GenTime,
		MinSample: opts.MinSample, MaxSample: opts.MaxSample,
	}

	truth := map[int]truthRow{}
	if opts.TruthPath != "" {
		var err error
		truth, err = readTruthLoci(opts.TruthPath)
		if err != nil {
			return nil, fmt.Errorf("reading truth file: %w", err)
		}
		if len(opts.LocusBits) == 0 {
			for b := range truth {
				opts.LocusBits = append(opts.LocusBits, b)
			}
		}
	}
	if len(opts.LocusBits) == 0 {
		return nil, fmt.Errorf("no focal loci: supply -argweaver-truth or -argweaver-loci")
	}
	sort.Ints(opts.LocusBits)

	// The exporter writes bp = (bit - contigBase)*bpPerBit + 1, 1-based
	// (pkg/analysis/argweaver.go:336). arg-summarize emits BED-style half-open
	// 0-based intervals, so compare against bp-1.
	type acc struct {
		vals   []float64
		tmrca  []float64
		exact  int
		skips  int
		bp     int
		bit    int
		offset int // 0-based coordinate used for interval containment
	}
	accs := make([]*acc, len(opts.LocusBits))
	for i, b := range opts.LocusBits {
		bp := b*opts.BpPerBit + 1
		accs[i] = &acc{bit: b, bp: bp, offset: bp - 1}
	}

	f, err := os.Open(opts.TreesPath)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var r io.Reader = f
	if strings.HasSuffix(strings.ToLower(opts.TreesPath), ".gz") {
		zr, err := gzip.NewReader(f)
		if err != nil {
			return nil, err
		}
		defer zr.Close()
		r = zr
	}

	sc := bufio.NewScanner(r)
	// Newick for 40 haplotypes with NHX annotation runs to several KB; the default
	// 64 KB token cap is not obviously enough on larger samples.
	sc.Buffer(make([]byte, 0, 1<<20), 64<<20)

	for sc.Scan() {
		line := sc.Text()
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		res.RowsRead++
		fields := strings.SplitN(line, "\t", 5)
		if len(fields) < 5 {
			// Tolerate space separation.
			fields = strings.Fields(line)
			if len(fields) < 5 {
				continue
			}
		}
		start, err1 := strconv.Atoi(strings.TrimSpace(fields[1]))
		end, err2 := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err1 != nil || err2 != nil {
			continue
		}

		// Burn-in filter, applied BEFORE the newick parse so excluded samples cost
		// nothing. Column 4 is arg-summarize's MCMC_sample. A row whose sample
		// column is unreadable is KEPT when no filter is set and DROPPED when one
		// is, because silently admitting unfilterable rows into a burned-in result
		// would defeat the point of asking for the filter.
		if opts.MinSample > 0 || opts.MaxSample > 0 {
			smp, err := strconv.Atoi(strings.TrimSpace(fields[3]))
			if err != nil ||
				(opts.MinSample > 0 && smp < opts.MinSample) ||
				(opts.MaxSample > 0 && smp > opts.MaxSample) {
				res.RowsFiltered++
				continue
			}
			if !res.sampleSeenOnce || smp < res.SampleLo {
				res.SampleLo = smp
			}
			if !res.sampleSeenOnce || smp > res.SampleHi {
				res.SampleHi = smp
			}
			res.sampleSeenOnce = true
		} else if smp, err := strconv.Atoi(strings.TrimSpace(fields[3])); err == nil {
			if !res.sampleSeenOnce || smp < res.SampleLo {
				res.SampleLo = smp
			}
			if !res.sampleSeenOnce || smp > res.SampleHi {
				res.SampleHi = smp
			}
			res.sampleSeenOnce = true
		}

		// Every tree is parsed, not just the ones covering a focal locus: the
		// span-weighted figure over ALL intervals is what cross-checks this parser
		// against arg-summarize --tmrca, and it is the only guard against a subtly
		// wrong newick reading that would still look plausible locus by locus.
		needed := false
		for _, a := range accs {
			if a.offset >= start && a.offset < end {
				needed = true
				break
			}
		}

		root, leaves, err := parseNewick(strings.TrimSpace(fields[4]))
		if err != nil {
			return nil, fmt.Errorf("line %d: %w", res.RowsRead, err)
		}
		res.TreesParsed++
		if res.NumLeaves == 0 {
			res.NumLeaves = leaves
		}

		t, exact, ok := tmrKA(root, leaves, opts.K)
		w := float64(end - start)
		if w <= 0 {
			w = 1
		}
		if ok {
			res.SpanWeighted += w * t
			res.SpanTotal += w
			res.SpanTMRCA += w * root.height
		}

		if !needed {
			continue
		}
		for _, a := range accs {
			if a.offset >= start && a.offset < end && ok {
				a.vals = append(a.vals, t)
				a.tmrca = append(a.tmrca, root.height)
				if exact {
					a.exact++
				} else {
					a.skips++
				}
			}
		}
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	if res.SpanTotal > 0 {
		res.SpanWeighted /= res.SpanTotal
		res.SpanTMRCA /= res.SpanTotal
	}

	for _, a := range accs {
		st := ARGweaverTreeStats{LocusBit: a.bit, LocusBp: a.bp, Samples: len(a.vals),
			ExactHits: a.exact, Skips: a.skips}
		if len(a.vals) > 0 {
			st.TMRKAMean, st.TMRKAStdev = meanStdev(a.vals)
			sorted := append([]float64(nil), a.vals...)
			sort.Float64s(sorted)
			st.TMRKAQ025 = quantile(sorted, 0.025)
			st.TMRKAQ50 = quantile(sorted, 0.5)
			st.TMRKAQ975 = quantile(sorted, 0.975)
			st.TMRCAMean, _ = meanStdev(a.tmrca)
		}
		if tr, ok := truth[a.bit]; ok {
			st.HasTruth = true
			st.TruthOutcome = tr.outcome
			st.TruthGens = tr.gens
			st.TruthCreatedNA = tr.createdAlleles
		}
		res.Loci = append(res.Loci, st)
	}
	return res, nil
}

type truthRow struct {
	outcome        string
	gens           float64 // -1 when censored
	createdAlleles int
}

func readTruthLoci(path string) (map[int]truthRow, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	rd := csv.NewReader(f)
	rd.FieldsPerRecord = -1
	recs, err := rd.ReadAll()
	if err != nil {
		return nil, err
	}
	if len(recs) < 2 {
		return nil, fmt.Errorf("%s: no data rows", path)
	}
	col := map[string]int{}
	for i, h := range recs[0] {
		col[strings.TrimSpace(h)] = i
	}
	need := []string{"locus", "outcome"}
	for _, n := range need {
		if _, ok := col[n]; !ok {
			return nil, fmt.Errorf("%s: missing column %q", path, n)
		}
	}
	out := map[int]truthRow{}
	for _, r := range recs[1:] {
		if len(r) <= col["outcome"] {
			continue
		}
		bit, err := strconv.Atoi(strings.TrimSpace(r[col["locus"]]))
		if err != nil {
			continue
		}
		tr := truthRow{outcome: strings.TrimSpace(r[col["outcome"]]), gens: -1, createdAlleles: -1}
		if i, ok := col["tmrka_gens_ago"]; ok && i < len(r) {
			if v, err := strconv.ParseFloat(strings.TrimSpace(r[i]), 64); err == nil {
				tr.gens = v
			}
		}
		if i, ok := col["created_alleles_surviving"]; ok && i < len(r) {
			if v, err := strconv.Atoi(strings.TrimSpace(r[i])); err == nil {
				tr.createdAlleles = v
			}
		}
		out[bit] = tr
	}
	return out, nil
}

func meanStdev(v []float64) (float64, float64) {
	if len(v) == 0 {
		return 0, 0
	}
	s := 0.0
	for _, x := range v {
		s += x
	}
	m := s / float64(len(v))
	if len(v) < 2 {
		return m, 0
	}
	ss := 0.0
	for _, x := range v {
		ss += (x - m) * (x - m)
	}
	return m, math.Sqrt(ss / float64(len(v)-1))
}

func quantile(sorted []float64, q float64) float64 {
	if len(sorted) == 0 {
		return 0
	}
	if len(sorted) == 1 {
		return sorted[0]
	}
	p := q * float64(len(sorted)-1)
	lo := int(math.Floor(p))
	hi := int(math.Ceil(p))
	if lo == hi {
		return sorted[lo]
	}
	return sorted[lo] + (p-float64(lo))*(sorted[hi]-sorted[lo])
}

// ---------------------------------------------------------------- reporting

// WriteARGweaverComparison writes the truth-vs-inference table — the paper's
// central object, with the created-allele floor beside the inferred depth.
func WriteARGweaverComparison(path string, res *ARGweaverParseResult) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	defer w.Flush()

	hdr := []string{"locus_bit", "locus_bp", "mcmc_samples",
		fmt.Sprintf("inferred_tmr%da_mean_gen", res.K), "inferred_stdev_gen",
		"inferred_q025", "inferred_q50", "inferred_q975",
		"inferred_tmrca_mean_gen", "exact_hits", "skipped_past_k",
		"truth_outcome", "truth_gens_ago", "truth_created_alleles",
		"inferred_over_truth", "inferred_years"}
	if err := w.Write(hdr); err != nil {
		return err
	}

	for _, l := range res.Loci {
		ratio := ""
		if l.HasTruth && l.TruthGens > 0 && l.Samples > 0 {
			ratio = fmt.Sprintf("%.2f", l.TMRKAMean/l.TruthGens)
		}
		outcome, tg, ca := "", "", ""
		if l.HasTruth {
			outcome = l.TruthOutcome
			if l.TruthGens >= 0 {
				tg = fmt.Sprintf("%.2f", l.TruthGens)
			}
			if l.TruthCreatedNA >= 0 {
				ca = strconv.Itoa(l.TruthCreatedNA)
			}
		}
		rec := []string{
			strconv.Itoa(l.LocusBit), strconv.Itoa(l.LocusBp), strconv.Itoa(l.Samples),
			fmt.Sprintf("%.1f", l.TMRKAMean), fmt.Sprintf("%.1f", l.TMRKAStdev),
			fmt.Sprintf("%.1f", l.TMRKAQ025), fmt.Sprintf("%.1f", l.TMRKAQ50),
			fmt.Sprintf("%.1f", l.TMRKAQ975),
			fmt.Sprintf("%.1f", l.TMRCAMean),
			strconv.Itoa(l.ExactHits), strconv.Itoa(l.Skips),
			outcome, tg, ca, ratio,
			fmt.Sprintf("%.0f", l.TMRKAMean*res.GenTime),
		}
		if err := w.Write(rec); err != nil {
			return err
		}
	}
	return nil
}

// PrintARGweaverComparison prints the console summary. The three DRIFT outcome
// classes are kept distinct all the way through, per §3.1 — a censored locus has
// NO true TMR-K-A, and reporting a ratio against one would be the exact error the
// module exists to expose.
func PrintARGweaverComparison(res *ARGweaverParseResult) {
	fmt.Printf("\nARGweaver TMR-%d-A (W7 parser): %d rows, %d trees parsed, %d haplotypes\n",
		res.K, res.RowsRead, res.TreesParsed, res.NumLeaves)
	if res.MinSample > 0 || res.MaxSample > 0 {
		fmt.Printf("  MCMC-sample filter: ")
		switch {
		case res.MinSample > 0 && res.MaxSample > 0:
			fmt.Printf("keeping %d-%d", res.MinSample, res.MaxSample)
		case res.MinSample > 0:
			fmt.Printf("burn-in, discarding samples below %d", res.MinSample)
		default:
			fmt.Printf("discarding samples above %d", res.MaxSample)
		}
		fmt.Printf(" — %d of %d rows dropped, %d retained (samples %d-%d)\n",
			res.RowsFiltered, res.RowsRead, res.RowsRead-res.RowsFiltered,
			res.SampleLo, res.SampleHi)
		if res.RowsRead-res.RowsFiltered == 0 {
			fmt.Printf("  WARNING: the filter matched NOTHING. Every row was discarded, so the\n")
			fmt.Printf("    figures below are empty rather than burned-in. Check the sample numbering.\n")
		}
	} else if res.sampleSeenOnce {
		fmt.Printf("  MCMC samples present: %d-%d (NO burn-in filter — pass -argweaver-min-sample\n"+
			"    if the chain's early iterations are not converged)\n", res.SampleLo, res.SampleHi)
	}
	fmt.Printf("  span-weighted over every interval: TMR-%d-A %.0f gen (%.0f y), TMRCA %.0f gen (%.0f y)\n",
		res.K, res.SpanWeighted, res.SpanWeighted*res.GenTime,
		res.SpanTMRCA, res.SpanTMRCA*res.GenTime)

	var coalesced, censored, missing int
	var ratioSum float64
	var ratioN int
	var skipTot, exactTot int
	for _, l := range res.Loci {
		skipTot += l.Skips
		exactTot += l.ExactHits
		if l.Samples == 0 {
			missing++
			continue
		}
		switch {
		case !l.HasTruth:
		case l.TruthOutcome == OutcomeCoalesced && l.TruthGens > 0:
			coalesced++
			ratioSum += l.TMRKAMean / l.TruthGens
			ratioN++
		case l.TruthOutcome == OutcomeCensored:
			censored++
		}
	}
	fmt.Printf("  focal loci covered: %d of %d", len(res.Loci)-missing, len(res.Loci))
	if missing > 0 {
		fmt.Printf("  (%d outside the sampled region — widen --region or move arg_loci)", missing)
	}
	fmt.Println()
	if exactTot+skipTot > 0 {
		fmt.Printf("  lineage count landed on exactly K in %d of %d locus-samples; skipped past K in %d\n",
			exactTot, exactTot+skipTot, skipTot)
	}
	if ratioN > 0 {
		fmt.Printf("  where DRIFT truth COALESCED (%d loci): mean inferred/true = %.1fx\n",
			coalesced, ratioSum/float64(ratioN))
	}
	if censored > 0 {
		fmt.Printf("  where DRIFT truth is CENSORED AT FOUNDING (%d loci): no true TMR-%d-A exists.\n",
			censored, res.K)
		fmt.Printf("    ARGweaver reports a depth anyway; that gap is the finding, and no ratio is\n")
		fmt.Printf("    computed for these loci because dividing by a nonexistent date would be the\n")
		fmt.Printf("    exact error this module exists to expose (TMR4A.md §0).\n")
	}
}
