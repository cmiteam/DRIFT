package analysis

// Haplotype & selection statistics — roadmap §6g. EHH, iHS, XP-EHH and cM-binned
// LD-decay: the standard selection-scan / haplotype-structure toolkit (Sabeti
// 2002; Voight 2006; Sabeti 2007). Serves goal 2 (testing selection claims tied
// to the OoA expansion) and goal 3 (comparability to the field's own tools —
// selscan / rehh conventions).
//
// SUBSTRATE = THE DE-NOVO MUTATION POOL (decided with Rob, from an empirical
// probe of the real engine). DRIFT has two genetic-bookkeeping systems and only
// one is a usable EHH substrate:
//
//   * the Chromosomes BITFIELD holds only FOUNDER standing variation. It receives
//     no de-novo input, drifts to loss/fixation (§6h), and — because the seeder
//     assigns INDEPENDENT per-site alleles — carries essentially no haplotype LD.
//     A probe found ~5-10 segregating sites per chromosome and near-flat EHH.
//     Critically, a beneficial de-novo mutation NEVER enters the bitfield, so a
//     selective sweep is structurally invisible there.
//
//   * the MUTATION POOL (pop.IndMutations + pop.MutationPool, gated by
//     track_mutations) receives continuous Poisson(mu) input; each mutation arises
//     on ONE ancestral strand and is inherited with its neighbours until
//     recombination separates them — genuine coalescent haplotype blocks, exactly
//     what EHH measures. The probe found ~1000 segregating sites per chromosome,
//     clean cM-scaled EHH decay, and — with selection cranked above the §6f drift
//     barrier (|s| >~ 1/(2*Ne), Ne/N ~= 0.31) — a ~4x extended-haplotype signature
//     on a swept haplotype vs a frequency-matched neutral core. This is where
//     selection acts, so this is what §6g reads.
//
// A pool "site" is one infinite-sites mutation lineage: a strand's DERIVED allele
// at that site = it carries the mutation id; ANCESTRAL = it does not. Two strands
// are identical over an interval iff they carry the same set of mutation ids in
// it. EHH is computed WITHIN a chromosome (recombination and the cM map are
// per-chromosome); distances are real centiMorgans from the §6a GeneticMap
// accessor (utils.EnsureGeneticMap), falling back to a uniform recomb_cM_per_bit
// map when no recombination map is configured.
//
// WIRING. track_haplostats (default off) gates the whole thing; a read-only
// end-of-run analysis that draws NO RNG on the full-sample path (haplostats_
// sample_size > 0 subsamples and does draw, opt-in — the §6d/§6f pattern), so a
// track_haplostats run is byte-identical to a track_haplostats-off run. Requires
// track_mutations (the pool). Output is shaped to selscan/rehh column conventions
// and, when haplostats_export_hap is set, the sampled pool haplotypes are written
// as selscan .hap/.map files so a critic can reproduce the numbers with their own
// tools.

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"math"
	"sort"
	"strconv"
	"strings"
)

// ---- tunables (all overridable via params; see HaploStatsConfig) ----

// HaploStatsConfig collects the knobs read from model params so the compute path
// is testable without a full Model.
type HaploStatsConfig struct {
	MinMAF     float64 // cores must have minor-allele freq >= this (default 0.05)
	EHHCutoff  float64 // stop integrating a side once EHH drops below this (default 0.05)
	NumFreqBin int     // derived-freq bins for iHS standardization (default 20)
	MaxCurves  int     // number of top-|iHS| cores to emit full EHH curves for (default 10)
	DemeA      int     // XP-EHH deme A (-1 = auto/none)
	DemeB      int     // XP-EHH deme B (-1 = auto/none)
}

// DefaultHaploConfig returns the shipped defaults.
func DefaultHaploConfig() HaploStatsConfig {
	return HaploStatsConfig{MinMAF: 0.05, EHHCutoff: 0.05, NumFreqBin: 20, MaxCurves: 10, DemeA: -1, DemeB: -1}
}

// ParseDemePair parses a "A,B" (or "A;B" / "A B") deme spec for XP-EHH into two
// deme ids. An empty/malformed spec returns (-1,-1) = auto-pick the first two
// demes present. Matches the separator-agnostic style of the other §6 params.
func ParseDemePair(spec string) (int, int) {
	f := func(r rune) bool { return r == ',' || r == ';' || r == ' ' || r == '\t' }
	toks := strings.FieldsFunc(spec, f)
	if len(toks) < 2 {
		return -1, -1
	}
	a, err1 := strconv.Atoi(toks[0])
	b, err2 := strconv.Atoi(toks[1])
	if err1 != nil || err2 != nil {
		return -1, -1
	}
	return a, b
}

// ---- result types ----

// CoreEHH is the per-core iHH / iHS record (the iHS scan headline). iHH is the
// integrated EHH (both flanks) in EHH*cM units. UnstdIHS = ln(iHHAnc/iHHDer); a
// derived-allele sweep drives iHHDer up, so UnstdIHS goes NEGATIVE. StdIHS is
// UnstdIHS standardized to mean 0 / sd 1 within its derived-frequency bin (the
// selscan/rehh normalization; NaN when its bin has too few cores).
type CoreEHH struct {
	MutID        int
	Chrom        int
	Position     int
	DerivedCount int
	Nhap         int
	DerivedFreq  float64
	IHHDerived   float64
	IHHAncestral float64
	UnstdIHS     float64
	StdIHS       float64
}

// EHHCurvePoint is one row of an EHH-vs-distance decay curve for a single core.
type EHHCurvePoint struct {
	MutID        int
	Chrom        int
	Position     int
	CM           float64 // signed cM offset from the core (negative = upstream)
	EHHDerived   float64
	EHHAncestral float64
}

// XPEHHSite is the cross-population unstandardized/standardized XP-EHH at a site
// segregating in the pooled two-deme sample. UnstdXPEHH = ln(iESA/iESB); StdXPEHH
// standardizes it genome-wide (mean 0 / sd 1).
type XPEHHSite struct {
	MutID      int
	Chrom      int
	Position   int
	IESA       float64
	IESB       float64
	UnstdXPEHH float64
	StdXPEHH   float64
}

// LDBin is one cM-distance bin of the LD-decay curve (mean r^2 over sampled
// same-chromosome site pairs whose cM separation falls in the bin).
type LDBin struct {
	LoCM     float64
	HiCM     float64
	MeanR2   float64
	NumPairs int
}

// HaploStatsResult is the full §6g product for one sample.
type HaploStatsResult struct {
	NumSamples int
	NumHaps    int
	NumCores   int
	Cores      []CoreEHH
	Curves     []EHHCurvePoint
	XPEHH      []XPEHHSite
	DemeA      int // -1 when XP-EHH not computed
	DemeB      int
	LDDecay    []LDBin
	Config     HaploStatsConfig
}

// ---- per-chromosome haplotype workspace ----

type hsSite struct {
	mutID   int
	pos     int
	cumCM   float64
	derived int // derived-allele count over all strands in the workspace
}

type chromHaps struct {
	chrom  int
	sites  []hsSite       // sorted by pos
	idxOf  map[int]int    // mutID -> index in sites
	member []map[int]bool // per strand: set of carried mutIDs on this chromosome
	deme   []int          // per strand: deme id
	nhap   int            // len(member)
}

// buildChromHaps extracts, per chromosome, the segregating pool sites and each
// sampled strand's carried-mutation set. gm supplies cumulative cM per site.
func buildChromHaps(pop *core.Pop, ids []int, arms map[int][2]int, gm *core.GeneticMap) map[int]*chromHaps {
	nhap := 2 * len(ids)
	// strand index s = 2*sampleIdx + copy.
	demeOf := make([]int, nhap)
	for i, id := range ids {
		d := 0
		if row, ok := pop.IndData[id]; ok && len(row) > int(individual.Deme) {
			d = row[individual.Deme]
		}
		demeOf[2*i], demeOf[2*i+1] = d, d
	}

	// First pass: derived counts per mutID (over the sample), and each strand's
	// carried set keyed to chromosome via the mutation's position.
	type acc struct {
		pos, chrom, derived int
	}
	info := map[int]*acc{}
	carried := make([]map[int]int, nhap) // strand -> mutID -> chrom (transient)
	for i, id := range ids {
		strands, ok := pop.IndMutations[id]
		if !ok {
			continue
		}
		for c := 0; c < 2; c++ {
			h := 2*i + c
			seen := map[int]bool{}
			for _, mid := range strands[c] {
				if seen[mid] {
					continue // a strand carries a lineage once (dedupe accidental repeats)
				}
				seen[mid] = true
				mut, ok := pop.MutationPool[mid]
				if !ok {
					continue
				}
				ch := chromOf(mut.Position, arms)
				if ch < 0 {
					continue
				}
				a := info[mid]
				if a == nil {
					a = &acc{pos: mut.Position, chrom: ch}
					info[mid] = a
				}
				a.derived++
				if carried[h] == nil {
					carried[h] = map[int]int{}
				}
				carried[h][mid] = ch
			}
		}
	}

	out := map[int]*chromHaps{}
	for ch, r := range arms {
		cw := &chromHaps{chrom: ch, idxOf: map[int]int{}, deme: demeOf, nhap: nhap}
		cw.member = make([]map[int]bool, nhap)
		for h := 0; h < nhap; h++ {
			cw.member[h] = map[int]bool{}
		}
		out[ch] = cw
		_ = r
	}
	// Populate sites (segregating only) and membership.
	for mid, a := range info {
		if a.derived <= 0 || a.derived >= nhap {
			continue // monomorphic in the sample
		}
		cw := out[a.chrom]
		if cw == nil {
			continue
		}
		cm := 0.0
		if gm != nil {
			cm = gm.CumCM(a.chrom, a.pos)
		}
		cw.sites = append(cw.sites, hsSite{mutID: mid, pos: a.pos, cumCM: cm, derived: a.derived})
	}
	for h := 0; h < nhap; h++ {
		for mid, ch := range carried[h] {
			if cw := out[ch]; cw != nil {
				if _, ok := info[mid]; ok && info[mid].derived < nhap {
					cw.member[h][mid] = true
				}
			}
		}
	}
	for _, cw := range out {
		sort.Slice(cw.sites, func(i, j int) bool {
			if cw.sites[i].pos != cw.sites[j].pos {
				return cw.sites[i].pos < cw.sites[j].pos
			}
			return cw.sites[i].mutID < cw.sites[j].mutID
		})
		for i, s := range cw.sites {
			cw.idxOf[s.mutID] = i
		}
	}
	return out
}

func chromOf(pos int, arms map[int][2]int) int {
	for ch, r := range arms {
		if pos >= r[0] && pos < r[1] {
			return ch
		}
	}
	return -1
}

// armRanges maps chromosome -> [minStart, maxEnd) over its arms.
func armRangesOf(model *core.Model) map[int][2]int {
	r := map[int][2]int{}
	for c, arms := range model.ChromosomeArms {
		lo, hi := 1<<30, 0
		for _, a := range arms {
			if len(a) < 2 {
				continue
			}
			if a[0] < lo {
				lo = a[0]
			}
			if a[0]+a[1] > hi {
				hi = a[0] + a[1]
			}
		}
		if hi > lo {
			r[c] = [2]int{lo, hi}
		}
	}
	return r
}

// ---- EHH core walk ----

// ehhSide integrates EHH away from core site index ci in direction dir (+1 =
// increasing position, -1 = decreasing) over the strands in `class`, using the
// standard incremental haplotype-partition split. Integration (trapezoid in cM)
// stops once EHH drops below cutoff or the chromosome end is reached. Returns the
// integrated iHH (EHH*cM). If curve != nil, appends (signed-cM, ehh) points via
// the supplied callback.
func ehhSide(cw *chromHaps, class []int, ci, dir int, cutoff float64, emit func(cm, ehh float64)) float64 {
	n := len(class)
	if n < 2 {
		return 0
	}
	typeOf := make([]int, n)
	counts := map[int]int{0: n}
	coreCM := cw.sites[ci].cumCM
	prevCM := coreCM
	prevEHH := 1.0
	if emit != nil {
		emit(0, 1.0)
	}
	iHH := 0.0
	for k := ci + dir; k >= 0 && k < len(cw.sites); k += dir {
		site := cw.sites[k]
		split := map[[2]int]int{}
		next := 0
		newType := make([]int, n)
		newCounts := map[int]int{}
		for j, h := range class {
			bit := 0
			if cw.member[h][site.mutID] {
				bit = 1
			}
			key := [2]int{typeOf[j], bit}
			nt, ok := split[key]
			if !ok {
				nt = next
				next++
				split[key] = nt
			}
			newType[j] = nt
			newCounts[nt]++
		}
		typeOf, counts = newType, newCounts
		ehh := homozygosity(counts, n)
		cm := site.cumCM
		dcm := math.Abs(cm - prevCM)
		iHH += (ehh + prevEHH) / 2 * dcm
		if emit != nil {
			emit(cm-coreCM, ehh)
		}
		prevCM, prevEHH = cm, ehh
		if ehh < cutoff {
			break
		}
	}
	return iHH
}

func homozygosity(counts map[int]int, n int) float64 {
	if n < 2 {
		return math.NaN()
	}
	num := 0
	for _, c := range counts {
		num += c * (c - 1)
	}
	return float64(num) / float64(n*(n-1))
}

// ---- top-level compute ----

// ComputeHaploStats runs the full §6g analysis over the sampled individuals.
// Returns nil when there is no map/genome or fewer than two samples.
func ComputeHaploStats(model *core.Model, pop *core.Pop, ids []int, cfg HaploStatsConfig) *HaploStatsResult {
	if len(ids) < 2 {
		return nil
	}
	arms := armRangesOf(model)
	if len(arms) == 0 {
		return nil
	}
	gm := haploGeneticMap(model)
	cw := buildChromHaps(pop, ids, arms, gm)

	res := &HaploStatsResult{
		NumSamples: len(ids),
		NumHaps:    2 * len(ids),
		DemeA:      -1,
		DemeB:      -1,
		Config:     cfg,
	}

	// --- iHS scan over all chromosomes ---
	chroms := sortedChromKeys(cw)
	for _, ch := range chroms {
		c := cw[ch]
		nhap := c.nhap
		// derived / ancestral strand index lists per core computed on the fly.
		for ci, site := range c.sites {
			der := make([]int, 0, site.derived)
			anc := make([]int, 0, nhap-site.derived)
			for h := 0; h < nhap; h++ {
				if c.member[h][site.mutID] {
					der = append(der, h)
				} else {
					anc = append(anc, h)
				}
			}
			freq := float64(site.derived) / float64(nhap)
			maf := math.Min(freq, 1-freq)
			if maf < cfg.MinMAF {
				continue
			}
			if len(der) < 2 || len(anc) < 2 {
				continue
			}
			iD := ehhSide(c, der, ci, +1, cfg.EHHCutoff, nil) + ehhSide(c, der, ci, -1, cfg.EHHCutoff, nil)
			iA := ehhSide(c, anc, ci, +1, cfg.EHHCutoff, nil) + ehhSide(c, anc, ci, -1, cfg.EHHCutoff, nil)
			if iD <= 0 || iA <= 0 {
				continue
			}
			res.Cores = append(res.Cores, CoreEHH{
				MutID: site.mutID, Chrom: ch, Position: site.pos,
				DerivedCount: site.derived, Nhap: nhap, DerivedFreq: freq,
				IHHDerived: iD, IHHAncestral: iA,
				UnstdIHS: math.Log(iA / iD),
			})
		}
	}
	res.NumCores = len(res.Cores)
	standardizeIHS(res.Cores, cfg.NumFreqBin)

	// --- full EHH decay curves for the top-|iHS| cores ---
	res.Curves = topCurves(cw, res.Cores, cfg)

	// --- XP-EHH (if a two-deme split is available) ---
	res.XPEHH, res.DemeA, res.DemeB = computeXPEHH(cw, chroms, cfg)

	// --- LD-decay in cM bins ---
	res.LDDecay = computeLDDecay(cw, chroms)

	return res
}

// haploGeneticMap returns the model's genetic map, or builds a uniform one from
// recomb_cM_per_bit (default 0.05 cM/bit) when none is configured, so EHH always
// has a cM axis. The locally built map is NOT stored on the model (no side effect
// on the recombination kernel).
func haploGeneticMap(model *core.Model) *core.GeneticMap {
	if gm := utils.EnsureGeneticMap(model); gm != nil {
		return gm
	}
	rate := model.Parameters["recomb_cM_per_bit"]
	if rate <= 0 {
		rate = 0.05
	}
	return core.BuildGeneticMap(model.ChromosomeArms, model.ArmCM, rate)
}

func sortedChromKeys(cw map[int]*chromHaps) []int {
	ks := make([]int, 0, len(cw))
	for k := range cw {
		ks = append(ks, k)
	}
	sort.Ints(ks)
	return ks
}

// standardizeIHS fills StdIHS: within each derived-frequency bin, subtract the
// bin mean and divide by the bin sd (selscan/rehh normalization). Bins with < 2
// finite cores leave StdIHS = NaN.
func standardizeIHS(cores []CoreEHH, nbin int) {
	if nbin < 1 {
		nbin = 1
	}
	type stat struct {
		sum, sumsq float64
		n          int
	}
	bins := make([]stat, nbin)
	binOf := func(f float64) int {
		b := int(f * float64(nbin))
		if b < 0 {
			b = 0
		}
		if b >= nbin {
			b = nbin - 1
		}
		return b
	}
	for i := range cores {
		if math.IsNaN(cores[i].UnstdIHS) || math.IsInf(cores[i].UnstdIHS, 0) {
			continue
		}
		b := binOf(cores[i].DerivedFreq)
		bins[b].sum += cores[i].UnstdIHS
		bins[b].sumsq += cores[i].UnstdIHS * cores[i].UnstdIHS
		bins[b].n++
	}
	means := make([]float64, nbin)
	sds := make([]float64, nbin)
	for b := range bins {
		if bins[b].n >= 2 {
			m := bins[b].sum / float64(bins[b].n)
			v := bins[b].sumsq/float64(bins[b].n) - m*m
			means[b] = m
			if v > 0 {
				sds[b] = math.Sqrt(v)
			}
		}
	}
	for i := range cores {
		cores[i].StdIHS = math.NaN()
		if math.IsNaN(cores[i].UnstdIHS) || math.IsInf(cores[i].UnstdIHS, 0) {
			continue
		}
		b := binOf(cores[i].DerivedFreq)
		if bins[b].n >= 2 && sds[b] > 0 {
			cores[i].StdIHS = (cores[i].UnstdIHS - means[b]) / sds[b]
		}
	}
}

// topCurves emits full EHH decay curves for the MaxCurves cores with the largest
// |StdIHS| (falling back to |UnstdIHS| when unstandardized), sorted for stable
// output.
func topCurves(cw map[int]*chromHaps, cores []CoreEHH, cfg HaploStatsConfig) []EHHCurvePoint {
	if cfg.MaxCurves <= 0 || len(cores) == 0 {
		return nil
	}
	ranked := make([]CoreEHH, len(cores))
	copy(ranked, cores)
	score := func(c CoreEHH) float64 {
		if !math.IsNaN(c.StdIHS) {
			return math.Abs(c.StdIHS)
		}
		return math.Abs(c.UnstdIHS)
	}
	sort.Slice(ranked, func(i, j int) bool {
		si, sj := score(ranked[i]), score(ranked[j])
		if si != sj {
			return si > sj
		}
		return ranked[i].MutID < ranked[j].MutID
	})
	if len(ranked) > cfg.MaxCurves {
		ranked = ranked[:cfg.MaxCurves]
	}
	var pts []EHHCurvePoint
	for _, core0 := range ranked {
		c := cw[core0.Chrom]
		ci, ok := c.idxOf[core0.MutID]
		if !ok {
			continue
		}
		der := make([]int, 0)
		anc := make([]int, 0)
		for h := 0; h < c.nhap; h++ {
			if c.member[h][core0.MutID] {
				der = append(der, h)
			} else {
				anc = append(anc, h)
			}
		}
		curve := map[float64][2]float64{} // cm -> {ehhD, ehhA}
		var order []float64
		seenCM := map[float64]bool{}
		add := func(idx int, cm, ehh float64) {
			v := curve[cm]
			v[idx] = ehh
			curve[cm] = v
			if !seenCM[cm] {
				seenCM[cm] = true
				order = append(order, cm)
			}
		}
		emitD := func(cm, ehh float64) { add(0, cm, ehh) }
		emitA := func(cm, ehh float64) { add(1, cm, ehh) }
		ehhSide(c, der, ci, +1, cfg.EHHCutoff, emitD)
		ehhSide(c, der, ci, -1, cfg.EHHCutoff, emitD)
		ehhSide(c, anc, ci, +1, cfg.EHHCutoff, emitA)
		ehhSide(c, anc, ci, -1, cfg.EHHCutoff, emitA)
		sort.Float64s(order)
		for _, cm := range order {
			v := curve[cm]
			pts = append(pts, EHHCurvePoint{
				MutID: core0.MutID, Chrom: core0.Chrom, Position: core0.Position,
				CM: cm, EHHDerived: v[0], EHHAncestral: v[1],
			})
		}
	}
	return pts
}

// computeXPEHH computes cross-population XP-EHH between two demes. iES (site EHH,
// unpolarized) is integrated per deme over each site segregating in the pooled
// two-deme sample; XP-EHH = ln(iESA/iESB), standardized genome-wide. Returns
// (nil,-1,-1) when fewer than two demes are represented.
func computeXPEHH(cw map[int]*chromHaps, chroms []int, cfg HaploStatsConfig) ([]XPEHHSite, int, int) {
	demeA, demeB := cfg.DemeA, cfg.DemeB
	if demeA < 0 || demeB < 0 {
		present := map[int]bool{}
		for _, ch := range chroms {
			for _, d := range cw[ch].deme {
				present[d] = true
			}
			break
		}
		ds := make([]int, 0, len(present))
		for d := range present {
			ds = append(ds, d)
		}
		sort.Ints(ds)
		if len(ds) < 2 {
			return nil, -1, -1
		}
		demeA, demeB = ds[0], ds[1]
	}
	var out []XPEHHSite
	for _, ch := range chroms {
		c := cw[ch]
		aIdx := demeStrands(c, demeA)
		bIdx := demeStrands(c, demeB)
		if len(aIdx) < 2 || len(bIdx) < 2 {
			continue
		}
		for ci, site := range c.sites {
			// segregating in the pooled A+B sample?
			dA := countCarry(c, aIdx, site.mutID)
			dB := countCarry(c, bIdx, site.mutID)
			tot := dA + dB
			if tot == 0 || tot == len(aIdx)+len(bIdx) {
				continue
			}
			iesA := iesSite(c, aIdx, ci, cfg.EHHCutoff)
			iesB := iesSite(c, bIdx, ci, cfg.EHHCutoff)
			if iesA <= 0 || iesB <= 0 {
				continue
			}
			out = append(out, XPEHHSite{
				MutID: site.mutID, Chrom: ch, Position: site.pos,
				IESA: iesA, IESB: iesB, UnstdXPEHH: math.Log(iesA / iesB),
			})
		}
	}
	standardizeXPEHH(out)
	return out, demeA, demeB
}

// iesSite integrates unpolarized site-EHH (EHHS) both flanks for a restricted set
// of strands. It reuses ehhSide over a local workspace whose strand list is the
// given deme's strands (the walk partitions all of them, ignoring the core
// allele).
func iesSite(c *chromHaps, idx []int, ci int, cutoff float64) float64 {
	// ehhSide walks over c.member[h] for h in class; iES uses the whole deme as
	// one class so the partition homozygosity is the site EHHS.
	return ehhSide(c, idx, ci, +1, cutoff, nil) + ehhSide(c, idx, ci, -1, cutoff, nil)
}

func demeStrands(c *chromHaps, deme int) []int {
	var out []int
	for h := 0; h < c.nhap; h++ {
		if c.deme[h] == deme {
			out = append(out, h)
		}
	}
	return out
}

func countCarry(c *chromHaps, idx []int, mutID int) int {
	n := 0
	for _, h := range idx {
		if c.member[h][mutID] {
			n++
		}
	}
	return n
}

func standardizeXPEHH(sites []XPEHHSite) {
	var sum, sumsq float64
	n := 0
	for _, s := range sites {
		if math.IsNaN(s.UnstdXPEHH) || math.IsInf(s.UnstdXPEHH, 0) {
			continue
		}
		sum += s.UnstdXPEHH
		sumsq += s.UnstdXPEHH * s.UnstdXPEHH
		n++
	}
	if n < 2 {
		for i := range sites {
			sites[i].StdXPEHH = math.NaN()
		}
		return
	}
	m := sum / float64(n)
	v := sumsq/float64(n) - m*m
	sd := 0.0
	if v > 0 {
		sd = math.Sqrt(v)
	}
	for i := range sites {
		if sd > 0 && !math.IsNaN(sites[i].UnstdXPEHH) && !math.IsInf(sites[i].UnstdXPEHH, 0) {
			sites[i].StdXPEHH = (sites[i].UnstdXPEHH - m) / sd
		} else {
			sites[i].StdXPEHH = math.NaN()
		}
	}
}

// computeLDDecay bins mean r^2 over same-chromosome site pairs by their cM
// separation. Pairs are enumerated within a bounded window of neighbouring sites
// per core so the pass is linear-ish rather than O(sites^2).
func computeLDDecay(cw map[int]*chromHaps, chroms []int) []LDBin {
	edges := []float64{0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64}
	sums := make([]float64, len(edges)-1)
	cnts := make([]int, len(edges)-1)
	const window = 40 // neighbours ahead per site
	for _, ch := range chroms {
		c := cw[ch]
		nhap := c.nhap
		if nhap < 2 {
			continue
		}
		for i := 0; i < len(c.sites); i++ {
			for j := i + 1; j <= i+window && j < len(c.sites); j++ {
				r2 := pairR2(c, c.sites[i].mutID, c.sites[j].mutID, nhap)
				if math.IsNaN(r2) {
					continue
				}
				dcm := math.Abs(c.sites[j].cumCM - c.sites[i].cumCM)
				b := binIndex(edges, dcm)
				if b < 0 {
					continue
				}
				sums[b] += r2
				cnts[b]++
			}
		}
	}
	out := make([]LDBin, 0, len(edges)-1)
	for b := 0; b < len(edges)-1; b++ {
		mean := 0.0
		if cnts[b] > 0 {
			mean = sums[b] / float64(cnts[b])
		}
		out = append(out, LDBin{LoCM: edges[b], HiCM: edges[b+1], MeanR2: mean, NumPairs: cnts[b]})
	}
	return out
}

func binIndex(edges []float64, x float64) int {
	for b := 0; b < len(edges)-1; b++ {
		if x >= edges[b] && x < edges[b+1] {
			return b
		}
	}
	return -1
}

// pairR2 computes the standard biallelic r^2 between two pool sites over all
// strands: pA/pB derived freqs, pAB co-carry freq, D=pAB-pA*pB,
// r^2 = D^2 / (pA(1-pA)pB(1-pB)). NaN when either site is monomorphic.
func pairR2(c *chromHaps, mA, mB, nhap int) float64 {
	var nA, nB, nAB int
	for h := 0; h < nhap; h++ {
		a := c.member[h][mA]
		b := c.member[h][mB]
		if a {
			nA++
		}
		if b {
			nB++
		}
		if a && b {
			nAB++
		}
	}
	pA := float64(nA) / float64(nhap)
	pB := float64(nB) / float64(nhap)
	if pA <= 0 || pA >= 1 || pB <= 0 || pB >= 1 {
		return math.NaN()
	}
	pAB := float64(nAB) / float64(nhap)
	D := pAB - pA*pB
	den := pA * (1 - pA) * pB * (1 - pB)
	if den <= 0 {
		return math.NaN()
	}
	return D * D / den
}
