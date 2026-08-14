package analysis

// TMR-K-A — the ground-truth lineage-count walker. TMR4A.md work item W5.
//
// WHAT THIS COMPUTES, AND WHY IT IS NOT AN ESTIMATE. TMR4A as published is the time at
// which the marginal genealogy at a locus has exactly four lineages, INFERRED by
// ARGweaver's MCMC from sequence data under a coalescent prior with an assumed effective
// population size. This module computes the same quantity from the RECORDED genealogy:
// DRIFT knows who every parent was (W1) and which parental strand every locus came from
// (W2), so the marginal tree at a focal locus is not reconstructed, it is read. No MCMC,
// no prior, no inference error. That is the whole point — it is the truth that an
// ARGweaver run on the same data (W6/W7) is then scored against.
//
// K IS A PARAMETER, NOT A CONSTANT. TMR4A.md §0: "four" is the ceiling only if each
// founder's germline is a faithful copy of a single diploid genome. If the created
// couple's reproductive cells carried variation among themselves, the founding pair can
// transmit MORE than four alleles per locus and 4 is simply the wrong bound to test. So
// the walker emits the full lineage-count-versus-year trajectory and TMR-K-A for any K
// is a lookup against it; tmrka_k selects which one goes in the headline column.
//
// THE THREE OUTCOMES, WHICH MUST NEVER BE COLLAPSED (TMR4A.md §0, §6 risk 5):
//
//   coalesced            the lineage count fell to K by descent, at a known year;
//   censored_at_founding the walk reached the founding generation with MORE than K
//                        lineages still distinct — there is no TMR-K-A, because those
//                        lineages never coalesce. Reporting this as a date would
//                        reproduce exactly the error the module exists to expose;
//   no_coalescence       the walk reached founding without a single coalescence: K was
//                        never approached at all.
//
// A created allele has no coalescence time with another created allele. There is no
// mutational path between them and no year at which they meet. The walk therefore
// TERMINATES at founders and says so, rather than forcing a root.
//
// THE HEADLINE FALSIFIABLE NUMBER is not a date. It is FounderStrandsSurviving: the
// number of distinct founder haplotypes still present in the sample at the locus. If
// that count is at or below the founding couple's transmissible allele count, a
// two-person origin is consistent with the data at that locus whatever an inferred
// TMR4A says. Note what this is exactly: distinct founder (individual, strand) SOURCES,
// which is the identity-by-DESCENT count. Once created alleles can coincide by chance
// (W3 labels, W4 created alleles) the identity-by-STATE count can be smaller, and
// CreatedAllelesSurviving will carry it; until W3 lands that field is reported as -1
// rather than silently conflated with this one.
//
// TIME CONVENTION. Two lineages coalesce when they arrive at the same (individual,
// strand). The event is dated to that ANCESTOR'S BIRTH YEAR. The alternative convention
// — the birth years of the two children whose gametes copied that strand — differs by
// less than one generation and would make the trajectory non-monotone in the walk order;
// the ancestor's birth year is the node's own time and is what a tree-sequence would
// record.
//
// COST AND DETERMINISM. Read-only over Pedigree + ARGDB. Consumes no RNG when
// tmrka_sample_size is 0 (a full-population walk), so a track_tmrka run is byte-identical
// to a track_tmrka-off run. Per locus the walk is a priority queue over 2n lineages,
// stepped one generation at a time in strictly decreasing birth-year order, so it is
// O(steps · log n) and independent of genome size. Loci are walked independently and in
// sorted order, which is what makes the locus-sweep strategy in TMR4A.md §5 safe: the
// result at a locus does not depend on which other loci were in the set.

import (
	"container/heap"
	"drift/pkg/core"
	"sort"
)

// TMRKA outcome classes. Kept as strings because they travel all the way into the CSV,
// the plots, and the prose, and a bare integer code is exactly how a censored walk gets
// silently reported as a date.
const (
	OutcomeCoalesced = "coalesced"
	OutcomeCensored  = "censored_at_founding"
	OutcomeNoCoal    = "no_coalescence"
)

// TMRKALocus is one focal locus's ground-truth genealogy summary.
type TMRKALocus struct {
	Position   int    // genome-bit position of the focal locus
	K          int    // the K in TMR-K-A for the headline columns
	Outcome    string // one of the three classes above
	SampleYear int    // simulated year the sample was taken (the walk's t=0)

	InitialLineages int // 2 × sampled individuals
	FinalLineages   int // lineages still distinct when the walk hit founding

	TMRKAYear    int  // year the count first reached K; only meaningful if HasTMRKA
	HasTMRKA     bool // false when the walk was censored above K
	TMRCAYear    int  // year the count first reached 1; only meaningful if HasTMRCA
	HasTMRCA     bool // false unless the whole sample coalesced to one lineage
	Coalescences int  // number of coalescence events observed

	// Distinct founder (individual, strand) sources the surviving lineages trace to —
	// the identity-by-descent count of founder haplotypes in the sample at this locus.
	FounderStrandsSurviving int
	// Distinct CREATED ALLELES surviving, i.e. the identity-by-state-at-founding count,
	// from the W3 labels. -1 when labelling is off; never silently equated with the
	// by-descent count above, because the whole argument turns on telling them apart.
	// It is ≤ FounderStrandsSurviving: several founder haplotypes can carry the same
	// created allele, and then their descendants' differences are purely mutational.
	CreatedAllelesSurviving int
	// Surviving lineages that froze somewhere other than a labelled founder (an
	// individual whose ancestry was never recorded). Each is counted as its own unknown
	// created allele above, since assuming shared identity is the error to avoid. 0 in a
	// well-formed run; non-zero means the labelled generation and the generation the
	// walk terminates at have come apart. -1 when labelling is off.
	UnlabelledSources int

	// Pairwise decomposition over the C(2n,2) sampled lineage pairs at this locus.
	//
	// CoalescedPairs — pairs that meet at a common ancestor, so any difference between
	// them accumulated by mutation since. Always computed.
	//
	// CreatedPairs — pairs whose lineages end on founder strands carrying DIFFERENT
	// created alleles. Their divergence at this locus was created, not accumulated: no
	// mutational path connects them and no year exists at which they coalesce. This is
	// the quantity §0 says a coalescent model converts into apparent age, so its share
	// of all pairs is the study's key per-locus number. -1 when labelling is off.
	//
	// The remainder — pairs that neither coalesce nor differ by created allele — are
	// distinct founder haplotypes that happen to carry the SAME created allele. Their
	// differences are mutational too, which is why CreatedPairs and not
	// (total − CoalescedPairs) is the created share.
	CoalescedPairs int
	CreatedPairs   int

	// Sequence-difference counts at the locus (§4). These are identity-by-state counts
	// over actual sequence and need W4's `founder_allele_divergence` to exist; -1 until
	// then. Deliberately NOT filled with the pair counts above, which measure something
	// different.
	MutationalDiffs int
	CreatedDiffs    int

	// Trajectory of the lineage count through time, newest first: one entry per
	// coalescence, holding the year and the count AFTER that coalescence. TMR-K-A for
	// any K is a lookup against this, which is what makes K a parameter rather than a
	// design centre.
	Trajectory []TMRKAPoint
}

// TMRKAPoint is one step of the lineage-count-versus-year trajectory.
type TMRKAPoint struct {
	Year     int // year of the coalescence (the ancestor's birth year)
	Lineages int // lineages remaining after it
}

// TMRKAResult is a whole run's walk over the focal locus set.
type TMRKAResult struct {
	K              int
	SampleYear     int
	NumSampled     int // individuals in the sample (2× this many starting lineages)
	Loci           []TMRKALocus
	MissingARG     int // loci skipped because no provenance was recorded
	GenerationTime float64
}

// LineagesAtYear returns the number of lineages remaining at (or immediately after) the
// given year — the trajectory read as a step function. Used by the output writer and by
// anyone asking for a K the run did not headline.
func (l *TMRKALocus) LineagesAtYear(year int) int {
	n := l.InitialLineages
	for _, p := range l.Trajectory {
		if p.Year > year {
			break
		}
		n = p.Lineages
	}
	return n
}

// TMRKAForK returns the year the count first fell to K or below, for any K, from the
// recorded trajectory — the "any K is a lookup" property. ok is false when the walk
// never got that low (a censored locus).
func (l *TMRKALocus) TMRKAForK(k int) (int, bool) {
	if l.InitialLineages <= k {
		return l.SampleYear, true
	}
	for _, p := range l.Trajectory {
		if p.Lineages <= k {
			return p.Year, true
		}
	}
	return 0, false
}

// lineage is one ancestral lineage's current position: a specific strand of a specific
// individual. Two lineages at the same position have coalesced.
// A CREATED terminus (TMR4A.md W4) is also a lineageKey: `created` set, `ind` the founder
// whose germline the allele came from, and `strand` carrying the created-allele index
// instead of 0/1. That is not a trick — it is the same object. Two lineages that reach the
// same created allele in the same founder have met (one created sequence, one germline,
// no mutational difference between the copies), so they must coalesce there, and reusing
// the key makes the existing collision test do that for free. Two lineages at DIFFERENT
// created alleles are different keys and never meet, which is the censoring §0 requires.
type lineageKey struct {
	ind     int
	strand  int
	created bool
}

// walkItem is a heap entry. The heap pops the YOUNGEST lineage first (largest birth
// year), because a lineage can only ever move to a strictly older node, so stepping the
// youngest first guarantees no meeting is missed.
type walkItem struct {
	key  lineageKey
	year int
}

type walkHeap []walkItem

func (h walkHeap) Len() int { return len(h) }
func (h walkHeap) Less(i, j int) bool {
	// Birth year descending; ties broken by id then strand so the walk order — and
	// therefore the trajectory — is fully determined by the data, not by map order.
	//
	// The id tie-break is load-bearing, not cosmetic. The walk's correctness rests on
	// never arriving at a node that has already been stepped past, which holds because
	// every arrival is strictly older than the node just popped. If a model ever allows
	// a parent and child to share a birth year, descending id preserves that: DRIFT
	// allocates ids from a monotonically increasing counter, so a child's id always
	// exceeds both parents' and the child is still popped first.
	if h[i].year != h[j].year {
		return h[i].year > h[j].year
	}
	if h[i].key.ind != h[j].key.ind {
		return h[i].key.ind > h[j].key.ind
	}
	return h[i].key.strand > h[j].key.strand
}
func (h walkHeap) Swap(i, j int)       { h[i], h[j] = h[j], h[i] }
func (h *walkHeap) Push(x interface{}) { *h = append(*h, x.(walkItem)) }
func (h *walkHeap) Pop() interface{} {
	old := *h
	n := len(old)
	it := old[n-1]
	*h = old[:n-1]
	return it
}

// ComputeTMRKA walks the true genealogy at every focal locus for the given sample.
// Returns nil when the pedigree or the focal locus set is missing — the walk has no
// substrate and inventing one would be worse than reporting nothing.
func ComputeTMRKA(model *core.Model, pop *core.Pop, ids []int) *TMRKAResult {
	loci := core.EnsureARGLoci(model)
	if pop == nil || pop.Pedigree == nil || len(loci) == 0 || len(ids) == 0 {
		return nil
	}
	k := int(model.Parameters["tmrka_k"])
	if k < 1 {
		k = 4
	}
	genTime := model.Parameters["generation_time"]
	if genTime <= 0 {
		genTime = 25
	}

	sampleYear := model.FreeParameters["year"]
	sorted := append([]int(nil), ids...)
	sort.Ints(sorted)

	res := &TMRKAResult{
		K:              k,
		SampleYear:     sampleYear,
		NumSampled:     len(sorted),
		GenerationTime: genTime,
	}
	for i, pos := range loci {
		loc := walkLocus(pop, sorted, i, len(loci), pos, k, sampleYear)
		if loc == nil {
			res.MissingARG++
			continue
		}
		res.Loci = append(res.Loci, *loc)
	}
	return res
}

// walkLocus is the walk itself, for one focal locus (index i into the locus set).
func walkLocus(pop *core.Pop, ids []int, locusIdx, nLoci, pos, k, sampleYear int) *TMRKALocus {
	// Seed one lineage per sampled strand. Distinct sampled individuals are distinct
	// nodes, so the initial set has exactly 2n members.
	active := make(map[lineageKey]bool, 2*len(ids))
	// weight[key] = how many of the ORIGINAL sampled lineages currently sit at this node.
	// It rides along with each lineage as it steps back and accumulates on coalescence,
	// so at the end every surviving source knows how much of the sample descends from it.
	// That is what turns the walk into a pairwise decomposition: two sampled lineages
	// share an ancestor exactly when they end up under the same source, and the pair
	// counts follow from the weights without re-walking anything.
	weight := make(map[lineageKey]int, 2*len(ids))
	h := &walkHeap{}
	for _, id := range ids {
		n, ok := pop.Pedigree[id]
		year := sampleYear
		if ok {
			year = n.BirthYear
		}
		for s := 0; s < 2; s++ {
			key := lineageKey{ind: id, strand: s}
			active[key] = true
			weight[key] = 1
			*h = append(*h, walkItem{key: key, year: year})
		}
	}
	heap.Init(h)

	loc := &TMRKALocus{
		Position:                pos,
		K:                       k,
		SampleYear:              sampleYear,
		InitialLineages:         len(active),
		CreatedAllelesSurviving: -1, // filled by summariseSources when W3 labels are on
		UnlabelledSources:       -1,
		CreatedPairs:            -1,
		MutationalDiffs:         -1, // sequence-difference counts; need W4 divergence
		CreatedDiffs:            -1,
	}
	// Coalescence YEARS, collected in walk order and sorted afterwards. They are not
	// collected in time order: the walk pops by the DESCENDANT's birth year but dates
	// the event to the ANCESTOR's birth year, and with DRIFT's long founder lifespans a
	// late birth to an ancient parent can produce an older event than one discovered
	// later. Sorting before building the trajectory is what makes "lineage count versus
	// time" mean what it says.
	var events []int

	// Lineages that can go no further: they sit on a founder strand. Held aside so they
	// remain available for a later lineage to coalesce WITH — arriving at a founder
	// strand another lineage already occupies is a genuine coalescence — while never
	// being stepped again.
	frozen := make(map[lineageKey]bool)

	// Safety bound: the pedigree is a DAG (a parent is always born before its child) so
	// the walk strictly decreases in time and must terminate. The cap only protects
	// against a corrupted pedigree turning that into a hang.
	maxSteps := 64 * (loc.InitialLineages + 1) * (len(pop.Pedigree) + 1)

	for steps := 0; h.Len() > 0 && steps < maxSteps; steps++ {
		it := heap.Pop(h).(walkItem)
		if !active[it.key] {
			continue // already merged away by an earlier coalescence
		}

		next, ok := stepBack(pop, it.key, locusIdx, nLoci)
		if !ok {
			// Founding generation reached on this lineage: terminate it HERE rather
			// than forcing a root. This is the censoring that TMR4A.md §0 insists must
			// be reportable.
			frozen[it.key] = true
			continue
		}
		pn, pok := pop.Pedigree[next.ind]
		if !pok {
			frozen[it.key] = true
			continue
		}

		delete(active, it.key)
		w := weight[it.key]
		delete(weight, it.key)
		if active[next] || frozen[next] {
			// Coalescence: two lineages have met in the same strand of the same
			// ancestor. Dated to the ancestor's birth year. The arriving lineage's
			// share of the sample merges into the node it met.
			events = append(events, pn.BirthYear)
			weight[next] += w
			continue
		}
		active[next] = true
		weight[next] = w
		heap.Push(h, walkItem{key: next, year: pn.BirthYear})
	}

	// Everything still standing is frozen at a founder strand (or was merged into one).
	// These are the distinct founder haplotypes the sample traces to at this locus —
	// the headline falsifiable number. Recomputing it from the surviving lineage
	// positions rather than from the event count is a deliberate cross-check: the two
	// must agree, and TestTMRKAFounderStrandsMatchLineageCount asserts they do.
	sources := make(map[lineageKey]bool, len(frozen))
	for key := range frozen {
		sources[key] = true
	}
	for key := range active {
		sources[key] = true
	}
	loc.FounderStrandsSurviving = len(sources)
	summariseSources(pop, loc, sources, weight, locusIdx)

	// Sort the events into true time order (most recent first) before turning them into
	// a trajectory, then read every derived quantity off that trajectory. SliceStable
	// keeps same-year events in deterministic walk order.
	sort.SliceStable(events, func(i, j int) bool { return events[i] > events[j] })
	loc.Coalescences = len(events)
	loc.FinalLineages = loc.InitialLineages - len(events)
	loc.Trajectory = make([]TMRKAPoint, 0, len(events))
	for i, year := range events {
		loc.Trajectory = append(loc.Trajectory, TMRKAPoint{Year: year, Lineages: loc.InitialLineages - i - 1})
	}
	loc.TMRKAYear, loc.HasTMRKA = loc.TMRKAForK(k)
	loc.TMRCAYear, loc.HasTMRCA = loc.TMRKAForK(1)

	switch {
	case loc.HasTMRKA:
		loc.Outcome = OutcomeCoalesced
	case loc.Coalescences == 0:
		loc.Outcome = OutcomeNoCoal
	default:
		loc.Outcome = OutcomeCensored
	}
	return loc
}

// summariseSources turns the surviving lineages and their sample weights into the
// created-allele counts and the pairwise decomposition (TMR4A.md W3).
//
// The pair arithmetic, stated once so the code below is checkable by eye. Write T for the
// sample's lineage count and w_i for the share of it descending from source i:
//
//	total pairs      C(T,2)
//	coalesced pairs  Σ C(w_i,2)                       — pairs meeting at a common ancestor
//	created pairs    Σ_{ℓ<ℓ'} W_ℓ·W_ℓ' = (T² − ΣW_ℓ²)/2  — pairs in different created alleles
//
// where W_ℓ sums w_i over the sources carrying created allele ℓ. Grouping by label rather
// than by source is the entire point: two sources with the SAME label contribute nothing
// to the created total, because their descendants' differences are mutational even though
// the two lineages never coalesce. An unlabelled source is given a private label so it
// can never be assumed to share an allele with anything.
func summariseSources(pop *core.Pop, loc *TMRKALocus, sources map[lineageKey]bool, weight map[lineageKey]int, locusIdx int) {
	// Deterministic order: map iteration is randomized and the unlabelled sources below
	// are assigned private labels by position.
	keys := make([]lineageKey, 0, len(sources))
	for key := range sources {
		keys = append(keys, key)
	}
	sort.Slice(keys, func(i, j int) bool {
		if keys[i].ind != keys[j].ind {
			return keys[i].ind < keys[j].ind
		}
		return keys[i].strand < keys[j].strand
	})

	for _, key := range keys {
		w := weight[key]
		loc.CoalescedPairs += w * (w - 1) / 2
	}

	if !pop.FounderLabelsOn() {
		loc.CreatedAllelesSurviving = -1
		loc.CreatedPairs = -1
		loc.UnlabelledSources = -1
		return
	}

	// Labelling is on, so the count of unlabelled sources is a real measurement of zero
	// or more — clear the "absent" sentinel before tallying rather than incrementing it.
	loc.UnlabelledSources = 0

	// Private labels for unlabelled sources start above every real label so they cannot
	// collide with one.
	byLabel := map[int]int{}
	nextPrivate := 0
	for _, key := range keys {
		if label, ok := sourceLabel(pop, key, locusIdx); ok {
			byLabel[label] += weight[key]
			if label >= nextPrivate {
				nextPrivate = label + 1
			}
		}
	}
	for _, key := range keys {
		if _, ok := sourceLabel(pop, key, locusIdx); !ok {
			loc.UnlabelledSources++
			byLabel[nextPrivate] = weight[key]
			nextPrivate++
		}
	}

	loc.CreatedAllelesSurviving = len(byLabel)

	total := 0
	sumSquares := 0
	for _, w := range byLabel {
		total += w
		sumSquares += w * w
	}
	loc.CreatedPairs = (total*total - sumSquares) / 2
}

// sourceLabel returns the created-allele identity a terminated lineage carries. For a
// created terminus (W4) the key already IS the allele — nothing to look up. For a founder
// strand it is the W3 label, which under the created model comes from that founder's pool
// so the two items answer with the same numbering and a run mixing both stays coherent.
func sourceLabel(pop *core.Pop, key lineageKey, locusIdx int) (int, bool) {
	if key.created {
		return key.strand, true
	}
	return pop.FounderLabel(key.ind, key.strand, locusIdx)
}

// stepBack maps one lineage position to the next one back at a focal locus — the single
// elementary move of the walk, and the exact place W1, W2 and W4 meet. Child strand 0 came
// from the paternal gamete and strand 1 from the maternal gamete (birth.go's convention),
// so the strand index doubles as the gamete index into the provenance record. ok=false
// means the walk stops here: a created allele, a founder, or an individual whose parent or
// provenance was never recorded.
func stepBack(pop *core.Pop, key lineageKey, locusIdx, nLoci int) (lineageKey, bool) {
	// A created allele has no ancestor by construction. This is the terminal case the
	// whole module is about: there is no earlier node and no year at which it meets
	// anything else, so the walk must stop rather than invent one.
	if key.created {
		return lineageKey{}, false
	}
	n, exists := pop.Pedigree[key.ind]
	if !exists {
		return lineageKey{}, false
	}
	parent := n.Dad
	if key.strand == 1 {
		parent = n.Mom
	}
	if parent == core.PedFounder {
		return lineageKey{}, false
	}
	// A gamete drawn from a founder's created pool (W4) resolves to a created ALLELE, not
	// to one of that founder's two somatic strands — checked first, because under the
	// created model the ordinary strand bit is also recorded but describes only which of
	// the germ cell's two alleles was taken, not which allele that is.
	if allele, ok := pop.CreatedProvenance(key.ind, key.strand, locusIdx, nLoci); ok {
		return lineageKey{ind: parent, strand: allele, created: true}, true
	}
	parentStrand, ok := pop.ARGStrand(key.ind, key.strand, locusIdx, nLoci)
	if !ok {
		return lineageKey{}, false
	}
	return lineageKey{ind: parent, strand: parentStrand}, true
}
