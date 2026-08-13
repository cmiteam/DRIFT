package analysis

import (
	"drift/pkg/core"
	"testing"
)

// TMR-K-A walker validation (TMR4A.md W5 / W8). The decisive tests here are against
// HAND-BUILT pedigrees whose coalescence times are known by construction — the same
// hand-built-panel discipline used for §6b/§6g — because the walker's whole claim is
// that it reports the truth rather than an estimate, and a test that only checks
// self-consistency could not tell a correct walk from a plausible one.

// tmrkaBuilder assembles a pedigree plus strand provenance for a fixed focal locus set.
type tmrkaBuilder struct {
	pop  *core.Pop
	loci []int
}

func newTMRKA(loci []int) *tmrkaBuilder {
	return &tmrkaBuilder{pop: &core.Pop{IndData: map[int][]int{}}, loci: loci}
}

func (b *tmrkaBuilder) founder(id, birthYear int) {
	b.pop.PedEnsure(id, birthYear, id%2)
}

// child records a birth and its provenance directly in strand terms: dadStrands[i] is
// the paternal strand inherited at focal locus i, momStrands[i] the maternal one. The
// masks are synthesised from those so the capture path under test is the real one.
func (b *tmrkaBuilder) child(id, dad, mom, birthYear int, dadStrands, momStrands []int) {
	b.pop.PedRecordBirth(id, dad, mom, birthYear, id%2)
	maxPos := 0
	for _, p := range b.loci {
		if p > maxPos {
			maxPos = p
		}
	}
	words := maxPos/64 + 1
	mkMask := func(strands []int) []uint64 {
		m := make([]uint64, words)
		for i, pos := range b.loci {
			// mask bit SET ⇒ meiosis takes parental copy 0 (see pkg/core/arg.go)
			if strands[i] == 0 {
				m[pos/64] |= 1 << (uint(pos) % 64)
			}
		}
		return m
	}
	b.pop.ARGCapture(id, b.loci, mkMask(dadStrands), mkMask(momStrands))
}

func (b *tmrkaBuilder) alive(ids ...int) {
	for _, id := range ids {
		b.pop.IndData[id] = []int{}
	}
}

func (b *tmrkaBuilder) model(k, year int) *core.Model {
	return &core.Model{
		Parameters:     map[string]float64{"tmrka_k": float64(k), "generation_time": 25},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"year": year, "genome_bits": 1024},
		ARGLoci:        b.loci,
	}
}

// Two full sibs. At a locus where both took the SAME grandparental strand from each
// parent, all four sampled lineages coalesce to two in the parents — never below two,
// because the two parents are unrelated founders. So with K=2 the walk coalesces at the
// parents' birth year, and with K=1 it is censored: there is no TMRCA.
func TestTMRKAKnownSibCoalescence(t *testing.T) {
	b := newTMRKA([]int{10})
	b.founder(1, -30) // dad
	b.founder(2, -30) // mom
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 5, []int{0}, []int{0})
	b.alive(3, 4)

	res := ComputeTMRKA(b.model(2, 100), b.pop, []int{3, 4})
	if res == nil || len(res.Loci) != 1 {
		t.Fatalf("expected one locus, got %+v", res)
	}
	l := res.Loci[0]
	if l.InitialLineages != 4 {
		t.Fatalf("initial lineages: got %d, want 4", l.InitialLineages)
	}
	if l.Coalescences != 2 {
		t.Fatalf("coalescences: got %d, want 2 (one in dad, one in mom)", l.Coalescences)
	}
	if l.Outcome != OutcomeCoalesced || !l.HasTMRKA {
		t.Fatalf("outcome: got %q (has=%v), want coalesced", l.Outcome, l.HasTMRKA)
	}
	if l.TMRKAYear != -30 {
		t.Errorf("TMR-2-A year: got %d, want -30 (the parents' birth year)", l.TMRKAYear)
	}
	if l.FinalLineages != 2 || l.FounderStrandsSurviving != 2 {
		t.Errorf("founder haplotypes surviving: got %d/%d, want 2/2 — two unrelated founders "+
			"contributed one strand each and those never coalesce",
			l.FinalLineages, l.FounderStrandsSurviving)
	}
	if l.HasTMRCA {
		t.Errorf("reported a TMRCA of %d; two created founder strands have no coalescence time", l.TMRCAYear)
	}
}

// The same two sibs at a locus where they took DIFFERENT grandparental strands from
// each parent: nothing coalesces at all. This is the `no_coalescence` class — K was
// never approached — and it must not be reported as a date.
func TestTMRKANoCoalescence(t *testing.T) {
	b := newTMRKA([]int{10})
	b.founder(1, -30)
	b.founder(2, -30)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 5, []int{1}, []int{1})
	b.alive(3, 4)

	l := ComputeTMRKA(b.model(2, 100), b.pop, []int{3, 4}).Loci[0]
	if l.Outcome != OutcomeNoCoal {
		t.Errorf("outcome: got %q, want %q", l.Outcome, OutcomeNoCoal)
	}
	if l.HasTMRKA {
		t.Errorf("censored walk reported TMR-K-A year %d", l.TMRKAYear)
	}
	if l.FinalLineages != 4 {
		t.Errorf("lineages at founding: got %d, want 4 (all four founder strands survive)", l.FinalLineages)
	}
}

// Censored ABOVE K: some coalescence happens, but the walk hits creation with more than
// K lineages left. Distinguishing this from `coalesced` is TMR4A.md §6 risk 5.
func TestTMRKACensoredAtFounding(t *testing.T) {
	b := newTMRKA([]int{10})
	b.founder(1, -30)
	b.founder(2, -30)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 5, []int{0}, []int{1}) // shares dad's strand 0, differs on mom
	b.alive(3, 4)

	l := ComputeTMRKA(b.model(1, 100), b.pop, []int{3, 4}).Loci[0]
	if l.Coalescences != 1 {
		t.Fatalf("coalescences: got %d, want 1", l.Coalescences)
	}
	if l.Outcome != OutcomeCensored {
		t.Errorf("outcome: got %q, want %q", l.Outcome, OutcomeCensored)
	}
	if l.HasTMRKA {
		t.Errorf("censored walk reported a TMR-1-A year of %d", l.TMRKAYear)
	}
	if l.FinalLineages != 3 {
		t.Errorf("lineages at founding: got %d, want 3", l.FinalLineages)
	}
}

// Three generations with a known, staged answer: the count falls 6 → 4 → 3 → 2 at
// birth years the test dictates, so every trajectory entry and every K lookup can be
// checked against arithmetic rather than against the walker's own output.
func TestTMRKAKnownTrajectoryThreeGenerations(t *testing.T) {
	b := newTMRKA([]int{0})
	// Founders (year -50): 1 (m), 2 (f), 5 (f)
	b.founder(1, -50)
	b.founder(2, -50)
	b.founder(5, -50)
	// Generation 1 (year 0): sibs 3 and 4 from 1x2, both taking dad's strand 0 and
	// mom's strand 0 — so the two paternal lineages meet in 1, the two maternal in 2.
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 0, []int{0}, []int{0})
	// Generation 2 (year 30): 6 and 7, children of 3x5 and 4x5. Both take their
	// father's strand 0 (which traces to founder 1) and their mother's strand 0.
	b.founder(8, -50) // unused male founder keeps ids from being suggestive
	b.child(6, 3, 5, 30, []int{0}, []int{0})
	b.child(7, 4, 5, 30, []int{0}, []int{0})
	b.alive(6, 7)

	res := ComputeTMRKA(b.model(4, 60), b.pop, []int{6, 7})
	l := res.Loci[0]
	if l.InitialLineages != 4 {
		t.Fatalf("initial lineages: got %d, want 4", l.InitialLineages)
	}
	// Maternal lineages of 6 and 7 both sit on strand 0 of founder 5 ⇒ coalesce at -50.
	// Paternal lineages run 6→3→1 and 7→4→1, both onto strand 0 of founder 1 ⇒ -50.
	if l.Coalescences != 2 {
		t.Fatalf("coalescences: got %d, want 2; trajectory=%+v", l.Coalescences, l.Trajectory)
	}
	for _, p := range l.Trajectory {
		if p.Year != -50 {
			t.Errorf("coalescence dated %d, want -50 (the founders' birth year)", p.Year)
		}
	}
	if y, ok := l.TMRKAForK(3); !ok || y != -50 {
		t.Errorf("TMR-3-A: got (%d,%v), want (-50,true)", y, ok)
	}
	if y, ok := l.TMRKAForK(2); !ok || y != -50 {
		t.Errorf("TMR-2-A: got (%d,%v), want (-50,true)", y, ok)
	}
	if _, ok := l.TMRKAForK(1); ok {
		t.Error("TMR-1-A must not exist: founders 1 and 5 are unrelated created individuals")
	}
	if l.FinalLineages != 2 {
		t.Errorf("founder haplotypes surviving: got %d, want 2", l.FinalLineages)
	}
	// K larger than the sample is degenerate but must still be answered, at t=0.
	if y, ok := l.TMRKAForK(10); !ok || y != 60 {
		t.Errorf("TMR-10-A on a 4-lineage sample: got (%d,%v), want (60,true)", y, ok)
	}
}

// The trajectory must be ordered by time, not by the order the walk happened to find
// events. With DRIFT's long founder lifespans a late birth to an ancient parent yields
// an older event than one discovered later, so an unsorted trajectory would attach the
// wrong lineage counts to the wrong years.
func TestTMRKATrajectoryIsTimeOrdered(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, -400) // ancient, long-lived
	b.founder(2, -400)
	b.founder(5, -10) // recent
	b.founder(6, -10)
	// A pair of sibs from the ancient couple, born LATE (year 50).
	b.child(10, 1, 2, 50, []int{0}, []int{0})
	b.child(11, 1, 2, 50, []int{0}, []int{0})
	// A pair of sibs from the recent couple, born EARLY (year 20).
	b.child(12, 5, 6, 20, []int{0}, []int{0})
	b.child(13, 5, 6, 20, []int{0}, []int{0})
	b.alive(10, 11, 12, 13)

	l := ComputeTMRKA(b.model(4, 100), b.pop, []int{10, 11, 12, 13}).Loci[0]
	if l.InitialLineages != 8 || l.Coalescences != 4 {
		t.Fatalf("got %d lineages, %d coalescences; want 8 and 4", l.InitialLineages, l.Coalescences)
	}
	prev := l.SampleYear
	for i, p := range l.Trajectory {
		if p.Year > prev {
			t.Errorf("trajectory entry %d is at year %d, newer than the previous %d", i, p.Year, prev)
		}
		prev = p.Year
		if want := l.InitialLineages - i - 1; p.Lineages != want {
			t.Errorf("trajectory entry %d: lineages %d, want %d", i, p.Lineages, want)
		}
	}
	// The two recent-founder coalescences (-10) must come before the ancient ones (-400).
	if l.Trajectory[0].Year != -10 || l.Trajectory[1].Year != -10 {
		t.Errorf("first two events at %d/%d, want -10/-10",
			l.Trajectory[0].Year, l.Trajectory[1].Year)
	}
	if l.Trajectory[2].Year != -400 || l.Trajectory[3].Year != -400 {
		t.Errorf("last two events at %d/%d, want -400/-400",
			l.Trajectory[2].Year, l.Trajectory[3].Year)
	}
	if y, ok := l.TMRKAForK(6); !ok || y != -10 {
		t.Errorf("TMR-6-A: got (%d,%v), want (-10,true)", y, ok)
	}
	if y, ok := l.TMRKAForK(4); !ok || y != -400 {
		t.Errorf("TMR-4-A: got (%d,%v), want (-400,true)", y, ok)
	}
}

// A sampled individual who is also an ancestor of another sampled individual must be
// handled: the descendant's lineage arrives at a node that is itself a starting lineage.
func TestTMRKASampleContainsItsOwnAncestor(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, -60)
	b.founder(2, -60)
	b.child(3, 1, 2, 0, []int{0}, []int{0}) // in the sample AND a parent of 5
	b.founder(4, -60)
	b.child(5, 3, 4, 40, []int{0}, []int{0}) // paternal strand 0 of 3
	b.alive(3, 5)

	l := ComputeTMRKA(b.model(2, 80), b.pop, []int{3, 5}).Loci[0]
	// 5's paternal lineage lands on (3, strand 0), which is a sampled lineage: they
	// coalesce in 3 itself, at 3's birth year.
	if l.Coalescences < 1 {
		t.Fatalf("no coalescence found; trajectory=%+v", l.Trajectory)
	}
	if l.Trajectory[0].Year != 0 {
		t.Errorf("first coalescence at %d, want 0 (individual 3's own birth year)", l.Trajectory[0].Year)
	}
}

// Loci are walked independently: the result at a locus must not depend on which other
// loci were in the set. This is what makes the locus-sweep escape hatch in TMR4A.md §5
// sound, and W8 requires it to be actually tested rather than assumed.
func TestTMRKALociAreIndependent(t *testing.T) {
	build := func(loci []int) *tmrkaBuilder {
		b := newTMRKA(loci)
		b.founder(1, -40)
		b.founder(2, -40)
		b.founder(5, -40)
		n := len(loci)
		alt := func(offset int) []int {
			s := make([]int, n)
			for i := range s {
				s[i] = (i + offset) % 2
			}
			return s
		}
		b.child(3, 1, 2, 0, alt(0), alt(1))
		b.child(4, 1, 2, 0, alt(0), alt(0))
		b.child(6, 3, 5, 30, alt(1), alt(0))
		b.child(7, 4, 5, 30, alt(1), alt(1))
		b.alive(6, 7)
		return b
	}

	full := []int{0, 1, 2, 3, 4, 5, 6, 7}
	resFull := ComputeTMRKA(build(full).model(4, 60), build(full).pop, []int{6, 7})
	// Rebuild with only a subset of the loci — a different locus SET, same genealogy.
	subset := []int{2, 5}
	resSub := ComputeTMRKA(build(subset).model(4, 60), build(subset).pop, []int{6, 7})

	byPos := map[int]TMRKALocus{}
	for _, l := range resFull.Loci {
		byPos[l.Position] = l
	}
	for _, l := range resSub.Loci {
		ref, ok := byPos[l.Position]
		if !ok {
			t.Fatalf("locus %d missing from the full run", l.Position)
		}
		if l.Outcome != ref.Outcome || l.Coalescences != ref.Coalescences ||
			l.FinalLineages != ref.FinalLineages || l.TMRKAYear != ref.TMRKAYear {
			t.Errorf("locus %d differs between locus sets: subset %+v vs full %+v", l.Position, l, ref)
		}
	}
}

// Surviving lineages and distinct founder strands are computed by two different routes
// (an event tally versus the positions the lineages actually froze at) and must agree.
func TestTMRKAFounderStrandsMatchLineageCount(t *testing.T) {
	b := newTMRKA([]int{0, 1, 2})
	b.founder(1, -40)
	b.founder(2, -40)
	b.founder(5, -40)
	b.child(3, 1, 2, 0, []int{0, 1, 0}, []int{1, 1, 0})
	b.child(4, 1, 2, 0, []int{0, 0, 1}, []int{1, 0, 0})
	b.child(6, 3, 5, 30, []int{0, 1, 1}, []int{0, 0, 1})
	b.child(7, 4, 5, 30, []int{1, 1, 0}, []int{0, 1, 0})
	b.alive(6, 7)

	res := ComputeTMRKA(b.model(4, 60), b.pop, []int{6, 7})
	for _, l := range res.Loci {
		if l.FinalLineages != l.FounderStrandsSurviving {
			t.Errorf("locus %d: %d lineages left but %d distinct founder strands — "+
				"the two must agree by construction", l.Position, l.FinalLineages, l.FounderStrandsSurviving)
		}
		if l.InitialLineages-l.Coalescences != l.FinalLineages {
			t.Errorf("locus %d: %d - %d != %d", l.Position, l.InitialLineages, l.Coalescences, l.FinalLineages)
		}
	}
}

// No substrate, no output: the walker refuses to invent a genealogy it was not given.
func TestTMRKARequiresSubstrate(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, 0)
	b.alive(1)
	m := b.model(4, 10)

	if got := ComputeTMRKA(m, &core.Pop{}, []int{1}); got != nil {
		t.Error("walked with no pedigree")
	}
	noLoci := b.model(4, 10)
	noLoci.ARGLoci = []int{}
	if got := ComputeTMRKA(noLoci, b.pop, []int{1}); got != nil {
		t.Error("walked with no focal loci")
	}
	if got := ComputeTMRKA(m, b.pop, nil); got != nil {
		t.Error("walked with an empty sample")
	}
}

// A founder in the living sample has no provenance of its own; its lineages must freeze
// where they are rather than being treated as missing data.
func TestTMRKAFounderInSample(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, -20)
	b.founder(2, -20)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.alive(1, 3)

	l := ComputeTMRKA(b.model(4, 40), b.pop, []int{1, 3}).Loci[0]
	if l.InitialLineages != 4 {
		t.Fatalf("initial lineages: got %d, want 4", l.InitialLineages)
	}
	// 3's paternal lineage lands on founder 1 strand 0, which is a sampled lineage.
	if l.Coalescences != 1 {
		t.Errorf("coalescences: got %d, want 1; trajectory=%+v", l.Coalescences, l.Trajectory)
	}
	if l.Trajectory[0].Year != -20 {
		t.Errorf("coalescence year: got %d, want -20", l.Trajectory[0].Year)
	}
}
