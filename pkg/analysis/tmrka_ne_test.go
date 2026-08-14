package analysis

import (
	"drift/pkg/core"
	"math"
	"testing"
)

// Unit tests for the analytic-anchor helpers (TMR4A.md W8). The end-to-end anchor lives
// in pkg/validation; these pin the arithmetic and the edge cases it rests on, on inputs
// small enough to check by hand.

func closeTo(got, want, tol float64) bool { return math.Abs(got-want) <= tol }

// The realised generation time is the mean parent-to-child birth-year gap over BOTH
// parent edges of every child — not the mean age of parents, which weights a childless
// adult the same as one with ten children.
func TestPedigreeGenerationTimeIsOffspringWeighted(t *testing.T) {
	pop := &core.Pop{}
	pop.PedEnsure(1, 0, 0)  // dad
	pop.PedEnsure(2, 10, 1) // mom, born later
	// One child at year 30: gaps of 30 (dad) and 20 (mom).
	pop.PedRecordBirth(3, 1, 2, 30, 0)
	got, edges := PedigreeGenerationTime(pop)
	if edges != 2 {
		t.Fatalf("edges: got %d, want 2 (one per parent)", edges)
	}
	if !closeTo(got, 25, 1e-9) {
		t.Errorf("generation time: got %.4f, want 25 = mean(30,20)", got)
	}

	// A second child of the SAME couple at year 40 adds gaps of 40 and 30, so the couple
	// now counts twice: mean(30,20,40,30) = 30.
	pop.PedRecordBirth(4, 1, 2, 40, 1)
	got, edges = PedigreeGenerationTime(pop)
	if edges != 4 {
		t.Fatalf("edges: got %d, want 4", edges)
	}
	if !closeTo(got, 30, 1e-9) {
		t.Errorf("generation time: got %.4f, want 30 — a second birth must count again", got)
	}
}

// Founders have no recorded parents, so they contribute nothing; a pedigree of founders
// alone has no generation time to report and must say so rather than return 0 as if it
// had measured one.
func TestPedigreeGenerationTimeNeedsEdges(t *testing.T) {
	pop := &core.Pop{}
	pop.PedEnsure(1, -20, 0)
	pop.PedEnsure(2, -20, 1)
	if got, edges := PedigreeGenerationTime(pop); edges != 0 || got != 0 {
		t.Errorf("founders-only pedigree: got %.2f over %d edges, want 0 over 0", got, edges)
	}
	if got, edges := PedigreeGenerationTime(nil); edges != 0 || got != 0 {
		t.Errorf("nil population: got %.2f over %d edges, want 0 over 0", got, edges)
	}
	// A pedigree that never had retention switched on.
	if got, edges := PedigreeGenerationTime(&core.Pop{}); edges != 0 || got != 0 {
		t.Errorf("no pedigree: got %.2f over %d edges, want 0 over 0", got, edges)
	}
}

// A parent whose node is gone (pruned, or never recorded) leaves the child's slot naming
// an id that resolves to nothing. That edge must be skipped, not counted as a zero gap,
// which would drag the mean down silently.
func TestPedigreeGenerationTimeSkipsUnresolvableParents(t *testing.T) {
	pop := &core.Pop{}
	pop.PedEnsure(1, 0, 0)
	pop.PedRecordBirth(2, 1, 99, 20, 1) // mom 99 has no node
	got, edges := PedigreeGenerationTime(pop)
	if edges != 1 || !closeTo(got, 20, 1e-9) {
		t.Errorf("got %.2f over %d edges, want 20 over 1 — the missing parent must be skipped", got, edges)
	}
}

// E[T(n->k)] = 4Ne(1/k - 1/n), and ImpliedCoalescentNe is its exact inverse. A round
// trip is the cheapest possible guard against a transposed constant.
func TestCoalescentExpectationRoundTrips(t *testing.T) {
	cases := []struct{ k, n int }{{1, 2}, {2, 48}, {4, 108}, {8, 96}}
	for _, c := range cases {
		gens, ok := ExpectedCoalescentGens(37.5, c.k, c.n)
		if !ok {
			t.Fatalf("k=%d n=%d: expectation reported undefined", c.k, c.n)
		}
		ne, ok := ImpliedCoalescentNe(gens, c.k, c.n)
		if !ok {
			t.Fatalf("k=%d n=%d: inverse reported undefined", c.k, c.n)
		}
		if !closeTo(ne, 37.5, 1e-9) {
			t.Errorf("k=%d n=%d: round trip gave Ne=%.6f, want 37.5", c.k, c.n, ne)
		}
	}
}

// The value TMR4A.md §1B leans on: at K=4 with the published sample of 108 haplotypes,
// the expected depth is ~0.96*Ne generations, so an inferred TMR4A tracks its assumed Ne
// before any data is consulted. If this number ever moves, §1B's argument moves with it.
func TestCoalescentExpectationMatchesTheS1BFigure(t *testing.T) {
	gens, ok := ExpectedCoalescentGens(1, 4, 108) // Ne=1 => the answer is the multiplier
	if !ok {
		t.Fatal("expectation reported undefined for k=4, n=108")
	}
	if !closeTo(gens, 0.9630, 1e-4) {
		t.Errorf("E[T(108->4)] = %.4f * Ne, want 0.9630 * Ne (TMR4A.md 1B says ~0.96)", gens)
	}
}

// Undefined cases must say so instead of returning a plausible-looking number: a sample
// already at or below K has no waiting time, and neither does a sample of one.
func TestCoalescentExpectationUndefinedCases(t *testing.T) {
	for _, c := range []struct {
		k, n int
		why  string
	}{
		{4, 4, "sample already at K"},
		{8, 4, "sample already below K"},
		{1, 1, "a single lineage"},
		{0, 48, "K below one"},
	} {
		if _, ok := ExpectedCoalescentGens(50, c.k, c.n); ok {
			t.Errorf("k=%d n=%d (%s): expectation should be undefined", c.k, c.n, c.why)
		}
		if _, ok := ImpliedCoalescentNe(100, c.k, c.n); ok {
			t.Errorf("k=%d n=%d (%s): implied Ne should be undefined", c.k, c.n, c.why)
		}
	}
	if _, ok := ExpectedCoalescentGens(0, 4, 48); ok {
		t.Error("Ne=0 should be undefined, not zero generations")
	}
}

// A censored locus has no depth. TMRKAGensAgo must refuse to produce one — averaging a
// censored walk into a mean depth is exactly the collapse TMR4A.md 6 risk 5 forbids.
func TestTMRKAGensAgoRefusesCensoredLoci(t *testing.T) {
	b := newTMRKA([]int{10})
	b.founder(1, -30)
	b.founder(2, -30)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 10, []int{0}, []int{0})
	b.alive(3, 4)

	res := ComputeTMRKA(b.model(2, 100), b.pop, []int{3, 4})
	if res == nil || len(res.Loci) != 1 {
		t.Fatalf("expected one locus, got %+v", res)
	}
	l := res.Loci[0]

	// K=2 is reached at the parents' birth year, 130 years before the sample.
	gens, ok := l.TMRKAGensAgo(2, 26)
	if !ok {
		t.Fatal("K=2 was reached; its depth must be reportable")
	}
	if !closeTo(gens, 5, 1e-9) {
		t.Errorf("TMR-2-A: got %.4f generations, want 5 = 130 years / 26", gens)
	}
	// K=1 is not reached: two unrelated founder strands never coalesce.
	if got, ok := l.TMRKAGensAgo(1, 26); ok {
		t.Errorf("TMR-1-A reported %.2f generations on a censored locus", got)
	}
	// A generation time of zero is not a licence to divide.
	if _, ok := l.TMRKAGensAgo(2, 0); ok {
		t.Error("a zero generation time must be refused, not divided by")
	}
}
