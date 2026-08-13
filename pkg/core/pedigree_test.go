package core

import (
	"sort"
	"testing"
)

// buildPed constructs a small pedigree by replaying births in order. Founders are
// ensured first, exactly as InitializePop does for the initial population.
func buildPed(founders []int, births [][5]int) *Pop {
	p := &Pop{}
	for _, f := range founders {
		p.PedEnsure(f, 0, f%2)
	}
	for _, b := range births {
		// b = {child, dad, mom, birthYear, sex}
		p.PedRecordBirth(b[0], b[1], b[2], b[3], b[4])
	}
	return p
}

func ids(p *Pop) []int {
	out := make([]int, 0, len(p.Pedigree))
	for id := range p.Pedigree {
		out = append(out, id)
	}
	sort.Ints(out)
	return out
}

// A pedigree records BOTH parents of every child — the gap W1 exists to close. The old
// MaleDB/FemaleDB pair records a son's father and a daughter's mother only, so a son's
// mother and a daughter's father are unrecoverable.
func TestPedigreeRecordsBothParents(t *testing.T) {
	//  1(m) x 2(f) -> 3 (son), 4 (daughter)
	p := buildPed([]int{1, 2}, [][5]int{
		{3, 1, 2, 10, 0},
		{4, 1, 2, 12, 1},
	})
	for _, child := range []int{3, 4} {
		n := p.Pedigree[child]
		if n == nil {
			t.Fatalf("child %d has no pedigree node", child)
		}
		if n.Dad != 1 || n.Mom != 2 {
			t.Errorf("child %d: got dad=%d mom=%d, want 1/2", child, n.Dad, n.Mom)
		}
	}
	if got := p.Pedigree[3].BirthYear; got != 10 {
		t.Errorf("child 3 birth year: got %d, want 10", got)
	}
	if p.Pedigree[1].Dad != PedFounder || p.Pedigree[2].Mom != PedFounder {
		t.Error("founders must have no recorded parents (walk terminates there)")
	}
}

// Refcounts: alive contributes 1, each retained child contributes 1.
func TestPedigreeRefcounts(t *testing.T) {
	p := buildPed([]int{1, 2}, [][5]int{
		{3, 1, 2, 10, 0},
		{4, 1, 2, 12, 1},
	})
	for _, id := range []int{1, 2} {
		if got := p.Pedigree[id].Refs; got != 3 { // alive + 2 children
			t.Errorf("founder %d refs: got %d, want 3", id, got)
		}
	}
	if got := p.Pedigree[3].Refs; got != 1 { // alive, childless
		t.Errorf("child 3 refs: got %d, want 1", got)
	}
}

// A dead individual who left no descendants is dropped; one with a living descendant is
// retained, along with every ancestor on the path — including ancestors the single-sex
// MaleDB/FemaleDB chains would have dropped.
func TestPedigreePrunesOnlyExtinctBranches(t *testing.T) {
	// 1(m) x 2(f) -> 3(m), 4(f, dies childless)
	// 3 x 5(f)    -> 6(f) who stays alive
	p := buildPed([]int{1, 2, 5}, [][5]int{
		{3, 1, 2, 10, 0},
		{4, 1, 2, 12, 1},
		{6, 3, 5, 30, 1},
	})

	p.PedRecordDeath(4) // childless: goes immediately
	if _, ok := p.Pedigree[4]; ok {
		t.Error("childless dead individual 4 should have been pruned")
	}

	// Everyone else dies; only 6 is left alive.
	for _, id := range []int{1, 2, 3, 5} {
		p.PedRecordDeath(id)
	}
	got := ids(p)
	want := []int{1, 2, 3, 5, 6} // 6's whole ancestry is retained
	if len(got) != len(want) {
		t.Fatalf("retained %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("retained %v, want %v", got, want)
		}
	}

	// 6 dies childless: the entire tree becomes unreachable and the cascade clears it.
	p.PedRecordDeath(6)
	if len(p.Pedigree) != 0 {
		t.Errorf("cascade left %v behind; the whole branch is extinct", ids(p))
	}
}

// The refcount discipline must be symmetric: every increment has exactly one matching
// decrement. This is the §6f mutation-pool failure mode restated — replay a whole
// birth/death history and assert the map ends empty once nobody is alive.
func TestPedigreeRefcountBalance(t *testing.T) {
	p := &Pop{}
	const nFounders = 6
	for i := 0; i < nFounders; i++ {
		p.PedEnsure(i, -i, i%2)
	}
	// A deterministic braided history: each new child takes an earlier male and an
	// earlier female, so lineages cross and nodes accumulate multiple references.
	next := nFounders
	born := []int{}
	for gen := 0; gen < 8; gen++ {
		for k := 0; k < 4; k++ {
			dad := (gen*4 + k) % next
			mom := (gen*4 + k + 3) % next
			if dad == mom {
				mom = (mom + 1) % next
			}
			p.PedRecordBirth(next, dad, mom, gen*20, k%2)
			born = append(born, next)
			next++
		}
	}
	// Kill everyone, oldest first.
	all := append([]int{}, born...)
	for i := 0; i < nFounders; i++ {
		all = append([]int{i}, all...)
	}
	for _, id := range all {
		p.PedRecordDeath(id)
	}
	if len(p.Pedigree) != 0 {
		t.Errorf("refcount leak: %d nodes retained with nobody alive: %v", len(p.Pedigree), ids(p))
	}
}

// A double RIP must not free a node that is still referenced.
func TestPedigreeDoubleDeathIsSafe(t *testing.T) {
	p := buildPed([]int{1, 2}, [][5]int{{3, 1, 2, 10, 0}})
	p.PedRecordDeath(1)
	p.PedRecordDeath(1)
	if _, ok := p.Pedigree[1]; !ok {
		t.Fatal("node 1 must be retained: its living child 3 descends through it")
	}
	if got := p.Pedigree[1].Refs; got != 1 {
		t.Errorf("node 1 refs after double death: got %d, want 1 (the child)", got)
	}
}

// Off by default: with no pedigree map, every method is inert and nothing is allocated.
func TestPedigreeOffIsInert(t *testing.T) {
	p := &Pop{}
	if p.PedigreeOn() {
		t.Error("pedigree must be off until something writes to it")
	}
	p.PedRecordDeath(7) // must not panic or allocate
	if p.Pedigree != nil {
		t.Error("PedRecordDeath allocated a pedigree map on the disabled path")
	}
	if got := p.PedAncestryPath(7, false); got != nil {
		t.Errorf("PedAncestryPath on a disabled pedigree: got %v, want nil", got)
	}
}

// The walk terminates at a founder rather than running off into unrecorded ancestry.
func TestPedigreeAncestryPathTerminatesAtFounder(t *testing.T) {
	p := buildPed([]int{1, 2, 5}, [][5]int{
		{3, 1, 2, 10, 0},
		{6, 3, 5, 30, 1},
	})
	got := p.PedAncestryPath(6, false) // paternal chain
	want := []int{6, 3, 1}
	if len(got) != len(want) {
		t.Fatalf("paternal path: got %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("paternal path: got %v, want %v", got, want)
		}
	}
}
