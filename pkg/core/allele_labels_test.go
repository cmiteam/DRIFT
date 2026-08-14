package core

import "testing"

// Founder allele labels (TMR4A.md W3). The property under test is not "labels round-trip"
// but the one the argument rests on: labels must let identity-by-STATE come apart from
// identity-by-DESCENT, and must never let two founder strands be assumed to share a
// created allele when they do not.

func TestAssignFounderLabelsRoundRobin(t *testing.T) {
	p := &Pop{}
	// Founders given out of order on purpose: the assignment must depend on sorted id,
	// not on the caller's slice order or on map iteration.
	p.AssignFounderLabels([]int{3, 1, 0, 2}, 4)

	// ranks 0..3 → global strand indices 0..7 → labels 0..3, 0..3
	want := map[int][2]int{
		0: {0, 1},
		1: {2, 3},
		2: {0, 1},
		3: {2, 3},
	}
	for id, w := range want {
		got, ok := p.FounderLabels[id]
		if !ok {
			t.Fatalf("founder %d unlabelled", id)
		}
		if got != w {
			t.Errorf("founder %d: got labels %v, want %v", id, got, w)
		}
	}
	if !p.FounderLabelsOn() {
		t.Error("labelling should report as on")
	}
}

func TestAssignFounderLabelsOrderIndependent(t *testing.T) {
	a, b := &Pop{}, &Pop{}
	a.AssignFounderLabels([]int{5, 9, 2, 7}, 3)
	b.AssignFounderLabels([]int{9, 7, 5, 2}, 3)
	for id, la := range a.FounderLabels {
		if lb := b.FounderLabels[id]; lb != la {
			t.Errorf("founder %d: %v vs %v — assignment depends on caller order", id, la, lb)
		}
	}
}

// A single created couple with A=4 gives one created allele per strand: identity by state
// and identity by descent coincide, which is the honest answer for the strict
// two-diploid-parents case that TMR4A's "four" assumes.
func TestFounderLabelsStrictDiploidCoupleIsOneToOne(t *testing.T) {
	p := &Pop{}
	p.AssignFounderLabels([]int{0, 1}, 4)
	seen := map[int]bool{}
	for _, id := range []int{0, 1} {
		for s := 0; s < 2; s++ {
			label, ok := p.FounderLabel(id, s, 0)
			if !ok {
				t.Fatalf("founder %d strand %d unlabelled", id, s)
			}
			if seen[label] {
				t.Errorf("label %d reused: A=4 over 4 founder strands must be one-to-one", label)
			}
			seen[label] = true
		}
	}
	if len(seen) != 4 {
		t.Errorf("got %d distinct created alleles, want 4", len(seen))
	}
}

// A < 2F is where the two identities come apart: several founder strands carry the SAME
// created allele, so their descendants' differences are mutational despite never
// coalescing. This is the case the label exists for.
func TestFounderLabelsFewerAllelesThanStrands(t *testing.T) {
	p := &Pop{}
	p.AssignFounderLabels([]int{0, 1, 2, 3, 4}, 2) // 10 strands, 2 created alleles
	counts := map[int]int{}
	for id := 0; id < 5; id++ {
		for s := 0; s < 2; s++ {
			label, ok := p.FounderLabel(id, s, 0)
			if !ok {
				t.Fatalf("founder %d strand %d unlabelled", id, s)
			}
			counts[label]++
		}
	}
	if len(counts) != 2 {
		t.Fatalf("got %d distinct created alleles, want 2", len(counts))
	}
	for label, n := range counts {
		if n != 5 {
			t.Errorf("created allele %d carried by %d founder strands, want 5 (balanced)", label, n)
		}
	}
}

// A=1 is the null arm: every founder strand carries the same created allele, so ALL
// divergence in the run is mutational. Worth having as a control for the study.
func TestFounderLabelsSingleAlleleNullArm(t *testing.T) {
	p := &Pop{}
	p.AssignFounderLabels([]int{0, 1, 2}, 1)
	for id := 0; id < 3; id++ {
		for s := 0; s < 2; s++ {
			if label, _ := p.FounderLabel(id, s, 0); label != 0 {
				t.Errorf("founder %d strand %d: label %d, want 0", id, s, label)
			}
		}
	}
	// A below 1 is clamped rather than producing a divide-by-zero.
	q := &Pop{}
	q.AssignFounderLabels([]int{0}, 0)
	if label, ok := q.FounderLabel(0, 0, 0); !ok || label != 0 {
		t.Errorf("A=0 clamp: got (%d,%v), want (0,true)", label, ok)
	}
}

// An unrecorded individual must report NO label rather than a default one. Defaulting
// would make two unrelated strands look like the same created allele — the exact error
// this module exists to prevent.
func TestFounderLabelUnknownIsNotDefaulted(t *testing.T) {
	p := &Pop{}
	if _, ok := p.FounderLabel(1, 0, 0); ok {
		t.Error("labelling off must report no label")
	}
	if p.FounderLabelsOn() {
		t.Error("labelling should be off until something is assigned")
	}
	p.AssignFounderLabels([]int{1}, 4)
	if _, ok := p.FounderLabel(99, 0, 0); ok {
		t.Error("an unlabelled individual must report no label")
	}
	if _, ok := p.FounderLabel(1, 2, 0); ok {
		t.Error("an out-of-range strand must report no label")
	}
}

// Labels share the pedigree's lifetime: a released node takes its label with it.
func TestFounderLabelsPrunedWithPedigree(t *testing.T) {
	p := &Pop{}
	p.PedEnsure(1, 0, 0)
	p.PedEnsure(2, 0, 1)
	p.AssignFounderLabels([]int{1, 2}, 4)
	p.PedRecordBirth(3, 1, 2, 10, 0)

	p.PedRecordDeath(1)
	if _, ok := p.FounderLabels[1]; !ok {
		t.Fatal("1 has a living descendant; its label must be retained")
	}
	p.PedRecordDeath(2)
	p.PedRecordDeath(3) // whole branch extinct
	if len(p.FounderLabels) != 0 {
		t.Errorf("labels outlived the pedigree: %v", p.FounderLabels)
	}
}
