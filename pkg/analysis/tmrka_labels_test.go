package analysis

import "testing"

// W3 in the walker: distinct CREATED alleles surviving, and the created-versus-accumulated
// pairwise decomposition. Every case here is arithmetic the test states independently, not
// a re-derivation of what the walker did.

// labels attaches created-allele labels to the builder's founders.
func (b *tmrkaBuilder) labels(alleles int, founderIDs ...int) {
	b.pop.AssignFounderLabels(founderIDs, alleles)
}

// The case W3 exists for: two founder strands carrying THE SAME created allele. The
// lineages never coalesce — there is no year at which they meet — but they are identical
// by state at founding, so any difference between their descendants is mutational. By
// descent the sample holds 2 founder haplotypes; by created identity it holds 1.
func TestTMRKACreatedAllelesBelowFounderStrands(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, -30) // dad
	b.founder(2, -30) // mom
	// A=1: every founder strand carries the same created allele.
	b.labels(1, 1, 2)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 5, []int{0}, []int{0})
	b.alive(3, 4)

	l := ComputeTMRKA(b.model(2, 100), b.pop, []int{3, 4}).Loci[0]
	if l.FounderStrandsSurviving != 2 {
		t.Fatalf("founder haplotypes: got %d, want 2", l.FounderStrandsSurviving)
	}
	if l.CreatedAllelesSurviving != 1 {
		t.Errorf("created alleles surviving: got %d, want 1 — both founder haplotypes carry "+
			"the same created allele", l.CreatedAllelesSurviving)
	}
	if l.CreatedPairs != 0 {
		t.Errorf("created pairs: got %d, want 0 — no pair is separated by created divergence",
			l.CreatedPairs)
	}
	if l.UnlabelledSources != 0 {
		t.Errorf("unlabelled sources: got %d, want 0", l.UnlabelledSources)
	}
}

// The same genealogy with A=4 (the strict two-diploid-parents case): the two surviving
// founder haplotypes now carry DIFFERENT created alleles, so every cross pair is created
// divergence. Same walk, same coalescences — only the identity model changed.
func TestTMRKACreatedAllelesDistinctUnderA4(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, -30)
	b.founder(2, -30)
	b.labels(4, 1, 2)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 5, []int{0}, []int{0})
	b.alive(3, 4)

	l := ComputeTMRKA(b.model(2, 100), b.pop, []int{3, 4}).Loci[0]
	if l.CreatedAllelesSurviving != 2 {
		t.Errorf("created alleles: got %d, want 2", l.CreatedAllelesSurviving)
	}
	// 4 sampled lineages: 2 descend from dad's strand 0, 2 from mom's strand 0, carrying
	// different created alleles. Pairs = C(4,2) = 6; coalesced = 2·C(2,2) = 2;
	// created = 2·2 = 4.
	if l.CoalescedPairs != 2 {
		t.Errorf("coalesced pairs: got %d, want 2", l.CoalescedPairs)
	}
	if l.CreatedPairs != 4 {
		t.Errorf("created pairs: got %d, want 4 (2 lineages × 2 lineages across the two "+
			"created alleles)", l.CreatedPairs)
	}
	if l.CoalescedPairs+l.CreatedPairs != 6 {
		t.Errorf("the 6 sampled pairs must be fully accounted for; got %d+%d",
			l.CoalescedPairs, l.CreatedPairs)
	}
}

// The three-way case that separates all three quantities at once: pairs that coalesce,
// pairs at distinct founder haplotypes sharing a created allele, and pairs separated by
// created divergence. Only the last is created divergence, and (total − coalesced) would
// overcount it — which is why CreatedPairs is computed by label, not by subtraction.
func TestTMRKAPairDecompositionSeparatesAllThree(t *testing.T) {
	b := newTMRKA([]int{0})
	// Four founders = 8 strands. A=2 ⇒ labels round-robin 0,1,0,1,0,1,0,1 by strand:
	// founder 1→{0,1}, 2→{0,1}, 3→{0,1}, 4→{0,1}.
	b.founder(1, -50)
	b.founder(2, -50)
	b.founder(3, -50)
	b.founder(4, -50)
	b.labels(2, 1, 2, 3, 4)
	// Two sibs from 1×2, both taking dad's strand 0 (label 0) and mom's strand 0 (label 0):
	// their paternal lineages coalesce in 1, their maternal in 2.
	b.child(5, 1, 2, 0, []int{0}, []int{0})
	b.child(6, 1, 2, 0, []int{0}, []int{0})
	// One child of 3×4 taking dad's strand 1 (label 1) and mom's strand 0 (label 0).
	b.child(7, 3, 4, 0, []int{1}, []int{0})
	b.alive(5, 6, 7)

	l := ComputeTMRKA(b.model(4, 60), b.pop, []int{5, 6, 7}).Loci[0]

	// Sources and weights: (1,s0) w=2 label 0; (2,s0) w=2 label 0; (3,s1) w=1 label 1;
	// (4,s0) w=1 label 0.
	if l.FounderStrandsSurviving != 4 {
		t.Fatalf("founder haplotypes: got %d, want 4", l.FounderStrandsSurviving)
	}
	if l.CreatedAllelesSurviving != 2 {
		t.Fatalf("created alleles: got %d, want 2", l.CreatedAllelesSurviving)
	}
	// 6 lineages ⇒ C(6,2) = 15 pairs. Coalesced = C(2,2)+C(2,2) = 2.
	if l.CoalescedPairs != 2 {
		t.Errorf("coalesced pairs: got %d, want 2", l.CoalescedPairs)
	}
	// Created: label 0 carries weight 2+2+1 = 5, label 1 carries 1 ⇒ 5×1 = 5.
	if l.CreatedPairs != 5 {
		t.Errorf("created pairs: got %d, want 5", l.CreatedPairs)
	}
	// The 8 remaining pairs are distinct founder haplotypes sharing a created allele:
	// mutational divergence with no coalescence. Subtracting coalesced from total would
	// have called all 13 of them created.
	if got := 15 - l.CoalescedPairs - l.CreatedPairs; got != 8 {
		t.Errorf("same-allele non-coalesced pairs: got %d, want 8", got)
	}
}

// Lineage weights must be conserved: every sampled lineage ends under exactly one source.
func TestTMRKAWeightsAreConserved(t *testing.T) {
	b := newTMRKA([]int{0, 1, 2})
	b.founder(1, -40)
	b.founder(2, -40)
	b.founder(5, -40)
	b.labels(3, 1, 2, 5)
	b.child(3, 1, 2, 0, []int{0, 1, 0}, []int{1, 1, 0})
	b.child(4, 1, 2, 0, []int{0, 0, 1}, []int{1, 0, 0})
	b.child(6, 3, 5, 30, []int{0, 1, 1}, []int{0, 0, 1})
	b.child(7, 4, 5, 30, []int{1, 1, 0}, []int{0, 1, 0})
	b.alive(6, 7)

	res := ComputeTMRKA(b.model(4, 60), b.pop, []int{6, 7})
	for _, l := range res.Loci {
		total := l.InitialLineages * (l.InitialLineages - 1) / 2
		if l.CoalescedPairs+l.CreatedPairs > total {
			t.Errorf("locus %d: %d coalesced + %d created exceeds %d total pairs",
				l.Position, l.CoalescedPairs, l.CreatedPairs, total)
		}
		if l.CreatedAllelesSurviving > l.FounderStrandsSurviving {
			t.Errorf("locus %d: %d created alleles from only %d founder haplotypes",
				l.Position, l.CreatedAllelesSurviving, l.FounderStrandsSurviving)
		}
		if l.CreatedAllelesSurviving < 1 {
			t.Errorf("locus %d: %d created alleles surviving", l.Position, l.CreatedAllelesSurviving)
		}
	}
}

// With labelling off the created columns are -1, never 0 — an absent quantity must not be
// mistaken for a measured zero. The by-descent columns are unaffected.
func TestTMRKACreatedColumnsAbsentWithoutLabels(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, -30)
	b.founder(2, -30)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 5, []int{0}, []int{0})
	b.alive(3, 4)

	l := ComputeTMRKA(b.model(2, 100), b.pop, []int{3, 4}).Loci[0]
	if l.CreatedAllelesSurviving != -1 || l.CreatedPairs != -1 || l.UnlabelledSources != -1 {
		t.Errorf("created columns with labelling off: got %d/%d/%d, want -1/-1/-1",
			l.CreatedAllelesSurviving, l.CreatedPairs, l.UnlabelledSources)
	}
	// CoalescedPairs is a by-descent quantity and is computed regardless.
	if l.CoalescedPairs != 2 {
		t.Errorf("coalesced pairs: got %d, want 2", l.CoalescedPairs)
	}
	if l.FounderStrandsSurviving != 2 {
		t.Errorf("founder haplotypes: got %d, want 2", l.FounderStrandsSurviving)
	}
}

// A lineage freezing at an UNLABELLED founder must be counted as its own unknown created
// allele, never folded in with a labelled one, and the fact must be reported.
func TestTMRKAUnlabelledSourceIsNotAssumedShared(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, -30)
	b.founder(2, -30)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 5, []int{1}, []int{1})
	b.alive(3, 4)
	// Label founder 1 only: founder 2's strands are unrecorded.
	b.labels(1, 1)

	l := ComputeTMRKA(b.model(4, 100), b.pop, []int{3, 4}).Loci[0]
	if l.FounderStrandsSurviving != 4 {
		t.Fatalf("founder haplotypes: got %d, want 4", l.FounderStrandsSurviving)
	}
	if l.UnlabelledSources != 2 {
		t.Errorf("unlabelled sources: got %d, want 2 (both of founder 2's strands)",
			l.UnlabelledSources)
	}
	// Founder 1's two strands share created allele 0 (A=1); founder 2's two strands are
	// two distinct unknowns ⇒ 3 distinct created alleles, not 1.
	if l.CreatedAllelesSurviving != 3 {
		t.Errorf("created alleles: got %d, want 3 — unknown alleles must not be assumed shared",
			l.CreatedAllelesSurviving)
	}
}
