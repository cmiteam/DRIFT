package core

import "testing"

// Created founder alleles (TMR4A.md W4), core half. The pool partition and the provenance
// record; the transmission itself is tested in pkg/simulation, where meiosis lives.

func TestCreatedModelOffByDefault(t *testing.T) {
	m := &Model{Parameters: map[string]float64{}, StringParams: map[string]string{}}
	if CreatedModelOn(m) {
		t.Error("created model must be off unless founder_allele_model=created")
	}
	m.StringParams["founder_allele_model"] = "diploid"
	if CreatedModelOn(m) {
		t.Error("founder_allele_model=diploid must be off")
	}
	m.StringParams["founder_allele_model"] = "created"
	if !CreatedModelOn(m) {
		t.Error("founder_allele_model=created must be on")
	}
	if CreatedModelOn(nil) {
		t.Error("nil model must be off")
	}
}

// A = 2F reproduces ordinary diploidy: two alleles per founder, all distinct. This is the
// strict two-diploid-parents case TMR4A's "four" assumes, and the created model must
// contain it rather than replace it.
func TestCreatedPoolsReduceToDiploidAtA2F(t *testing.T) {
	p := &Pop{}
	p.AssignCreatedPools([]int{0, 1}, 4) // one couple, A=4
	seen := map[int]bool{}
	for _, id := range []int{0, 1} {
		pool := p.CreatedPool(id)
		if len(pool) != 2 {
			t.Fatalf("founder %d pool %v, want 2 alleles", id, pool)
		}
		for _, a := range pool {
			if seen[a] {
				t.Errorf("allele %d appears in two pools; A=2F must partition", a)
			}
			seen[a] = true
		}
	}
	if len(seen) != 4 {
		t.Errorf("got %d created alleles over the couple, want 4", len(seen))
	}
}

// A > 2F is the item's whole point: ONE couple carrying more than four alleles between
// them, which strict diploidy cannot express.
func TestCreatedPoolsExceedFourForOneCouple(t *testing.T) {
	p := &Pop{}
	p.AssignCreatedPools([]int{0, 1}, 10)
	total := 0
	seen := map[int]bool{}
	for _, id := range []int{0, 1} {
		pool := p.CreatedPool(id)
		if len(pool) != 5 {
			t.Errorf("founder %d pool %v, want 5 alleles", id, pool)
		}
		total += len(pool)
		for _, a := range pool {
			if seen[a] {
				t.Errorf("allele %d shared between founders; the partition must be disjoint", a)
			}
			seen[a] = true
		}
	}
	if total != 10 {
		t.Errorf("the couple carries %d created alleles, want 10 — more than four is the point", total)
	}
}

// Order-independent and never empty, including the degenerate A < F case.
func TestCreatedPoolsRobustness(t *testing.T) {
	a, b := &Pop{}, &Pop{}
	a.AssignCreatedPools([]int{7, 2, 5}, 7)
	b.AssignCreatedPools([]int{5, 7, 2}, 7)
	for id, pa := range a.FounderAllelePool {
		pb := b.FounderAllelePool[id]
		if len(pa) != len(pb) {
			t.Fatalf("founder %d: %v vs %v — partition depends on caller order", id, pa, pb)
		}
		for i := range pa {
			if pa[i] != pb[i] {
				t.Fatalf("founder %d: %v vs %v", id, pa, pb)
			}
		}
	}
	// Fewer alleles than founders: everyone must still be able to reproduce.
	c := &Pop{}
	c.AssignCreatedPools([]int{0, 1, 2, 3}, 2)
	for id := 0; id < 4; id++ {
		if len(c.CreatedPool(id)) == 0 {
			t.Errorf("founder %d has an empty pool and cannot reproduce", id)
		}
	}
	// A below 1 is clamped rather than producing an empty allele table.
	d := &Pop{}
	d.AssignCreatedPools([]int{0, 1}, 0)
	for id := 0; id < 2; id++ {
		if len(d.CreatedPool(id)) == 0 {
			t.Errorf("founder %d has an empty pool under A=0", id)
		}
	}
}

// Provenance records the created allele the mask actually selected, using the same
// polarity as meiosis: a SET mask bit takes the FIRST of the germ cell's two alleles.
func TestCreatedCaptureFollowsMaskPolarity(t *testing.T) {
	loci := []int{0, 5, 64}
	mask := make([]uint64, 2)
	mask[0] |= 1 << 5 // locus 5 takes germ[0]; loci 0 and 64 take germ[1]
	p := &Pop{}
	p.CreatedCapture(1, loci,
		[2][]uint64{mask, mask},
		[2][2]int{{7, 3}, {2, 9}},
		[2]bool{true, true})

	want := [2][3]int{{3, 7, 3}, {9, 2, 9}}
	for g := 0; g < 2; g++ {
		for i := range loci {
			got, ok := p.CreatedProvenance(1, g, i, len(loci))
			if !ok {
				t.Fatalf("gamete %d locus %d: no created provenance", g, loci[i])
			}
			if got != want[g][i] {
				t.Errorf("gamete %d locus %d: got allele %d, want %d", g, loci[i], got, want[g][i])
			}
		}
	}
}

// A gamete that was NOT a created draw must report so, distinctly from "created allele 0"
// — otherwise a mixed birth would silently attribute an ordinary strand to a created one.
func TestCreatedProvenanceDistinguishesAbsenceFromAlleleZero(t *testing.T) {
	loci := []int{0, 1}
	mask := make([]uint64, 1)
	p := &Pop{}
	// Paternal gamete is a created draw of allele 0; maternal gamete is ordinary.
	p.CreatedCapture(1, loci, [2][]uint64{mask, mask},
		[2][2]int{{0, 0}, {0, 0}}, [2]bool{true, false})

	for i := range loci {
		got, ok := p.CreatedProvenance(1, 0, i, len(loci))
		if !ok || got != 0 {
			t.Errorf("paternal locus %d: got (%d,%v), want (0,true)", i, got, ok)
		}
		if _, ok := p.CreatedProvenance(1, 1, i, len(loci)); ok {
			t.Errorf("maternal locus %d reported a created draw where there was none", i)
		}
	}
	// A birth with no created gamete at all stores nothing.
	q := &Pop{}
	q.CreatedCapture(2, loci, [2][]uint64{mask, mask}, [2][2]int{}, [2]bool{false, false})
	if q.CreatedARGDB != nil {
		t.Error("stored a created record for a birth with no created gametes")
	}
	if _, ok := q.CreatedProvenance(2, 0, 0, 2); ok {
		t.Error("unrecorded individual reported a created draw")
	}
}

// Under the created model the pool is the authority on what a founder's somatic strands
// carry, so W3 and W4 cannot disagree about founder identity.
func TestFounderLabelPrefersCreatedPool(t *testing.T) {
	p := &Pop{}
	p.AssignFounderLabels([]int{0, 1}, 4) // W3 table: 0→{0,1}, 1→{2,3}
	p.AssignCreatedPools([]int{0, 1}, 10) // W4 pools: 0→{0,2,4,6,8}, 1→{1,3,5,7,9}

	if l, _ := p.FounderLabel(0, 0, 0); l != 0 {
		t.Errorf("founder 0 strand 0: got %d, want 0 (pool[0])", l)
	}
	if l, _ := p.FounderLabel(0, 1, 0); l != 2 {
		t.Errorf("founder 0 strand 1: got %d, want 2 (pool[1]) — the pool must win over the W3 table", l)
	}
	if l, _ := p.FounderLabel(1, 1, 0); l != 3 {
		t.Errorf("founder 1 strand 1: got %d, want 3", l)
	}
	// A founder homozygous for a single created allele answers with it on both strands.
	q := &Pop{}
	q.AssignCreatedPools([]int{0, 1, 2, 3}, 2)
	pool := q.CreatedPool(0)
	for s := 0; s < 2; s++ {
		if l, ok := q.FounderLabel(0, s, 0); !ok || l != pool[0] {
			t.Errorf("one-allele pool, strand %d: got (%d,%v), want (%d,true)", s, l, ok, pool[0])
		}
	}
}

// Created provenance shares the pedigree's lifetime, like the W2 record it parallels.
func TestCreatedProvenancePrunedWithPedigree(t *testing.T) {
	p := &Pop{}
	p.PedEnsure(1, 0, 0)
	p.PedEnsure(2, 0, 1)
	p.AssignCreatedPools([]int{1, 2}, 6)
	p.PedRecordBirth(3, 1, 2, 10, 0)
	p.CreatedCapture(3, []int{0}, [2][]uint64{make([]uint64, 1), make([]uint64, 1)},
		[2][2]int{{1, 3}, {2, 4}}, [2]bool{true, true})

	p.PedRecordDeath(1)
	p.PedRecordDeath(2)
	if _, ok := p.CreatedARGDB[3]; !ok {
		t.Fatal("3 is alive; its created provenance must be retained")
	}
	p.PedRecordDeath(3)
	if _, ok := p.CreatedARGDB[3]; ok {
		t.Error("created provenance outlived the pedigree node it indexes")
	}
	// Pools are the standing definition of what was created and are NOT pruned.
	if len(p.FounderAllelePool) != 2 {
		t.Errorf("founder allele pools were pruned: %v", p.FounderAllelePool)
	}
}
