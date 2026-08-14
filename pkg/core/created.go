package core

import "sort"

// Created founder alleles with germline heterogeneity — TMR4A.md work item W4
// (`founder_allele_model=created`).
//
// THE OBJECTION THIS ENCODES. TMR4A is the time at which a locus has four lineages, and
// reading that as a bound on a founding couple assumes four is the ceiling — which holds
// only if each founder's germline is a faithful copy of a single diploid genome. If the
// REPRODUCTIVE CELLS of the created couple carried variation among themselves, the pair
// transmits more than four alleles per locus and "4" is the wrong bound to test
// (TMR4A.md §0). DRIFT could not represent that: founders are strictly diploid and
// meiosis draws from exactly two strands. This module is what lets it.
//
// THE MODEL, AND WHY IT IS THIS ONE. A created allele is a whole HAPLOTYPE — a sequence
// over the genome, not a value at a point — so *A* created alleles are *A* haplotypes,
// generated once at seeding, any two differing at an expected fraction
// `founder_allele_divergence` of sites. Each founder carries a POOL of them. Then:
//
//	a founder's SOMATIC genome is two of its pool alleles — it is an ordinary diploid;
//	a founder's GERMLINE is heterogeneous — each germ cell carries two alleles drawn
//	from the pool, and ordinary meiosis recombines that cell's two.
//
// That is "created germline heterogeneity, not somatic diploidy" taken literally. It also
// makes the change small and safe where it must be: recombination, the mask, and the
// meiosis polarity are all untouched — the ONLY difference is which two sequences the
// mask is applied to. A founder can therefore transmit every allele in its pool across
// its offspring while every other genetic mechanism behaves exactly as it always did.
//
// The alternative — drawing an allele independently at each locus — was rejected: it
// destroys linkage within a gamete, so the transmitted haplotype would be a mosaic of
// arbitrarily many created alleles and no recombination model would apply to it. A germ
// cell carrying two alleles is both more physical and less invasive.
//
// COALESCENCE AT A CREATED ALLELE. Two lineages that reach the SAME created allele in the
// same founder have met: they are copies of one created sequence, taken from one
// germline, with no mutational difference between them. They coalesce, dated to that
// founder. Two lineages at DIFFERENT created alleles never coalesce at all — no
// mutational path connects them and no year exists at which they meet, which is exactly
// the censoring §0 insists must be reportable rather than converted into a date.
//
// GATING. `founder_allele_model` defaults to "diploid", under which nothing here is
// allocated, no RNG is drawn, and every path is the pre-W4 path byte-for-byte. The
// created path is reached only when the model is set AND the created seeder has built a
// pool for the parent in question; drift.go warns when one is configured without the
// other rather than silently taking half the model.

// CreatedNone marks a gamete that was NOT drawn from a created-allele pool, so a reader
// can tell "no created draw here" from "created allele 0".
const CreatedNone = uint16(0xFFFF)

// CreatedModelOn reports whether the created founder-allele model is selected. This is
// the master switch; the presence of a pool for a given founder is the second condition
// every created path also checks.
func CreatedModelOn(m *Model) bool {
	return m != nil && m.StringParam("founder_allele_model", "diploid") == "created"
}

// CreatedPool returns the created-allele indices founder `id` carries, or nil.
func (p *Pop) CreatedPool(id int) []int {
	if p == nil || p.FounderAllelePool == nil {
		return nil
	}
	return p.FounderAllelePool[id]
}

// AssignCreatedPools apportions *A* created alleles among the founders — the "distribution
// across gametes" W4 calls for, in its minimum viable form: a stride partition, so
// founder of rank r carries every allele a with a ≡ r (mod F). Deterministic, RNG-free,
// and balanced. With A = 2F each founder gets exactly two and the model reduces to
// ordinary diploidy; A > 2F is the germline-heterogeneity case the whole item exists for.
// A founder that would otherwise receive nothing (A < F) is given one allele so it can
// still reproduce.
func (p *Pop) AssignCreatedPools(founderIDs []int, alleles int) {
	if alleles < 1 {
		alleles = 1
	}
	ids := append([]int(nil), founderIDs...)
	sort.Ints(ids) // map order is randomized; the partition must be reproducible
	if p.FounderAllelePool == nil {
		p.FounderAllelePool = make(map[int][]int, len(ids))
	}
	f := len(ids)
	for rank, id := range ids {
		var pool []int
		for a := rank; a < alleles; a += f {
			pool = append(pool, a)
		}
		if len(pool) == 0 {
			pool = []int{rank % alleles}
		}
		p.FounderAllelePool[id] = pool
	}
}

// CreatedCapture records, for one birth, which created allele each focal locus was drawn
// from on each gamete whose parent presented a created germ cell (TMR4A.md W4 provenance).
//
// Held in a SEPARATE record from the W2 strand bits on purpose. W2's hot path is one bit
// per locus per gamete and covers every birth in the run; created draws happen only where
// a parent is a founder — one generation in an Eden run — and need a whole allele index.
// Folding the wide, rare case into the narrow, universal one would cost every birth in
// the simulation to serve a handful of them.
//
// germ[g] is the two created alleles the parent's germ cell carried; the mask decides
// which of the two each locus took, using exactly the meiosis convention documented in
// arg.go — a SET mask bit means the gamete took the first of the pair.
func (p *Pop) CreatedCapture(child int, loci []int, masks [2][]uint64, germ [2][2]int, present [2]bool) {
	if len(loci) == 0 || (!present[0] && !present[1]) {
		return
	}
	if p.CreatedARGDB == nil {
		p.CreatedARGDB = make(map[int][]uint16)
	}
	rec := make([]uint16, 2*len(loci))
	for g := 0; g < 2; g++ {
		base := g * len(loci)
		if !present[g] {
			for i := range loci {
				rec[base+i] = CreatedNone
			}
			continue
		}
		mask := masks[g]
		for i, pos := range loci {
			word := pos / 64
			pick := 1
			if word < len(mask) && mask[word]&(1<<(uint(pos)%64)) != 0 {
				pick = 0
			}
			rec[base+i] = uint16(germ[g][pick])
		}
	}
	p.CreatedARGDB[child] = rec
}

// CreatedProvenance returns the created allele individual `child` inherited at focal locus
// `locusIdx` on gamete `gamete`, and whether that gamete was a created draw at all.
// ok=false means the ordinary strand provenance (W2) applies.
func (p *Pop) CreatedProvenance(child, gamete, locusIdx, nLoci int) (int, bool) {
	if p == nil || p.CreatedARGDB == nil || gamete < 0 || gamete > 1 || locusIdx < 0 || locusIdx >= nLoci {
		return 0, false
	}
	rec, ok := p.CreatedARGDB[child]
	if !ok || len(rec) < 2*nLoci {
		return 0, false
	}
	v := rec[gamete*nLoci+locusIdx]
	if v == CreatedNone {
		return 0, false
	}
	return int(v), true
}
