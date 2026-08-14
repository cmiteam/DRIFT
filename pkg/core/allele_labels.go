package core

import "sort"

// Founder allele identity labels — TMR4A.md work item W3 (`track_allele_labels`).
//
// WHY A LABEL AND NOT A BIT. The founder bitfield holds one bit per site: two states, no
// identity. Two created alleles that happen to carry the same bit at a site are
// indistinguishable from two copies of one allele. That is the identity-by-STATE versus
// identity-by-DESCENT confusion, and the whole TMR4A argument turns on telling them
// apart (TMR4A.md §0): the walker's FounderStrandsSurviving counts distinct founder
// haplotypes BY DESCENT, which is an upper bound on the number of distinct CREATED
// ALLELES the sample actually carries. A label makes the second number exact.
//
// WHAT A LABEL IS. Each founder strand is assigned one of *A* created alleles, where
// *A* = `founder_alleles_per_locus`. The label is the created allele's identity, so two
// founder strands sharing a label carry THE SAME created sequence — they are distinct by
// descent but identical by state at founding, and any difference between their
// descendants is necessarily mutational. Two strands with different labels carry created
// divergence that no amount of elapsed time produced, which is precisely the quantity
// §0 says a coalescent model silently converts into apparent age.
//
// WHY LABELS ARE A PROPERTY OF THE FOUNDER STRAND, NOT OF EVERY INDIVIDUAL. TMR4A.md's
// W3 sketch has the label "inherited through meiosis alongside the bitfield", i.e. stored
// per individual per locus. That is unnecessary once W2 exists: the walk already resolves
// (sampled individual, strand, locus) → (founder individual, founder strand), so the
// label is a lookup at the end of the walk rather than a payload carried through every
// birth. The result is identical and the cost is a 2F-entry table instead of a per-birth
// allocation. This holds because, under the diploid founder model, a founder strand IS
// one created haplotype at every locus. W4's created-gamete model breaks that — a founder
// there transmits more than two alleles and the mapping becomes per-locus — which is why
// FounderLabel takes a locus index it currently ignores: the call sites are already
// correct for W4.
//
// ASSIGNMENT. Deterministic round-robin over founder strands in sorted-id order: global
// strand index g = 2·rank + strand, label = g mod A. No RNG (the §2 gate), balanced (each
// created allele appears on ⌊2F/A⌋ or ⌈2F/A⌉ founder strands), and every one of the A
// alleles is guaranteed present. With A ≥ 2F the assignment is one label per strand and
// created-allele identity collapses onto by-descent identity — which is the honest answer
// for the strict two-diploid-parents case, where 4 strands carry 4 created alleles and the
// two counts genuinely coincide. Modelling A > 4 for a single created couple is exactly
// what W4 is for; it cannot be expressed here.
//
// WHICH GENERATION GETS LABELLED. The pedigree's founder generation — the individuals
// present at setup, which is where the backward walk terminates. If a model seeds genomes
// at a `seed_year` after t=0 the two do not coincide; for the TMR4A study use seed_year=0
// so the pedigree founders ARE the created generation.
//
// NEVER KEY IDENTITY ON POSITION. TMR4A.md W3 and §6 risk 4: `GenerateNewMutations` draws
// `RandIntn(genome_bits)` and positions recur, so a coordinate does not identify an
// allele. Founder alleles are identified by these labels and de-novo mutations by their
// pool IDs; nothing in this module or in the walker keys on a genome coordinate.

// AssignFounderLabels gives every member of the initial population a created-allele label
// per strand. Called once, after population setup. Consumes no RNG.
//
// alleles is *A*; values below 1 are treated as 1 (a single created allele shared by
// every founder strand, i.e. no created divergence anywhere — a useful null arm for the
// study, since it isolates mutational divergence).
func (p *Pop) AssignFounderLabels(founderIDs []int, alleles int) {
	if alleles < 1 {
		alleles = 1
	}
	ids := append([]int(nil), founderIDs...)
	sort.Ints(ids) // map order is randomized; the assignment must be reproducible
	if p.FounderLabels == nil {
		p.FounderLabels = make(map[int][2]int, len(ids))
	}
	for rank, id := range ids {
		p.FounderLabels[id] = [2]int{
			(2 * rank) % alleles,
			(2*rank + 1) % alleles,
		}
	}
}

// FounderLabel returns the created-allele label carried by (ind, strand) at the given
// focal locus, and whether one is recorded. ok=false means labelling is off, or the
// lineage terminated somewhere other than a labelled founder — the caller must treat that
// as an unknown allele rather than defaulting it to a shared one, since assuming shared
// identity is the exact error this module exists to prevent.
//
// locusIdx is accepted and currently ignored: under the diploid founder model a strand
// carries one created haplotype at every locus. W4 makes the mapping locus-dependent and
// will use it.
func (p *Pop) FounderLabel(ind, strand, locusIdx int) (int, bool) {
	if p == nil || strand < 0 || strand > 1 {
		return 0, false
	}
	// Under the created model (W4) the labels are not a separate table: a founder's
	// somatic genome IS the first two alleles of its created pool, so the pool is the
	// authority and consulting it keeps the two items from disagreeing about what
	// founder strand s carries. A one-allele pool means a founder homozygous for that
	// created allele, so both strands answer with it.
	if pool := p.CreatedPool(ind); len(pool) > 0 {
		if strand < len(pool) {
			return pool[strand], true
		}
		return pool[0], true
	}
	if p.FounderLabels == nil {
		return 0, false
	}
	l, ok := p.FounderLabels[ind]
	if !ok {
		return 0, false
	}
	return l[strand], true
}

// FounderLabelsOn reports whether founder-allele identity is available at all — either
// from the W3 label table or from a W4 created-allele pool.
func (p *Pop) FounderLabelsOn() bool {
	return p != nil && (p.FounderLabels != nil || p.FounderAllelePool != nil)
}
