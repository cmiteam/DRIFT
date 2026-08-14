package analysis

// Per-lineage terminus resolution — the piece the W6 merged VCF export needs from the
// W5 walker, and the only thing in this package that answers "where did THIS strand
// come from" rather than "how many lineages were there at time t".
//
// WHY IT IS NOT THE WALKER. walkLocus (tmrka.go) merges lineages as they coalesce and
// carries a sample WEIGHT through the merge, which is exactly right for the pairwise
// decomposition and destroys the mapping from an individual sampled strand back to its
// source. The VCF annotation needs that mapping: one FL value per sample per strand.
// So this walks each lineage independently to its terminus, using the SAME elementary
// move (stepBack) and the SAME terminus rule (sourceLabel) as the walker, which is what
// keeps the two from drifting apart. It is deliberately NOT a second implementation of
// the coalescence logic — there is none here to get wrong.
//
// COST. O(2n · depth) per focal locus with no heap, against the walker's
// O(steps · log n). Focal locus sets are small (TMR4A.md §1A), so this is cheap; it
// consumes no RNG and reads nothing but Pedigree + ARGDB + CreatedARGDB.

import "drift/pkg/core"

// FounderAlleleLabelsAt resolves every sampled strand back to the created-allele label
// it carries at focal locus index locusIdx, by walking the recorded genealogy to its
// terminus. The result maps individual id -> {strand 0 label, strand 1 label}, with -1
// where no label is recoverable.
//
// -1 IS NOT A DEFAULT, IT IS AN ADMISSION. A strand whose lineage freezes somewhere
// other than a labelled founder carries an UNKNOWN created allele, and the whole point
// of W3 is that assuming shared identity in the absence of evidence is the error to
// avoid. Callers must treat two -1 strands as possibly different alleles, never as the
// same one; the walker's UnlabelledSources counts exactly this class and is 0 in a
// well-formed run.
//
// Returns nil when there is no substrate to walk.
func FounderAlleleLabelsAt(pop *core.Pop, ids []int, locusIdx, nLoci int) map[int][2]int {
	if pop == nil || pop.Pedigree == nil || locusIdx < 0 || locusIdx >= nLoci {
		return nil
	}
	// The pedigree is a DAG in birth-year order so every walk terminates; the cap only
	// keeps a corrupted pedigree from turning that into a hang. Same discipline as
	// walkLocus's maxSteps.
	maxSteps := 64 * (len(pop.Pedigree) + 1)

	out := make(map[int][2]int, len(ids))
	for _, id := range ids {
		var pair [2]int
		for s := 0; s < 2; s++ {
			key := lineageKey{ind: id, strand: s}
			for steps := 0; steps < maxSteps; steps++ {
				next, ok := stepBack(pop, key, locusIdx, nLoci)
				if !ok {
					break // founder, created allele, or unrecorded provenance
				}
				if _, pok := pop.Pedigree[next.ind]; !pok {
					break // the parent left no node: this is as far back as truth goes
				}
				key = next
			}
			pair[s] = -1
			if label, ok := sourceLabel(pop, key, locusIdx); ok {
				pair[s] = label
			}
		}
		out[id] = pair
	}
	return out
}
