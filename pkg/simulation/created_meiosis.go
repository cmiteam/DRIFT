package simulation

import (
	"drift/pkg/core"
	"drift/pkg/utils"
)

// Created-founder germline transmission — TMR4A.md work item W4, the meiosis half.
//
// THIS FILE IS SEPARATE ON PURPOSE. TMR4A.md §6 risk 2: "W4 has no precedent in the
// codebase and touches meiosis, the most correctness-critical path in the engine. Gate it
// hard; keep the created-gamete path structurally separate." So the created path lives
// here rather than as a branch inside birthStandard, and birth.go calls into it at two
// clearly-marked points.
//
// WHAT IT CHANGES, AND WHAT IT DELIBERATELY DOES NOT. Under the created model a founder's
// germ cell carries TWO alleles drawn from that founder's created pool, instead of the
// founder's two somatic strands. Everything downstream is untouched: the same crossover
// mask, the same recombination model, the same meiosis polarity. The gamete is still a
// two-parent-strand mosaic — it is only the identity of those two strands that varies
// from one germ cell to the next, which is what lets a single couple transmit more than
// four alleles at a locus while every other genetic mechanism behaves as before.
//
// RNG. Exactly two draws per created gamete (the germ cell's two alleles), taken in a
// fixed order — paternal gamete then maternal — regardless of how many focal loci are
// configured. Drawing per locus instead would make the stream depend on `arg_loci`, which
// would destroy the locus-sweep reproducibility W8 tests. Nothing here executes unless
// founder_allele_model="created" AND the parent has a created pool, so the default path
// draws nothing and is byte-identical.

// createdGermCell draws the two created alleles a founder's germ cell carries for this
// gamete. Returns ok=false when `parent` has no created pool — it is not a founder, or
// the created seeder never ran — in which case the caller takes the ordinary path.
//
// The two alleles are drawn independently and with replacement: a germ cell may carry the
// same created allele twice, which is the founder being homozygous for it in that cell.
// Excluding that would forbid a real state rather than model one.
func createdGermCell(model *core.Model, pop *core.Pop, parent int) ([2]int, bool) {
	if !core.CreatedModelOn(model) {
		return [2]int{}, false
	}
	pool := pop.CreatedPool(parent)
	if len(pool) == 0 {
		return [2]int{}, false
	}
	if len(pool) == 1 {
		// A single-allele pool is deterministic; drawing would consume RNG to no effect.
		return [2]int{pool[0], pool[0]}, true
	}
	return [2]int{pool[utils.RandIntn(len(pool))], pool[utils.RandIntn(len(pool))]}, true
}

// createdMeiosis writes the child's gamete from a germ cell's two created alleles, using
// the mask exactly as ordinary meiosis does. Delegates to meiosisFrom so the polarity —
// a SET mask bit takes the FIRST source — is defined in one place for both paths and
// cannot drift between them.
func createdMeiosis(pop *core.Pop, mask []uint64, germ [2]int, child, copyIdx int) {
	seqs := pop.CreatedAlleleSeqs
	if germ[0] >= len(seqs) || germ[1] >= len(seqs) {
		return // a pool referencing an allele that was never built; leave the gamete empty
	}
	pop.Chromosomes[child][copyIdx] = meiosisFrom(mask, seqs[germ[0]], seqs[germ[1]])
}
