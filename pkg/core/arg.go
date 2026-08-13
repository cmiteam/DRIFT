package core

import (
	"fmt"
	"sort"
	"strconv"
	"strings"
)

// Strand provenance — TMR4A.md work item W2 (`track_arg`, `arg_loci`), focal-loci mode.
//
// WHY THIS EXISTS. A pedigree (W1) says who your parents were. It does not say which of
// a parent's two strands you inherited at a given position, and without that the
// autosomal genealogy at a locus is unrecoverable: a child has two parents but only one
// ancestor per strand per locus, and which one is decided by the meiosis mask. Today
// those masks are locals in birth.go — built, consumed by meiosis and InheritMutations,
// and discarded. This module keeps the answer for a configured set of focal loci.
//
// WHAT IS STORED. Per birth, per gamete, one BIT per focal locus: the index of the
// parental strand that locus was drawn from (0 or 1). Two gametes per birth, so a run
// with |L| focal loci pays 2·⌈|L|/64⌉ words per retained individual — for the handful of
// loci a focal study uses, one word per birth. Contrast the breakpoint mode TMR4A.md
// also specifies (store the crossover list, ~23-40 ints per gamete, which yields EVERY
// locus and is the tskit-edge form §7 wants): richer, and the eventual target, but this
// mode is the cheap one that makes the walker testable first.
//
// THE MASK CONVENTION, IN ONE PLACE. meiosis() computes
//
//	child = (mask & parent[0]) | (^mask & parent[1])
//
// so a SET mask bit means the child took parental copy 0 at that position, and a CLEAR
// mask bit means copy 1. That inversion is the single most error-prone fact in this
// module, so it is applied in exactly one function (ARGCapture) and asserted directly
// against meiosis() in the tests rather than being re-derived at each read site.
//
// GAMETE INDEXING follows birth.go: gamete 0 is the paternal gamete (mask from
// createMask(model,0), written to child strand 0), gamete 1 is the maternal gamete
// (child strand 1). So a lineage sitting on child strand g at locus ℓ steps back to
// (parent, ARGStrand(child,g,ℓ)) where the parent is Dad for g=0 and Mom for g=1 — which
// is exactly the step the W5 walker takes.
//
// LIFETIME. Provenance records are pruned by the SAME refcount cascade as the pedigree
// (see pedRelease): a node the walker can no longer reach takes its provenance with it.
// track_arg therefore requires track_pedigree — without it the records would be both
// unreadable (no parent links) and unbounded (nothing to prune them). birth.go enforces
// this and drift.go warns at startup.
//
// DETERMINISM. Capture is a read of masks that were drawn anyway; it consumes no RNG of
// its own. Enabling track_arg on a run that already tracks DNA or mutations is therefore
// byte-identical. Enabling it on a run that tracks NEITHER does change the stream, because
// masks are only drawn when there is something to inherit — that is opt-in behaviour,
// documented here and covered by TestARGCaptureIsFreeWhenMasksAlreadyDrawn.

// ARGFocalOn reports whether focal-loci strand provenance should be captured.
func ARGFocalOn(m *Model) bool {
	return m != nil && m.Parameters["track_arg"] >= 1 && len(EnsureARGLoci(m)) > 0
}

// EnsureARGLoci parses and caches the focal locus set from the `arg_loci` parameter.
// Built once (no RNG) and cached on the model, like GeneticMap. Returns a sorted,
// deduplicated, in-range slice; the LOCUS INDEX used everywhere else in this module is
// an index into this slice, not a genome-bit position, so the packed bitsets stay dense
// however sparse the positions are.
func EnsureARGLoci(m *Model) []int {
	if m == nil {
		return nil
	}
	if m.ARGLoci != nil {
		return m.ARGLoci
	}
	m.ARGLoci = ParseARGLoci(m.StringParam("arg_loci", ""), m.FreeParameters["genome_bits"])
	return m.ARGLoci
}

// ParseARGLoci turns an `arg_loci` spec into genome-bit positions.
//
// Grammar — comma- or whitespace-separated items, each one of:
//
//	P        a single genome-bit position
//	A-B      every position in [A, B] inclusive
//	A-B:S    every S-th position in [A, B] inclusive (even spacing across an arm)
//
// Positions outside [0, genomeBits) are dropped, the result is sorted and deduplicated,
// and an unparseable item is skipped with a warning rather than aborting the run — a
// typo in a locus list should not throw away a long simulation. Returns a non-nil empty
// slice for an empty spec so the caller's cache is populated and reparsing is skipped.
func ParseARGLoci(spec string, genomeBits int) []int {
	out := []int{}
	fields := strings.FieldsFunc(spec, func(r rune) bool {
		return r == ',' || r == ' ' || r == '\t' || r == ';'
	})
	seen := map[int]bool{}
	add := func(p int) {
		if p < 0 || (genomeBits > 0 && p >= genomeBits) || seen[p] {
			return
		}
		seen[p] = true
		out = append(out, p)
	}

	for _, f := range fields {
		item := strings.TrimSpace(f)
		if item == "" {
			continue
		}
		step := 1
		if i := strings.LastIndex(item, ":"); i >= 0 {
			s, err := strconv.Atoi(strings.TrimSpace(item[i+1:]))
			if err != nil || s < 1 {
				fmt.Printf("WARNING: arg_loci: bad step in %q, skipped\n", item)
				continue
			}
			step = s
			item = strings.TrimSpace(item[:i])
		}
		if i := strings.Index(item, "-"); i > 0 { // i>0: a leading '-' is a bad position, not a range
			lo, errLo := strconv.Atoi(strings.TrimSpace(item[:i]))
			hi, errHi := strconv.Atoi(strings.TrimSpace(item[i+1:]))
			if errLo != nil || errHi != nil || hi < lo {
				fmt.Printf("WARNING: arg_loci: bad range %q, skipped\n", f)
				continue
			}
			for p := lo; p <= hi; p += step {
				add(p)
			}
			continue
		}
		p, err := strconv.Atoi(item)
		if err != nil {
			fmt.Printf("WARNING: arg_loci: bad position %q, skipped\n", f)
			continue
		}
		add(p)
	}

	sort.Ints(out)
	return out
}

// argWords is the number of uint64 words one gamete's provenance bitset needs.
func argWords(nLoci int) int { return (nLoci + 63) / 64 }

// ARGCapture records, for one birth, which parental strand each focal locus was drawn
// from on each gamete. maskDad / maskMom are the meiosis masks birth.go has just built
// and is about to hand to meiosis() and InheritMutations() — the same masks, so the
// record is the truth about this birth rather than a reconstruction of it.
//
// Storage is one flat allocation per birth: gamete 0 occupies words [0,w), gamete 1
// occupies [w,2w). A SET bit means the locus came from the parent's strand 1; this is
// the inverse of the mask bit, per the convention documented at the top of this file.
func (p *Pop) ARGCapture(child int, loci []int, maskDad, maskMom []uint64) {
	if len(loci) == 0 {
		return
	}
	if p.ARGDB == nil {
		p.ARGDB = make(map[int][]uint64)
	}
	w := argWords(len(loci))
	rec := make([]uint64, 2*w)
	for g, mask := range [2][]uint64{maskDad, maskMom} {
		base := g * w
		for i, pos := range loci {
			word := pos / 64
			if word >= len(mask) {
				continue // out of the mask's range: leave as strand 0
			}
			if mask[word]&(1<<(uint(pos)%64)) == 0 {
				// mask bit clear ⇒ meiosis took parental copy 1
				rec[base+i/64] |= 1 << (uint(i) % 64)
			}
		}
	}
	p.ARGDB[child] = rec
}

// ARGStrand returns which parental strand (0 or 1) individual `child` inherited at focal
// locus index `locusIdx` on gamete `gamete` (0 = paternal, 1 = maternal), and whether a
// record exists. ok=false means the walk has reached a founder — or an individual born
// before capture was enabled — and must terminate rather than guess.
func (p *Pop) ARGStrand(child, gamete, locusIdx, nLoci int) (int, bool) {
	if p == nil || p.ARGDB == nil || gamete < 0 || gamete > 1 || locusIdx < 0 || locusIdx >= nLoci {
		return 0, false
	}
	rec, ok := p.ARGDB[child]
	if !ok {
		return 0, false
	}
	w := argWords(nLoci)
	if len(rec) < 2*w {
		return 0, false
	}
	idx := gamete*w + locusIdx/64
	if rec[idx]&(1<<(uint(locusIdx)%64)) != 0 {
		return 1, true
	}
	return 0, true
}
