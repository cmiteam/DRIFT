package analysis

// The analytic anchor for the TMR-K-A walker — TMR4A.md work item W8, fourth bullet:
// "a neutral large-N run where the true TMR4A is known analytically-ish, as a sanity
// anchor."
//
// WHAT THERE IS TO ANCHOR AGAINST. Under Kingman's coalescent with diploid effective
// size Ne, i ancestral lineages coalesce at rate C(i,2)/(2Ne) per generation, so the
// expected wait with i lineages is 4Ne/(i(i−1)) and the expected time for n sampled
// lineages to fall to k is the telescoping sum
//
//	E[T(n→k)] = Σ_{i=k+1}^{n} 4Ne/(i(i−1)) = 4Ne(1/k − 1/n) generations.
//
// That is the same formula TMR4A.md §1B uses to make its central point about the
// published analysis — at k=4 and n=108 it gives ≈0.96·Ne generations, so an inferred
// TMR4A tracks its own assumed Ne before any data is consulted. Here it is turned around
// and used as a test: DRIFT knows the true genealogy, so the walker's TMR-K-A can be
// converted back into an implied Ne and compared against an Ne measured a completely
// different way (θ_W/2μ from the mutation pool). Two independent routes to the same
// number is the anchor.
//
// WHY "ANALYTICALLY-ISH" AND NOT "ANALYTICALLY". DRIFT is not Wright-Fisher. It is
// monogamous, overlapping-generations, age-structured, and birth-then-cull, and §6h
// records what that does: a reproducible excess of rare variants (Tajima's D ≈ −0.66,
// θ_π/θ_W ≈ 0.81) which is exactly a statement that its genealogies are NOT Kingman-
// shaped. So the formula is the right shape and the right scale, not an identity, and
// the checks built on it belong in a characterized band — the same discipline
// CharacterizedTolerances applies to Tajima's D. What survives the non-Kingman skew, and
// is asserted tightly, is the ORDERING: TMR-K-A must fall as K rises, at every locus, for
// any genealogy whatever.
//
// UNITS. The walker reports simulated YEARS because that is what DRIFT's overlapping
// generations make meaningful. Coalescent theory is in generations. The bridge is
// PedigreeGenerationTime, which measures the realised mean parent-to-child birth-year gap
// from the pedigree itself rather than trusting the `generation_time` parameter, which is
// a reporting convenience and need not match what the life history actually did.

import "drift/pkg/core"

// PedigreeGenerationTime returns the realised generation time in simulated years: the
// mean parent-to-child birth-year gap over every recorded parent-child edge in the
// pedigree, together with the number of edges it averaged over.
//
// Both parents of every child contribute an edge, so this is the mean age at reproduction
// weighted by offspring — which is the generation time a coalescent rescaling wants, not
// the mean age of parents. Founders contribute nothing (they have no recorded parents and
// their birth years are pre-simulation). ok is false when the pedigree holds no edge at
// all, in which case there is nothing to measure and a caller must not silently divide by
// a made-up number.
func PedigreeGenerationTime(pop *core.Pop) (float64, int) {
	if pop == nil || pop.Pedigree == nil {
		return 0, 0
	}
	sum, edges := 0, 0
	for _, n := range pop.Pedigree {
		for _, parent := range [2]int{n.Dad, n.Mom} {
			if parent == core.PedFounder {
				continue
			}
			pn, ok := pop.Pedigree[parent]
			if !ok {
				continue
			}
			gap := n.BirthYear - pn.BirthYear
			if gap <= 0 {
				continue // a parent must predate its child; anything else is not a generation
			}
			sum += gap
			edges++
		}
	}
	if edges == 0 {
		return 0, 0
	}
	return float64(sum) / float64(edges), edges
}

// coalescentFactor is (1/k − 1/n), the dimensionless part of E[T(n→k)] = 4Ne(1/k − 1/n).
// ok is false when the quantity is not defined: fewer than two lineages, or a k that the
// sample is already at or below, where no waiting time exists to predict.
func coalescentFactor(k, n int) (float64, bool) {
	if k < 1 || n < 2 || k >= n {
		return 0, false
	}
	return 1/float64(k) - 1/float64(n), true
}

// ExpectedCoalescentGens returns E[T(n→k)] = 4Ne(1/k − 1/n) generations: the expected
// time for n sampled lineages to fall to k under a constant-Ne Kingman coalescent.
// n is a count of LINEAGES (2 × sampled diploids), not of individuals.
func ExpectedCoalescentGens(ne float64, k, n int) (float64, bool) {
	f, ok := coalescentFactor(k, n)
	if !ok || ne <= 0 {
		return 0, false
	}
	return 4 * ne * f, true
}

// ImpliedCoalescentNe inverts ExpectedCoalescentGens: given an OBSERVED TMR-K-A of `gens`
// generations from n sampled lineages, the Ne a constant-Ne coalescent would need to
// produce it. This is precisely the inference TMR4A.md §1B says is running backwards in
// the published analysis — there Ne is assumed and the depth follows; here the depth is
// known truth and the Ne it implies is the thing being checked.
func ImpliedCoalescentNe(gens float64, k, n int) (float64, bool) {
	f, ok := coalescentFactor(k, n)
	if !ok || gens <= 0 {
		return 0, false
	}
	return gens / (4 * f), true
}

// TMRKAGensAgo converts a locus's TMR-K-A for the given K into generations before the
// sample year, using a realised generation time in years. ok is false when the walk never
// reached K lineages — a censored locus, which has no depth to report and must never be
// averaged into one (TMR4A.md §6 risk 5).
func (l *TMRKALocus) TMRKAGensAgo(k int, genTimeYears float64) (float64, bool) {
	if genTimeYears <= 0 {
		return 0, false
	}
	year, ok := l.TMRKAForK(k)
	if !ok {
		return 0, false
	}
	return float64(l.SampleYear-year) / genTimeYears, true
}
