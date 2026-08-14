package validation

// The TMR-K-A analytic anchor — TMR4A.md work item W8, fourth bullet: "a neutral
// large-N run where the true TMR4A is known analytically-ish, as a sanity anchor."
//
// WHAT THIS IS ANCHORING. W5's walker reads the true genealogy off the recorded pedigree
// (W1) and strand provenance (W2). Nothing in the rest of W8 checks its ABSOLUTE
// timescale: the hand-built-pedigree tests pin exact coalescence years on genealogies
// small enough to verify by eye, the sweep tests pin reproducibility, and the byte-
// identity tests pin inertness. All of those would still pass if the walk were, say,
// systematically reporting half the depth. This anchor is the one that would not.
//
// HOW. Run the §6h neutral scenario — the real engine, strictly neutral, one panmictic
// population at known census N — with genealogy capture on, then measure the effective
// population size TWICE from the SAME run, by two routes that share no code:
//
//	Ne_pool  = θ_W / 2μ            from the mutation pool's segregating sites (§6h)
//	Ne_geneal = TMR-K-A / 4(1/K − 1/n)   from the recorded genealogy (W5)
//
// The first is a statement about total branch length inferred from sequence; the second
// is a statement about the depth at which the sample falls to K lineages, read directly
// off who actually descended from whom. They are computed by different packages from
// different substrates and they must land on the same population.
//
// WHAT IS ASSERTED TIGHTLY AND WHAT IS ASSERTED IN A BAND. Two properties hold for ANY
// genealogy, Kingman-shaped or not, and are asserted exactly:
//
//   - TMR-K-A is non-increasing in K at every locus. More lineages is never deeper.
//   - Ne_geneal(K) computed at different K must agree, because the 4Ne(1/K − 1/n) factor
//     is the only K-dependence a coalescent admits. Spread across K is a direct,
//     Ne-free measurement of how far DRIFT's genealogy departs from Kingman's shape.
//
// The absolute agreement Ne_geneal ≈ Ne_pool is asserted in a CHARACTERIZED band rather
// than as the textbook identity, and what the measurement actually shows is that it is
// very close to one: at the default settings the ratio is 0.95, stable at 0.95-1.05 over
// four independent base seeds. The residual is a KNOWN and directional bias, not noise —
// see the burn-in note below.
//
// WHY THE BURN-IN IS FAR LONGER THAN §6H'S, and what that exposed. Two reasons, one
// expected and one not:
//
//  1. Expected: the walk terminates at the founding generation and reports
//     `censored_at_founding` rather than inventing a root (TMR4A.md §0). For the anchor
//     to measure a coalescence depth at all the sample has to actually coalesce before
//     the run runs out of history; a run shorter than a few times 4·Ne generations
//     censors instead, which is a true answer to a different question. The censored
//     fraction is therefore reported and bounded, never quietly dropped from the average.
//
//  2. Unexpected: at §6h's 1500-year burn-in the mutation pool has NOT reached
//     mutation-drift equilibrium, and Ne_pool is biased low as a result. Measured on the
//     §6h scenario at R=6, holding the sample size fixed so only the burn-in moves:
//
//     burn-in    1500y     3000y     6000y
//     θ_W        76.7      84.6      92.3
//     Ne/N       0.320     0.352     0.385
//     Tajima's D −0.683    −0.266    −0.279
//     θ_π/θ_W    0.802     0.923     0.919
//
//     The SFS SHAPE converges by ~3000y; the SCALE is still creeping at 6000y, which is
//     about 3.4 × 4·Ne generations. So Ne_pool remains a slight underestimate here, and
//     the anchor's ratio duly approaches 1 from below as the burn-in lengthens: 0.926 at
//     3000y, 0.947 at 4500y, 0.951 at 6000y. That is the residual, and it is in the
//     direction the mechanism predicts.
//
//     This also says something about the §6h gate itself that is worth stating plainly:
//     its characterized Tajima's D ≈ −0.66 is measured before equilibrium, and most of
//     that rare-variant excess is a burn-in artefact rather than the life-history property
//     the package comment attributes it to — by 3000y the same scenario sits at ≈ −0.27
//     with θ_π/θ_W ≈ 0.92. Nothing here changes the gate, whose whole value is being a
//     fixed regression tripwire, but the number should not be read as DRIFT's equilibrium
//     neutral baseline.

import (
	"drift/pkg/analysis"
	"math"
)

// TMRKAAnchorConfig parameterizes an anchor run: a neutral scenario plus the set of K
// values to check the coalescent prediction at.
type TMRKAAnchorConfig struct {
	Neutral NeutralConfig
	Ks      []int // K values in TMR-K-A; must be ≥1 and < the sampled lineage count
}

// DefaultTMRKAAnchorConfig returns a calibrated anchor run.
//
// The deviations from DefaultNeutralConfig are all forced by what the walk needs:
//
//   - a much longer burn-in, because the sample must coalesce to K lineages BEFORE the
//     walk reaches the founders (see the censoring note above). At N=120 the implied
//     Ne is ≈0.31·N ≈ 37 and the realised generation time is ≈12 years, so the ≈4·Ne
//     generations to a common ancestor is ≈1800 years; 6000 gives a comfortable margin.
//   - a smaller sample, because TMR-K-A depends on n only through 1/K − 1/n and a modest
//     n already puts that within a few percent of its 1/K limit, while the walk cost and
//     the pedigree retained both grow with it.
//   - fewer replicates, since each is far longer; the loci within a replicate supply
//     most of the averaging (one focal locus per synthetic chromosome, and the legacy
//     recombination model resolves each chromosome independently, so they are close to
//     independent draws of the marginal genealogy).
func DefaultTMRKAAnchorConfig() TMRKAAnchorConfig {
	n := DefaultNeutralConfig()
	n.BurnInYears = 6000
	n.SampleSize = 24
	n.Replicates = 6
	n.BaseSeed = 20260814
	n.TrackPedigree = true
	// Focal loci spread UNIFORMLY over genome positions, at a stride coprime with the
	// chromosome size so successive loci fall at different within-chromosome offsets.
	// This is load-bearing and is the subject of its own test — see the locus-placement
	// note below.
	n.ARGLoci = "0-16383:997"
	return TMRKAAnchorConfig{
		Neutral: n,
		Ks:      []int{2, 4, 8},
	}
}

// WHY THE FOCAL LOCI MUST BE SPREAD UNIFORMLY OVER POSITIONS, measured rather than
// assumed. DRIFT's default meiosis (`recombination_model="legacy"`) builds each
// chromosome's mask as ONE contiguous interior segment: a coin decides whether any
// segment is taken at all, and if it is, its endpoints are drawn uniformly in the p- and
// q-arms. The probability that a given position is inherited from parental copy 0 is
// therefore ≈1/2 at the centromere and ≈0 at a chromosome end — a position at the very
// start of an arm is taken from the same parental copy almost every meiosis, so its
// lineage barely switches strands at all.
//
// That makes the marginal genealogy strongly position-dependent, and this harness
// measures how strongly. Under the default anchor scenario, the Ne implied by the
// walker's own depths comes out:
//
//	chromosome starts   Ne_geneal ≈ 17     (lineages almost never switch strand)
//	centromeres         Ne_geneal ≈ 68     (free strand switching, deepest genealogy)
//	uniform positions   Ne_geneal ≈ 45     (matches Ne_pool ≈ 44)
//
// A four-fold spread, and only the uniform placement is comparable with the mutation
// pool, because mutations are dropped at uniformly-drawn positions and so average over
// exactly this heterogeneity. TestTMRKAAnchorDepthDependsOnLocusPlacement pins the
// contrast so that "simplifying" the locus spec to one locus per chromosome cannot
// quietly break the anchor.
//
// This is a property of the ENGINE's legacy recombination model, not of the walker — the
// walker reports the genealogy it is given, and reports it correctly at all three
// placements. It matters downstream: W7 feeds DRIFT output to ARGweaver, whose model
// assumes position-homogeneous recombination, so a W7 run must not inherit this.

// TMRKAAnchorK is the anchor's finding at one value of K, pooled over every locus of
// every replicate.
type TMRKAAnchorK struct {
	K         int
	Reached   int     // locus-genealogies that fell to K lineages
	Censored  int     // locus-genealogies that hit the founders with more than K
	MeanGens  float64 // mean TMR-K-A depth, in generations, over the Reached loci
	SDGens    float64
	ImpliedNe float64 // MeanGens / 4(1/K − 1/n): the Ne this depth implies
}

// TMRKAAnchorResult is a whole anchor run.
type TMRKAAnchorResult struct {
	Config  TMRKAAnchorConfig
	Neutral *NeutralRunResult

	Lineages       int     // n: sampled lineages = 2 × sampled individuals
	GenerationTime float64 // realised mean parent-to-child birth-year gap, in years
	CensusN        int
	NeFromPool     float64 // θ_W/2μ, the §6h implied Ne
	MeanNeFromGen  float64 // mean of ImpliedNe over the K values that had data
	CVNeFromGen    float64 // spread of ImpliedNe across K: the Ne-free shape check
	NeRatio        float64 // NeFromPool / MeanNeFromGen

	TotalLoci       int // locus-genealogies walked, over all replicates
	OrderViolations int // loci where TMR-K-A rose with K — structurally impossible
	MissingARG      int // loci the walker had no provenance for

	ByK []TMRKAAnchorK
}

// RunTMRKAAnchor runs the neutral scenario with genealogy capture on and compares the
// walker's depths against the constant-Ne coalescent prediction.
func RunTMRKAAnchor(cfg TMRKAAnchorConfig) *TMRKAAnchorResult {
	ncfg := cfg.Neutral
	ncfg.TrackPedigree = true // the anchor is meaningless without it
	neutral := RunNeutral(ncfg)

	res := &TMRKAAnchorResult{
		Config:  cfg,
		Neutral: neutral,
		CensusN: ncfg.N,
	}
	if ncfg.Mu > 0 {
		res.NeFromPool = neutral.MeanThetaW / (2 * ncfg.Mu)
	}

	// depths[K] accumulates one entry per locus-genealogy that reached K, in generations.
	depths := map[int][]float64{}
	censored := map[int]int{}
	genSum, genWeight := 0.0, 0

	for _, rr := range neutral.Replicates {
		if rr.Pop == nil || len(rr.SampleIDs) == 0 {
			continue
		}
		// Realised generation time from this replicate's own pedigree. Weighted by the
		// number of parent-child edges it averaged over, so a replicate that happened to
		// leave a smaller pedigree does not get equal say.
		gt, edges := analysis.PedigreeGenerationTime(rr.Pop)
		if edges == 0 || gt <= 0 {
			continue
		}
		genSum += gt * float64(edges)
		genWeight += edges

		walk := analysis.ComputeTMRKA(rr.Model, rr.Pop, rr.SampleIDs)
		if walk == nil {
			continue
		}
		res.MissingARG += walk.MissingARG
		if res.Lineages == 0 {
			res.Lineages = 2 * walk.NumSampled
		}

		for i := range walk.Loci {
			loc := &walk.Loci[i]
			res.TotalLoci++
			prev := math.Inf(1)
			for _, k := range cfg.Ks {
				gens, ok := loc.TMRKAGensAgo(k, gt)
				if !ok {
					censored[k]++
					continue
				}
				depths[k] = append(depths[k], gens)
				// Structural: raising K can only make the target shallower. A rise here
				// is not a tolerance question, it is a broken trajectory.
				if gens > prev {
					res.OrderViolations++
				}
				prev = gens
			}
		}
	}

	if genWeight > 0 {
		res.GenerationTime = genSum / float64(genWeight)
	}

	neSum, neSumSq, neCount := 0.0, 0.0, 0
	for _, k := range cfg.Ks {
		row := TMRKAAnchorK{K: k, Reached: len(depths[k]), Censored: censored[k]}
		if row.Reached > 0 {
			sum, sumSq := 0.0, 0.0
			for _, d := range depths[k] {
				sum += d
				sumSq += d * d
			}
			row.MeanGens = sum / float64(row.Reached)
			row.SDGens = stddev(sum, sumSq, float64(row.Reached))
			if ne, ok := analysis.ImpliedCoalescentNe(row.MeanGens, k, res.Lineages); ok {
				row.ImpliedNe = ne
				neSum += ne
				neSumSq += ne * ne
				neCount++
			}
		}
		res.ByK = append(res.ByK, row)
	}
	if neCount > 0 {
		res.MeanNeFromGen = neSum / float64(neCount)
		if res.MeanNeFromGen > 0 {
			res.CVNeFromGen = stddev(neSum, neSumSq, float64(neCount)) / res.MeanNeFromGen
			res.NeRatio = res.NeFromPool / res.MeanNeFromGen
		}
	}
	return res
}

// CensoredFraction is the share of locus-genealogies that hit the founding generation
// still holding more than K lineages, at the SMALLEST K in the run — the hardest target
// and therefore the binding one. A high value means the run was too short for the sample
// to coalesce and the anchor is measuring the founding date rather than a coalescence
// depth; it is reported rather than silently absorbed into the averages, per TMR4A.md §6
// risk 5.
func (r *TMRKAAnchorResult) CensoredFraction() float64 {
	if r.TotalLoci == 0 || len(r.ByK) == 0 {
		return 0
	}
	worst := r.ByK[0]
	for _, row := range r.ByK {
		if row.K < worst.K {
			worst = row
		}
	}
	return float64(worst.Censored) / float64(r.TotalLoci)
}
