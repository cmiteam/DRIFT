package analysis

// The realised generation time, measured rather than declared — the companion to
// argweaver_scaling.go's two rate anchors, and the third number every year-denominated
// claim in this study silently depends on.
//
// WHY THIS IS NOT THE `generation_time` PARAMETER. Nothing in the engine reads
// `generation_time`. Individuals mature at `maturity`, conceive with probability
// 1/`birth_prob` per year subject to `spacing`, and stop at `menopause` × their own
// lifespan — and that lifespan is itself inherited and decayed by `lifespan_drop` every
// generation. The years per generation that results is an EMERGENT property of that life
// history. `generation_time` is a reporting convenience that no simulated individual
// consults, so it can be — and on the patriarchal models is — wrong by a large factor.
//
// WHY THE MEAN PARENT-CHILD GAP IS THE RIGHT CONSTANT, AND NOT THE MEDIAN. ARGweaver's
// model is in GENERATIONS: it is told mutations and recombinations per site per
// generation and returns depths in generations. A "generation" in that accounting is one
// parent→child transmission, because GenerateNewMutations fires exactly once per child
// (pkg/simulation/birth.go). So the conversion constant is the mean years per
// transmission event. A lineage's elapsed time is a SUM of such gaps, so by the CLT the
// MEAN is what accumulates, whatever shape the per-gap distribution has. On a
// right-skewed distribution the median is materially smaller, and converting with it
// under-dates everything.
//
// BOTH PARENTS COUNT, AND THE PEDIGREE IS ALREADY THE RIGHT SAMPLE. Each child
// contributes two edges, so this is the mean age at reproduction weighted by offspring —
// the quantity a coalescent rescaling wants, not the mean age of parents. And the
// pedigree is refcount-pruned to ancestors of the living (pkg/core/pedigree.go), so the
// edges measured here are the ones an ARG actually traverses rather than the full birth
// register.
//
// WHAT A SINGLE CONSTANT CANNOT SAY. On a decaying-lifespan model the generation time
// FALLS over the run, so years and generations are related by a cumulative map, not a
// multiplier. ARGweaver has no way to be told this — its time grid is in generations and
// its likelihood only ever sees the per-generation rates, which ARE stationary. That is
// the saving grace: the non-stationarity lives entirely in the OUTPUT conversion, so it
// is fixable after the fact. Bins below is that map, and Drift is the diagnostic that
// says whether ignoring it matters for this model.

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
)

// GenerationBin is one generation since founding: the local years-per-generation at that
// depth and the calendar year the generation sits at. Together the bins are the
// generations <-> years dictionary that replaces a single constant when the life history
// is not stationary.
type GenerationBin struct {
	Gen           int     // 1 = children of founders
	Edges         int     // parent-child edges whose CHILD sits at this depth
	MeanGap       float64 // years per generation, locally
	SDGap         float64
	MeanBirthYear float64 // mean birth year of the children at this depth
}

// GenTimeStats is the realised generation time of one run, with the spread and the
// non-stationarity profile that the point estimate hides.
type GenTimeStats struct {
	// The point estimate and its spread. SEM is the standard error of the MEAN and is
	// the error bar on the conversion constant; SD is the spread of individual gaps and
	// is a description of the life history, NEVER an error bar on a date.
	Edges                               int
	Mean, SD, CV, Skew, SEM             float64
	Min, P5, Q25, Median, Q75, P95, Max float64

	Paternal, Maternal           float64
	PaternalEdges, MaternalEdges int

	// Non-stationarity, measured by GENERATION rather than by calendar year — a
	// founding couple that reproduces for centuries is three pedigree generations and
	// almost no edges, so a year-split averages it away and reports a stationary life
	// history that is not one.
	//
	// Plateau is the settled rate the run converged to and PlateauFrom the generation it
	// settles at; PrePlateau is the founder transient before it. Drift is their ratio —
	// the one number that says whether a single constant is defensible for this model.
	//
	// Plateau, not Mean, is the constant to hand a deep inference: arg-sample returns
	// depths far past the real history, no measured generation time exists out there,
	// and the settled rate is the only defensible characterisation of a generation in
	// this population. Mean folds in a founder transient that no inferred deep
	// generation corresponds to.
	Bins                       []GenerationBin
	Plateau, PrePlateau, Drift float64
	PlateauFrom                int
	PlateauEdges               int
	PrePlateauEdges            int
	// CensoredFrom is the first generation whose children were born less than one
	// generation before the sample, so their own children are not in the record yet and
	// their apparent gap is an artifact of the run ending. 0 when none are.
	CensoredFrom int

	// Generational depth of the living population: the truth depth in GENERATIONS,
	// which needs no conversion constant at all and is therefore the unit a method-error
	// ratio should be quoted in.
	DepthN             int
	MeanDepth, SDDepth float64
	MinDepth, MaxDepth float64

	Declared   float64 // the generation_time parameter, for comparison only
	SampleYear int
	FirstYear  int
}

// PedigreeGenerationTime returns the realised generation time in simulated years and the
// number of parent-child edges it averaged over. Preserved as the narrow entry point the
// TMR-K-A anchor uses; PedigreeGenerationStats is the same walk with the spread kept.
func PedigreeGenerationTime(pop *core.Pop) (float64, int) {
	s := PedigreeGenerationStats(pop, nil)
	if s == nil {
		return 0, 0
	}
	return s.Mean, s.Edges
}

// PedigreeGenerationStats measures the generation time and its distribution from the
// pedigree. model may be nil; it supplies only the declared parameter and the sample year
// used for reporting. Returns nil when there is no pedigree edge to measure, because an
// invented constant is worse than an absent one.
func PedigreeGenerationStats(pop *core.Pop, model *core.Model) *GenTimeStats {
	if pop == nil || pop.Pedigree == nil {
		return nil
	}
	depth := pedigreeDepths(pop)

	var gaps, pat, mat []float64
	type acc struct {
		gaps  []float64
		years []float64
	}
	byGen := map[int]*acc{}
	minYear, maxYear := math.Inf(1), math.Inf(-1)

	for id, n := range pop.Pedigree {
		for slot, parent := range [2]int{n.Dad, n.Mom} {
			if parent == core.PedFounder {
				continue
			}
			pn, ok := pop.Pedigree[parent]
			if !ok {
				continue
			}
			gap := float64(n.BirthYear - pn.BirthYear)
			if gap <= 0 {
				continue // a parent must predate its child; anything else is not a generation
			}
			gaps = append(gaps, gap)
			if slot == 0 {
				pat = append(pat, gap)
			} else {
				mat = append(mat, gap)
			}
			y := float64(n.BirthYear)
			minYear = math.Min(minYear, y)
			maxYear = math.Max(maxYear, y)

			g := int(math.Round(depth[id]))
			a := byGen[g]
			if a == nil {
				a = &acc{}
				byGen[g] = a
			}
			a.gaps = append(a.gaps, gap)
			a.years = append(a.years, y)
		}
	}
	if len(gaps) == 0 {
		return nil
	}

	s := &GenTimeStats{Edges: len(gaps)}
	s.Mean, s.SD, s.Skew = momentStats(gaps)
	if s.Mean > 0 {
		s.CV = s.SD / s.Mean
	}
	s.SEM = s.SD / math.Sqrt(float64(len(gaps)))
	s.Min, s.P5, s.Q25, s.Median, s.Q75, s.P95, s.Max = gapQuantiles(gaps)
	s.Paternal, _, _ = momentStats(pat)
	s.Maternal, _, _ = momentStats(mat)
	s.PaternalEdges, s.MaternalEdges = len(pat), len(mat)

	// The per-generation map, ordered shallowest (nearest the founders) first.
	var gens []int
	for g := range byGen {
		gens = append(gens, g)
	}
	sort.Ints(gens)
	for _, g := range gens {
		a := byGen[g]
		m, sd, _ := momentStats(a.gaps)
		my, _, _ := momentStats(a.years)
		s.Bins = append(s.Bins, GenerationBin{
			Gen: g, Edges: len(a.gaps), MeanGap: m, SDGap: sd, MeanBirthYear: my,
		})
	}

	s.FirstYear = int(minYear)
	if model != nil {
		s.Declared = model.Parameters["generation_time"]
		s.SampleYear = model.FreeParameters["year"]
	}
	if s.SampleYear <= 0 {
		// No model to ask: the sample can be no earlier than the last recorded birth,
		// and the right-censoring test needs SOME present to measure back from.
		s.SampleYear = int(maxYear)
	}
	s.findPlateau()

	// Generational depth of the living population.
	var ds []float64
	for id := range pop.IndData {
		if _, ok := pop.Pedigree[id]; ok {
			ds = append(ds, depth[id])
		}
	}
	if len(ds) > 0 {
		s.DepthN = len(ds)
		s.MeanDepth, s.SDDepth, _ = momentStats(ds)
		sort.Float64s(ds)
		s.MinDepth, s.MaxDepth = ds[0], ds[len(ds)-1]
	}

	return s
}

// plateauTolerance is how far a generation's rate may sit from the settled mean and
// still count as settled. 20% is wide enough to absorb the sampling noise of a few
// hundred edges per generation and narrow enough that a founder transient — which on the
// shipped patriarchal models runs 2-4x the settled rate — cannot hide inside it.
const plateauTolerance = 0.20

// findPlateau splits the per-generation map into a founder transient and the settled rate
// the run converged to, and measures how far apart they are.
//
// TWO THINGS HAVE TO BE HANDLED, and they sit at opposite ends of the map.
//
// RIGHT CENSORING at the present end. The youngest generations' children were born too
// recently to have had their own children recorded, so only their EARLIEST births are in
// the pedigree and their apparent gap collapses towards `maturity`. That is an artifact
// of the run ending, not a life history that sped up, and averaging it in would drag the
// settled rate down. A generation is censored when its children were born less than one
// mean generation before the sample.
//
// THE FOUNDER TRANSIENT at the other end. On a decaying-lifespan model the founders
// reproduce for centuries, so generations 1-3 can run several times the settled rate.
// They carry very few edges, which is exactly why an edge-weighted or calendar-year
// diagnostic reports a stationary life history that is not one — and why the split has to
// be made over GENERATIONS.
func (s *GenTimeStats) findPlateau() {
	if len(s.Bins) == 0 {
		return
	}
	// Trim the right-censored tail.
	usable := len(s.Bins)
	for usable > 0 {
		b := s.Bins[usable-1]
		if float64(s.SampleYear)-b.MeanBirthYear >= s.Mean {
			break
		}
		usable--
	}
	if usable < len(s.Bins) {
		s.CensoredFrom = s.Bins[usable].Gen
	}
	if usable == 0 {
		// Everything is censored: there is no settled rate to report, so report none
		// rather than one computed from the artifact.
		return
	}
	bins := s.Bins[:usable]

	// The plateau is the longest SUFFIX of the map whose generations all sit within
	// tolerance of that suffix's own edge-weighted mean. Scanning from generation 1
	// forward takes the earliest such suffix, so only generations that genuinely have
	// not settled are called a transient.
	for start := 0; start < len(bins); start++ {
		mean := weightedGap(bins[start:])
		if mean <= 0 {
			continue
		}
		settled := true
		for _, b := range bins[start:] {
			if math.Abs(b.MeanGap-mean)/mean > plateauTolerance {
				settled = false
				break
			}
		}
		if !settled {
			continue
		}
		s.Plateau, s.PlateauFrom = mean, bins[start].Gen
		for _, b := range bins[start:] {
			s.PlateauEdges += b.Edges
		}
		if start > 0 {
			s.PrePlateau = weightedGap(bins[:start])
			for _, b := range bins[:start] {
				s.PrePlateauEdges += b.Edges
			}
			if s.Plateau > 0 {
				s.Drift = s.PrePlateau / s.Plateau
			}
		}
		return
	}
	// Never settles within tolerance: the whole map is a transient. Say so by leaving
	// Plateau at the deepest usable generation alone, which is the best available guess.
	last := bins[len(bins)-1]
	s.Plateau, s.PlateauFrom, s.PlateauEdges = last.MeanGap, last.Gen, last.Edges
	if len(bins) > 1 {
		s.PrePlateau = weightedGap(bins[:len(bins)-1])
		for _, b := range bins[:len(bins)-1] {
			s.PrePlateauEdges += b.Edges
		}
		if s.Plateau > 0 {
			s.Drift = s.PrePlateau / s.Plateau
		}
	}
}

// weightedGap is the edge-weighted mean gap over a run of generations — the mean an
// individual transmission sees, not the mean of the per-generation means.
func weightedGap(bins []GenerationBin) float64 {
	var sum float64
	var n int
	for _, b := range bins {
		sum += b.MeanGap * float64(b.Edges)
		n += b.Edges
	}
	if n == 0 {
		return 0
	}
	return sum / float64(n)
}

// Recommended is the single constant to hand ARGweaver's output when only one number can
// be used: the settled rate where the run reached one, else the overall mean. On a
// stationary life history the two agree, which is the whole reason a single constant is
// ever defensible.
func (s *GenTimeStats) Recommended() float64 {
	if s == nil {
		return 0
	}
	if s.Plateau > 0 {
		return s.Plateau
	}
	return s.Mean
}

// ConversionConstant is the years-per-generation to apply to an inferred depth of g
// generations. Inside the real history the per-generation map is exact; past it — which
// is where an over-deep ARG inference always lands — there is no measured value and the
// settled rate is the defensible extrapolation. deep reports which case applied.
func (s *GenTimeStats) ConversionConstant(g float64) (rate float64, deep bool) {
	if s == nil {
		return 0, false
	}
	if len(s.Bins) > 0 {
		maxG := float64(s.Bins[len(s.Bins)-1].Gen - s.Bins[0].Gen)
		if g <= maxG {
			if years, _ := s.YearsForGens(g); g > 0 && years > 0 {
				return years / g, false
			}
		}
	}
	if s.Plateau > 0 {
		return s.Plateau, true
	}
	return s.Mean, true
}

// ResolveGenerationTime returns the years-per-generation that this run's year <->
// generation conversions should use, together with the measured statistics (nil when
// there is no pedigree to measure).
//
// DEFAULT IS THE DECLARED PARAMETER, deliberately: switching the constant silently would
// rewrite the `gens_ago` column of every CSV this package emits and break comparison
// against already-archived runs. Set `generation_time_source=pedigree` to convert with
// the value the life history actually realised. Either way the measured number and the
// disagreement are reported — the default is quiet, not hidden.
func ResolveGenerationTime(model *core.Model, pop *core.Pop) (float64, *GenTimeStats) {
	stats := PedigreeGenerationStats(pop, model)

	declared := 0.0
	if model != nil {
		declared = model.Parameters["generation_time"]
	}
	if declared <= 0 {
		declared = 25
	}
	if model != nil && model.StringParam("generation_time_source", "declared") == "pedigree" {
		if stats != nil && stats.Mean > 0 {
			return stats.Mean, stats
		}
		fmt.Println("WARNING: generation_time_source=pedigree but no pedigree edges were " +
			"recorded (is track_pedigree on?); falling back to the declared generation_time")
	}
	return declared, stats
}

// pedigreeDepths returns each pedigree member's generational distance from founding: 0
// for a founder, otherwise one more than the MEAN depth of its recorded parents. The mean
// (not the min) is used because pedigree collapse makes the two disagree, and a lineage's
// expected transmission count is what the years <-> generations map needs.
func pedigreeDepths(pop *core.Pop) map[int]float64 {
	memo := make(map[int]float64, len(pop.Pedigree))
	done := make(map[int]bool, len(pop.Pedigree))
	var walk func(id int) float64
	walk = func(id int) float64 {
		if done[id] {
			return memo[id]
		}
		n, ok := pop.Pedigree[id]
		if !ok {
			return 0
		}
		done[id] = true // cycle guard: a pedigree loop resolves at the founder depth
		memo[id] = 0
		var sum float64
		var cnt int
		for _, p := range [2]int{n.Dad, n.Mom} {
			if p == core.PedFounder {
				continue
			}
			if _, ok := pop.Pedigree[p]; !ok {
				continue
			}
			sum += walk(p)
			cnt++
		}
		if cnt == 0 {
			if n.Dad == core.PedFounder && n.Mom == core.PedFounder {
				memo[id] = 0 // a founder
			} else {
				memo[id] = 1 // parents recorded but pruned away
			}
			return memo[id]
		}
		memo[id] = 1 + sum/float64(cnt)
		return memo[id]
	}
	for id := range pop.Pedigree {
		walk(id)
	}
	return memo
}

// YearsForGens converts a depth in generations before the sample into simulated years,
// through the measured per-generation map rather than by a constant. extrapolated is true
// when the depth reaches past the pedigree, where no generation time was ever measured
// and the late-run rate is assumed — which is every inferred depth that exceeds the real
// history, and is exactly the regime an over-deep ARG inference lands in.
func (s *GenTimeStats) YearsForGens(g float64) (years float64, extrapolated bool) {
	if s == nil || g <= 0 {
		return 0, false
	}
	if len(s.Bins) == 0 {
		return g * s.Mean, true
	}
	deepest := s.Bins[len(s.Bins)-1] // nearest the present
	shallowest := s.Bins[0]          // nearest the founders
	maxG := float64(deepest.Gen - shallowest.Gen)
	if g <= maxG {
		target := float64(deepest.Gen) - g
		return float64(s.SampleYear) - interpBinYear(s.Bins, target), false
	}
	rate := s.Recommended()
	base := float64(s.SampleYear) - shallowest.MeanBirthYear
	return base + (g-maxG)*rate, true
}

// interpBinYear linearly interpolates the mean birth year between adjacent bins.
func interpBinYear(bins []GenerationBin, gen float64) float64 {
	if gen <= float64(bins[0].Gen) {
		return bins[0].MeanBirthYear
	}
	last := bins[len(bins)-1]
	if gen >= float64(last.Gen) {
		return last.MeanBirthYear
	}
	for i := 1; i < len(bins); i++ {
		lo, hi := bins[i-1], bins[i]
		if gen <= float64(hi.Gen) {
			span := float64(hi.Gen - lo.Gen)
			if span <= 0 {
				return hi.MeanBirthYear
			}
			f := (gen - float64(lo.Gen)) / span
			return lo.MeanBirthYear + f*(hi.MeanBirthYear-lo.MeanBirthYear)
		}
	}
	return last.MeanBirthYear
}

// Warnings reports every way the declared constant misrepresents what the life history
// actually did. Same contract as RateAnchors.check: say it in the file that ships with
// the data, not in a footnote somebody has to find.
func (s *GenTimeStats) Warnings() []string {
	if s == nil {
		return []string{"generation time could NOT be measured: no pedigree edges. Set " +
			"track_pedigree=1 to verify the generation_time parameter instead of trusting it — " +
			"every year-denominated depth from this export rests on it."}
	}
	var w []string
	rec := s.Recommended()
	if s.Declared > 0 && rec > 0 {
		f := rec / s.Declared
		if f < 0.9 || f > 1.1 {
			w = append(w, fmt.Sprintf("generation_time is declared %.4g y but the pedigree realised "+
				"%.4g y (%.3gx). ARGweaver returns GENERATIONS, so every year-denominated depth built "+
				"on the declared value is wrong by that same %.3gx. Pass -argweaver-gen-time %.4g "+
				"when parsing this run's output.", s.Declared, rec, f, f, rec))
		}
	}
	if s.Drift > 1.15 || (s.Drift > 0 && s.Drift < 0.87) {
		w = append(w, fmt.Sprintf("generation time is NOT stationary: generations 1-%d run %.4g y "+
			"against %.4g y from generation %d on (%.3gx), so years and generations are related by a "+
			"MAP, not a constant. ARGweaver cannot be told this. Those generations are only %d of %d "+
			"edges but they are the DEEPEST part of every genealogy, so a flat constant mistimes "+
			"exactly the oldest coalescences. Convert shallow (truth-side) depths through the "+
			"companion _generation_time.csv, use the settled %.4g y for depths past the real history, "+
			"and prefer quoting method error in GENERATIONS, where no constant is needed at all.",
			s.PlateauFrom-1, s.PrePlateau, s.Plateau, s.PlateauFrom, s.Drift,
			s.PrePlateauEdges, s.Edges, s.Plateau))
	}
	if s.CensoredFrom > 0 {
		w = append(w, fmt.Sprintf("generations %d and later are right-censored: their children were "+
			"born less than one generation before the sample, so only the earliest of their own "+
			"births are on record and their apparent gap is short. They are excluded from the "+
			"settled rate and must not be read as the life history speeding up.", s.CensoredFrom))
	}
	if math.Abs(s.Skew) > 1 && s.Median > 0 {
		w = append(w, fmt.Sprintf("the age-at-reproduction distribution is skewed (skew %+.3g): mean "+
			"%.4g y against median %.4g y, %.1f%% apart. The MEAN is the conversion constant — a "+
			"lineage's elapsed time is a sum of gaps — so do not quote the median as a 'typical' "+
			"generation and convert with it.", s.Skew, s.Mean, s.Median,
			100*(s.Mean-s.Median)/s.Median))
	}
	if s.Edges < 100 {
		w = append(w, fmt.Sprintf("only %d parent-child edges: the generation time is poorly "+
			"determined (SEM %.4g y). Lengthen the run or raise the population before relying on it.",
			s.Edges, s.SEM))
	}
	return w
}

// WriteReport renders the generation-time section of the scaling rationale.
func (s *GenTimeStats) WriteReport(w io.Writer) {
	fmt.Fprintln(w, "THE REALISED GENERATION TIME")
	if s == nil {
		fmt.Fprintln(w, "  NOT MEASURED — no pedigree. Set track_pedigree=1; until then the")
		fmt.Fprintln(w, "  generation_time parameter is an assertion, not a measurement.")
		fmt.Fprintln(w)
		return
	}
	fmt.Fprintln(w, "  arg-sample reports depths in GENERATIONS. One generation is one parent->child")
	fmt.Fprintln(w, "  transmission — the same event that carries one round of mutations — so the")
	fmt.Fprintln(w, "  constant that converts its output to years is the mean parent-child birth-year")
	fmt.Fprintln(w, "  gap, measured here from the pedigree rather than read off generation_time.")
	fmt.Fprintf(w, "    mean over the whole pedigree        %12.6g y  +/- %.4g (SEM over %d edges)\n",
		s.Mean, s.SEM, s.Edges)
	fmt.Fprintf(w, "    SD of individual gaps               %12.6g y  (CV %.3f, skew %+.3g)\n",
		s.SD, s.CV, s.Skew)
	fmt.Fprintf(w, "    median / q25 / q75                  %12.6g / %.4g / %.4g y\n",
		s.Median, s.Q25, s.Q75)
	fmt.Fprintf(w, "    min / p5 / p95 / max                %12.6g / %.4g / %.4g / %.4g y\n",
		s.Min, s.P5, s.P95, s.Max)
	fmt.Fprintf(w, "    paternal / maternal mean            %12.6g / %.4g y\n", s.Paternal, s.Maternal)
	if s.Plateau > 0 {
		fmt.Fprintf(w, "    settled rate, generation %2d on      %12.6g y  (%d edges)\n",
			s.PlateauFrom, s.Plateau, s.PlateauEdges)
	}
	if s.Declared > 0 {
		fmt.Fprintf(w, "    generation_time parameter           %12.6g y  (realised is %.3gx this)\n",
			s.Declared, s.Recommended()/s.Declared)
	}
	fmt.Fprintf(w, "  ==> PASS -argweaver-gen-time %.6g WHEN PARSING THIS RUN'S OUTPUT.\n", s.Recommended())
	fmt.Fprintln(w, "  The SD is the spread of the life history, NOT the error bar on a date: a depth")
	fmt.Fprintln(w, "  of g generations is a SUM of g gaps, so the conversion's own uncertainty falls")
	fmt.Fprintln(w, "  as SD/sqrt(g). Use the SEM above.")
	fmt.Fprintln(w)

	fmt.Fprintln(w, "  IS IT STATIONARY? ARGweaver assumes one generation time for all of history, and")
	fmt.Fprintln(w, "  cannot be told otherwise — its time grid is in generations and its likelihood")
	fmt.Fprintln(w, "  only ever sees the per-generation rates, which ARE constant here. So a moving")
	fmt.Fprintln(w, "  generation time is not a modelling error in the run; it is a conversion problem")
	fmt.Fprintln(w, "  in the OUTPUT, and it is fixable after the fact. Measured by generation:")
	if s.PrePlateau > 0 {
		fmt.Fprintf(w, "    generations 1-%-2d (founder transient)%12.6g y  (%d edges)\n",
			s.PlateauFrom-1, s.PrePlateau, s.PrePlateauEdges)
		fmt.Fprintf(w, "    generation %2d onward (settled)      %12.6g y  (%d edges)\n",
			s.PlateauFrom, s.Plateau, s.PlateauEdges)
		fmt.Fprintf(w, "    ratio                               %12.3f\n", s.Drift)
		fmt.Fprintln(w, "  The transient is few edges but it is the DEEPEST part of every genealogy, so a")
		fmt.Fprintln(w, "  flat constant mistimes exactly the oldest coalescences. Convert truth-side")
		fmt.Fprintln(w, "  depths through the map below; use the settled rate past the real history.")
	} else {
		fmt.Fprintf(w, "    stationary from generation %d on — a single constant is defensible here.\n",
			s.PlateauFrom)
	}
	if s.CensoredFrom > 0 {
		fmt.Fprintf(w, "    generations %d+ are right-censored (children born under one generation before\n",
			s.CensoredFrom)
		fmt.Fprintln(w, "    the sample, so their own births are not all on record) and are excluded.")
	}
	if s.DepthN > 0 {
		fmt.Fprintf(w, "  Living individuals sit %.3f +/- %.3f generations from founding (range %.0f-%.0f).\n",
			s.MeanDepth, s.SDDepth, s.MinDepth, s.MaxDepth)
		fmt.Fprintln(w, "  That is the TRUE depth in GENERATIONS, with no conversion constant involved.")
		fmt.Fprintln(w, "  Quote a method-error ratio against it and this whole question goes away.")
	}
	fmt.Fprintln(w, "  The per-generation map is in the companion _generation_time.csv.")
	fmt.Fprintln(w)
}

// WriteGenTimeCSV writes the per-generation map: the local years-per-generation at each
// depth and the calendar year it sits at. This is what a converter needs when one
// constant will not do.
func WriteGenTimeCSV(path string, s *GenTimeStats) error {
	if s == nil {
		return nil
	}
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create generation-time map: %w", err)
	}
	defer f.Close()
	cw := csv.NewWriter(f)
	defer cw.Flush()

	if err := cw.Write([]string{"generation_from_founding", "edges", "mean_gap_yr",
		"sd_gap_yr", "mean_birth_year", "years_before_sample",
		"generations_before_sample"}); err != nil {
		return err
	}
	deepest := 0
	if len(s.Bins) > 0 {
		deepest = s.Bins[len(s.Bins)-1].Gen
	}
	for _, b := range s.Bins {
		if err := cw.Write([]string{
			fmt.Sprintf("%d", b.Gen),
			fmt.Sprintf("%d", b.Edges),
			fmt.Sprintf("%.6g", b.MeanGap),
			fmt.Sprintf("%.6g", b.SDGap),
			fmt.Sprintf("%.6g", b.MeanBirthYear),
			fmt.Sprintf("%.6g", float64(s.SampleYear)-b.MeanBirthYear),
			fmt.Sprintf("%d", deepest-b.Gen),
		}); err != nil {
			return err
		}
	}
	return cw.Error()
}

// PrintGenTime renders the one-screen console summary.
func PrintGenTime(s *GenTimeStats) {
	if s == nil {
		fmt.Println("Generation time: NOT MEASURED (no pedigree; set track_pedigree=1)")
		return
	}
	fmt.Printf("Generation time (realised, over %d pedigree edges): %.3f +/- %.3f y "+
		"(SD %.3f, CV %.3f, skew %+.2f)\n", s.Edges, s.Mean, s.SEM, s.SD, s.CV, s.Skew)
	fmt.Printf("  median %.1f y, q25-q75 %.1f-%.1f, p5/p95 %.1f/%.1f, max %.1f\n",
		s.Median, s.Q25, s.Q75, s.P5, s.P95, s.Max)
	if s.PrePlateau > 0 {
		fmt.Printf("  NOT stationary: generations 1-%d run %.2f y, generation %d on %.2f y (%.3fx)\n",
			s.PlateauFrom-1, s.PrePlateau, s.PlateauFrom, s.Plateau, s.Drift)
	} else if s.Plateau > 0 {
		fmt.Printf("  stationary from generation %d on at %.2f y\n", s.PlateauFrom, s.Plateau)
	}
	if s.DepthN > 0 {
		fmt.Printf("  living population is %.2f +/- %.2f generations from founding\n",
			s.MeanDepth, s.SDDepth)
	}
	for _, msg := range s.Warnings() {
		fmt.Printf("WARNING: generation time: %s\n", msg)
	}
}

// momentStats returns the mean, population SD and skewness of x.
func momentStats(x []float64) (mean, sd, skew float64) {
	if len(x) == 0 {
		return 0, 0, 0
	}
	n := float64(len(x))
	for _, v := range x {
		mean += v
	}
	mean /= n
	var m2, m3 float64
	for _, v := range x {
		d := v - mean
		m2 += d * d
		m3 += d * d * d
	}
	m2 /= n
	m3 /= n
	sd = math.Sqrt(m2)
	if m2 > 0 {
		skew = m3 / math.Pow(m2, 1.5)
	}
	return mean, sd, skew
}

// gapQuantiles returns min, p5, q25, median, q75, p95, max of x. x is not modified.
func gapQuantiles(x []float64) (mn, p5, q25, med, q75, p95, mx float64) {
	if len(x) == 0 {
		return 0, 0, 0, 0, 0, 0, 0
	}
	s := append([]float64(nil), x...)
	sort.Float64s(s)
	q := func(p float64) float64 { return s[int(p*float64(len(s)-1))] }
	return s[0], q(0.05), q(0.25), q(0.50), q(0.75), q(0.95), s[len(s)-1]
}
