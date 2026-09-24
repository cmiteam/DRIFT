package analysis

import (
	"drift/pkg/core"
	"math"
	"strings"
	"testing"
)

// Unit tests for the measured generation time (gentime.go). The existing offspring-
// weighting and no-pedigree contracts are pinned in tmrka_ne_test.go against
// PedigreeGenerationTime, which is now a wrapper; these pin the spread, the
// non-stationarity map and the two ways the constant is allowed to be chosen.

// A three-generation chain with KNOWN gaps, so every moment can be checked by hand.
// Founders 1 (y0) and 2 (y0); child 3 at y20 (gaps 20, 20); child 4 of 3 and a second
// founder at y60 (gap 40 from 3).
func chainPop() *core.Pop {
	pop := &core.Pop{IndData: map[int][]int{}}
	pop.PedEnsure(1, 0, 0)
	pop.PedEnsure(2, 0, 1)
	pop.PedRecordBirth(3, 1, 2, 20, 1)
	pop.PedEnsure(5, 10, 0)
	pop.PedRecordBirth(4, 5, 3, 60, 0)
	return pop
}

// The SD reported is the spread of individual gaps; the SEM is that divided by sqrt(n).
// Conflating them is the mistake the whole struct exists to prevent, so both are pinned.
func TestGenTimeSpreadIsGapSDAndSEMIsSmaller(t *testing.T) {
	s := PedigreeGenerationStats(chainPop(), nil)
	if s == nil {
		t.Fatal("no stats from a pedigree with edges")
	}
	// Gaps: 20, 20 (child 3) and 50, 40 (child 4, from founder 5 at y10 and parent 3 at y20).
	if s.Edges != 4 {
		t.Fatalf("edges: got %d, want 4", s.Edges)
	}
	if !closeTo(s.Mean, 32.5, 1e-9) {
		t.Errorf("mean: got %.4f, want 32.5 = mean(20,20,50,40)", s.Mean)
	}
	wantSD := math.Sqrt((12.5*12.5 + 12.5*12.5 + 17.5*17.5 + 7.5*7.5) / 4)
	if !closeTo(s.SD, wantSD, 1e-9) {
		t.Errorf("SD: got %.4f, want %.4f", s.SD, wantSD)
	}
	if !closeTo(s.SEM, wantSD/2, 1e-9) {
		t.Errorf("SEM: got %.4f, want SD/sqrt(4) = %.4f", s.SEM, wantSD/2)
	}
	if s.SEM >= s.SD {
		t.Error("SEM must be smaller than the gap SD; they are not interchangeable")
	}
}

// The per-generation map is the whole answer to a life history whose generation time
// moves. Generation 1 must hold the founders' children and generation 2 their
// grandchildren, each with its own local rate.
func TestGenTimeBinsAreTheGenerationsToYearsMap(t *testing.T) {
	s := PedigreeGenerationStats(chainPop(), nil)
	if len(s.Bins) != 2 {
		t.Fatalf("bins: got %d, want 2 (children and grandchildren)", len(s.Bins))
	}
	if s.Bins[0].Gen != 1 || !closeTo(s.Bins[0].MeanGap, 20, 1e-9) {
		t.Errorf("generation 1: got gen=%d gap=%.4f, want 1 and 20", s.Bins[0].Gen, s.Bins[0].MeanGap)
	}
	if s.Bins[1].Gen != 2 || !closeTo(s.Bins[1].MeanGap, 45, 1e-9) {
		t.Errorf("generation 2: got gen=%d gap=%.4f, want 2 and 45 = mean(50,40)",
			s.Bins[1].Gen, s.Bins[1].MeanGap)
	}
	if !closeTo(s.Bins[1].MeanBirthYear, 60, 1e-9) {
		t.Errorf("generation 2 sits at year %.4f, want 60", s.Bins[1].MeanBirthYear)
	}
}

// Beyond the pedigree there is no measured generation time, and YearsForGens must SAY so
// rather than quietly extrapolate — that regime is exactly where an over-deep ARG
// inference lands, and a caller has to know the conversion left the data behind.
func TestYearsForGensFlagsExtrapolationPastThePedigree(t *testing.T) {
	pop := chainPop()
	s := PedigreeGenerationStats(pop, nil)
	s.SampleYear = 60

	if _, extrap := s.YearsForGens(1); extrap {
		t.Error("one generation back is inside the pedigree; must not be flagged")
	}
	years, extrap := s.YearsForGens(500)
	if !extrap {
		t.Fatal("500 generations is far past the pedigree and must be flagged extrapolated")
	}
	if years <= 500 {
		t.Errorf("extrapolated depth %.2f y is not even 1 y per generation", years)
	}
}

// genPop builds a pedigree of discrete generations — four individuals each, both parents
// drawn from the generation before — whose successive parent-child gaps are exactly the
// years given. Every member of generation g sits at depth exactly g, which a lineage
// married into fresh founders each generation would NOT: pedigree collapse drags the
// mean-of-parents depth back towards 1 and smears the generations together.
//
// sampleYear is left far past the last birth so nothing is right-censored unless a test
// asks for it.
func genPop(gaps ...int) *core.Pop {
	const width = 4
	pop := &core.Pop{IndData: map[int][]int{}}
	prev := make([]int, width)
	for j := 0; j < width; j++ {
		prev[j] = j + 1
		pop.PedEnsure(prev[j], 0, j%2)
	}
	next := make([]int, width)
	id, year := 100, 0
	for _, gap := range gaps {
		year += gap
		for j := 0; j < width; j++ {
			id++
			pop.PedRecordBirth(id, prev[(2*j)%width], prev[(2*j+1)%width], year, j%2)
			next[j] = id
		}
		copy(prev, next)
	}
	return pop
}

// A decaying life history must be caught, and a constant one must not be flagged. The
// split has to be made over GENERATIONS: the founder transient carries very few edges, so
// a calendar-year or edge-weighted split averages it away and calls a patriarchal model
// stationary when it is not.
func TestNonStationarityIsDetectedByGenerationNotByYear(t *testing.T) {
	steady := PedigreeGenerationStats(genPop(20, 20, 20, 20, 20, 20, 20, 20), nil)
	if steady.PrePlateau > 0 {
		t.Errorf("a constant life history has no founder transient; got %.3f y over %d edges",
			steady.PrePlateau, steady.PrePlateauEdges)
	}
	if steady.PlateauFrom != 1 {
		t.Errorf("a constant life history is settled from generation 1; got %d", steady.PlateauFrom)
	}
	if hasGenWarning(steady.Warnings(), "NOT stationary") {
		t.Errorf("a constant life history must not warn: %v", steady.Warnings())
	}

	// Patriarchal decay: a long founder transient collapsing to a settled 20 years.
	decay := PedigreeGenerationStats(genPop(200, 100, 50, 20, 20, 20, 20, 20), nil)
	if decay.Drift <= 1.15 {
		t.Fatalf("a decaying life history must register as non-stationary; drift %.3f", decay.Drift)
	}
	if decay.PlateauFrom < 3 {
		t.Errorf("the transient runs at least two generations; plateau starts at %d", decay.PlateauFrom)
	}
	if !closeTo(decay.Plateau, 20, 1e-9) {
		t.Errorf("settled rate: got %.4f, want 20", decay.Plateau)
	}
	if !hasGenWarning(decay.Warnings(), "NOT stationary") {
		t.Errorf("no non-stationarity warning: %v", decay.Warnings())
	}
	// And the whole point of splitting by generation: the plain mean is dragged well
	// above the settled rate by a transient that is a small minority of the edges.
	if decay.Mean <= decay.Plateau {
		t.Errorf("mean %.3f should exceed the settled %.3f on a decaying history",
			decay.Mean, decay.Plateau)
	}
	if decay.Recommended() != decay.Plateau {
		t.Error("the recommended constant must be the settled rate, not the transient-inflated mean")
	}
}

// Right-censoring at the present end: the youngest generations' children have not had
// their own children yet, so their apparent gap is short. Counting it would drag the
// settled rate down and make the life history look as though it sped up.
func TestRightCensoredGenerationsAreExcludedFromTheSettledRate(t *testing.T) {
	pop := genPop(40, 40, 40, 40, 40)
	s := PedigreeGenerationStats(pop, nil)
	// Sample immediately after the last birth: that generation is censored.
	s.SampleYear = 200
	s.findPlateau()
	if s.CensoredFrom == 0 {
		t.Fatal("the generation born at the sample year must be flagged right-censored")
	}
	if !hasGenWarning(s.Warnings(), "right-censored") {
		t.Errorf("no censoring warning: %v", s.Warnings())
	}
	if !closeTo(s.Plateau, 40, 1e-9) {
		t.Errorf("settled rate: got %.4f, want 40 — censored bins must not drag it", s.Plateau)
	}
}

// The declared parameter disagreeing with the measurement is the headline failure mode,
// and the warning has to name the factor and the flag that fixes it.
func TestDeclaredGenerationTimeDisagreementIsWarned(t *testing.T) {
	// A stationary 30-year life history, sampled long after the last birth so nothing is
	// censored: the realised constant is unambiguously 30.
	pop := genPop(30, 30, 30, 30, 30, 30)
	model := &core.Model{
		Parameters:     map[string]float64{"generation_time": 25},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"year": 1000},
	}
	s := PedigreeGenerationStats(pop, model)
	if !closeTo(s.Recommended(), 30, 1e-9) {
		t.Fatalf("realised constant: got %.4f, want 30", s.Recommended())
	}
	if !hasGenWarning(s.Warnings(), "argweaver-gen-time") {
		t.Errorf("30 y realised against a 25 y declared must warn and name the flag: %v",
			s.Warnings())
	}

	model.Parameters["generation_time"] = 30
	s = PedigreeGenerationStats(pop, model)
	if hasGenWarning(s.Warnings(), "argweaver-gen-time") {
		t.Errorf("an agreeing declaration must not warn: %v", s.Warnings())
	}
}

// The conversion constant stays the DECLARED one unless asked, because switching it
// silently would rewrite the gens_ago column of every archived run.
func TestResolveGenerationTimeDefaultsToDeclared(t *testing.T) {
	model := &core.Model{
		Parameters:     map[string]float64{"generation_time": 25},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"year": 60},
	}
	pop := chainPop()

	got, stats := ResolveGenerationTime(model, pop)
	if !closeTo(got, 25, 1e-9) {
		t.Errorf("default: got %.4f, want the declared 25", got)
	}
	if stats == nil || !closeTo(stats.Mean, 32.5, 1e-9) {
		t.Error("the measurement must be reported even when it is not used")
	}

	model.StringParams["generation_time_source"] = "pedigree"
	if got, _ = ResolveGenerationTime(model, pop); !closeTo(got, 32.5, 1e-9) {
		t.Errorf("opted in: got %.4f, want the realised 32.5", got)
	}

	// Opting in without a pedigree must fall back, not divide by zero.
	if got, _ = ResolveGenerationTime(model, &core.Pop{}); !closeTo(got, 25, 1e-9) {
		t.Errorf("no pedigree: got %.4f, want the declared 25 as a fallback", got)
	}
}

// A nil GenTimeStats is the honest answer when track_pedigree was off, and every
// consumer must survive it — including the report writer, which has to say so in the
// file rather than omit the section.
func TestNilGenTimeStatsIsUsableEverywhere(t *testing.T) {
	var s *GenTimeStats
	if got, edges := PedigreeGenerationTime(nil); got != 0 || edges != 0 {
		t.Error("a nil pop must report nothing, not a made-up constant")
	}
	if w := s.Warnings(); len(w) == 0 || !strings.Contains(w[0], "track_pedigree") {
		t.Errorf("a missing measurement must warn and name the flag: %v", w)
	}
	var b strings.Builder
	s.WriteReport(&b)
	if !strings.Contains(b.String(), "NOT MEASURED") {
		t.Errorf("the report must say the section is unmeasured, got %q", b.String())
	}
	if err := WriteGenTimeCSV("", s); err != nil {
		t.Errorf("writing a nil map must be a no-op, got %v", err)
	}
}

func hasGenWarning(ws []string, substr string) bool {
	for _, w := range ws {
		if strings.Contains(w, substr) {
			return true
		}
	}
	return false
}
