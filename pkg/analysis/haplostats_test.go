package analysis

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"math"
	"testing"
)

// ---- fixtures ----

type indSpec struct {
	deme    int
	strand0 []int
	strand1 []int
}

// buildHaploPop builds a Pop with the given individuals (each two strands of
// mutation ids) and a MutationPool from mutPos/mutEff. IDs are 1..n.
func buildHaploPop(inds []indSpec, mutPos map[int]int, mutEff map[int]float64) *core.Pop {
	pop := &core.Pop{
		IndData:      map[int][]int{},
		IndMutations: map[int]map[int][]int{},
		MutationPool: map[int]core.Mutation{},
	}
	for i, s := range inds {
		id := i + 1
		row := individual.MakeIndData()
		row[individual.Deme] = s.deme
		pop.IndData[id] = row
		pop.IndMutations[id] = map[int][]int{0: append([]int{}, s.strand0...), 1: append([]int{}, s.strand1...)}
	}
	for mid, pos := range mutPos {
		pop.MutationPool[mid] = core.Mutation{Id: mid, Position: pos, Effect: mutEff[mid]}
	}
	return pop
}

// oneChromModel: chromosome 0 = arm0 [0,50) cM=5, arm1 [50,100) cM=5 (so
// cM/bit = 0.1 everywhere: CumCM(50)=5, CumCM(60)=6, CumCM(70)=7).
func oneChromModel() *core.Model {
	m := &core.Model{
		Parameters:     map[string]float64{"recomb_cM_per_bit": 0.1},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"genome_bits": 100},
		ChromosomeArms: map[int]map[int][]int{0: {0: {0, 50}, 1: {50, 50}}},
		ArmCM:          map[int]map[int]float64{0: {0: 5.0, 1: 5.0}},
	}
	return m
}

func findCore(cores []CoreEHH, mutID int) *CoreEHH {
	for i := range cores {
		if cores[i].MutID == mutID {
			return &cores[i]
		}
	}
	return nil
}

func approxHS(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s = %.6f, want %.6f (tol %.4g)", name, got, want, tol)
	}
}

// ---- tests ----

func TestHomozygosity(t *testing.T) {
	// 2 groups of 2 out of 4 => (2*1 + 2*1)/(4*3) = 4/12.
	got := homozygosity(map[int]int{0: 2, 1: 2}, 4)
	approxHS(t, "homozygosity", got, 4.0/12.0, 1e-12)
	// all identical => 1.0
	approxHS(t, "homozygosity-all", homozygosity(map[int]int{0: 5}, 5), 1.0, 1e-12)
}

// Hand-computed EHH/iHH on a designed fixture (see oneChromModel geometry).
//
// Core M @ pos50 (cM 5). Derived carriers = 4 strands; site A @ pos60 (cM 6)
// carried by ALL 4 derived; site B @ pos70 (cM 7) carried by 2 of the 4.
// Downstream derived EHH: 1.0 at core, 1.0 at A (no split), 4/12=0.333 at B.
//
//	iHH_D(+) = (1+1)/2*(6-5) + (1+0.333..)/2*(7-6) = 1.0 + 0.6667 = 1.6667
//
// Ancestral 4 strands carry neither A nor B => EHH stays 1.0:
//
//	iHH_A(+) = 1.0 + 1.0 = 2.0
//
// No sites upstream of the core => both upstream integrals 0.
//
//	UnstdIHS = ln(2.0/1.6667) = 0.18232
func TestEHHHandComputed(t *testing.T) {
	// ids: 1..4 carry M (derived), 5..8 do not (ancestral). One id = one diploid
	// (2 strands). To make counts clean, put both a strand's copies distinct.
	mutPos := map[int]int{100: 50, 200: 60, 300: 70}
	mutEff := map[int]float64{}
	inds := []indSpec{
		// derived individuals: each strand carries M(100). strand structure below
		// gives 4 derived strands carrying A(200); 2 of them also carry B(300).
		{0, []int{100, 200, 300}, []int{100, 200, 300}}, // id1: both strands carry M,A,B
		{0, []int{100, 200}, []int{100, 200}},           // id2: both carry M,A (not B)
		// ancestral individuals: carry neither M,A,B
		{0, []int{}, []int{}}, // id3
		{0, []int{}, []int{}}, // id4
	}
	// That is 4 derived strands (id1x2, id2x2), of which 2 (id1's) carry B.
	// Ancestral strands = 4 (id3,id4). nhap=8. M derived=4 (0.5), A derived=4, B=2.
	pop := buildHaploPop(inds, mutPos, mutEff)
	m := oneChromModel()
	ids := []int{1, 2, 3, 4}
	cfg := DefaultHaploConfig()
	cfg.NumFreqBin = 1
	res := ComputeHaploStats(m, pop, ids, cfg)
	if res == nil {
		t.Fatal("nil result")
	}
	c := findCore(res.Cores, 100)
	if c == nil {
		t.Fatalf("core M(100) not found among %d cores", len(res.Cores))
	}
	approxHS(t, "iHH_derived", c.IHHDerived, 1.6667, 1e-3)
	approxHS(t, "iHH_ancestral", c.IHHAncestral, 2.0, 1e-6)
	approxHS(t, "UnstdIHS", c.UnstdIHS, 0.18232, 1e-3)
}

// Distinguishable outcome: a swept core (all derived carriers share a long
// identical downstream block) shows much larger iHH_derived and a more NEGATIVE
// iHS than a neutral core (derived carriers diverse) at the same frequency.
func TestIHSSweepVsNeutral(t *testing.T) {
	// Two chromosomes. Chrom 0 [0,100): a SWEEP. Chrom 1 [100,200): NEUTRAL.
	m := &core.Model{
		Parameters:     map[string]float64{"recomb_cM_per_bit": 0.1},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"genome_bits": 200},
		ChromosomeArms: map[int]map[int][]int{
			0: {0: {0, 50}, 1: {50, 50}},
			1: {0: {100, 50}, 1: {150, 50}},
		},
		ArmCM: map[int]map[int]float64{
			0: {0: 5.0, 1: 5.0},
			1: {0: 5.0, 1: 5.0},
		},
	}
	// mutation positions
	mutPos := map[int]int{}
	// Sweep core S=1000 @55; long shared block of markers 1001..1005 at 60..80 ALL
	// carried by every derived carrier (identity by descent => high EHH far out).
	mutPos[1000] = 55
	block := []int{1001, 1002, 1003, 1004, 1005}
	for i, mid := range block {
		mutPos[mid] = 60 + i*5
	}
	// Neutral core N=2000 @155; derived carriers each carry DISTINCT flanking
	// markers (no shared block => EHH collapses immediately).
	mutPos[2000] = 155
	// distinct neutral flanks per carrier
	nflank := map[int][]int{}
	fid := 3000
	// build 6 derived carriers for each core (freq 6/12 = 0.5)
	var inds []indSpec
	swDerived := [][]int{}
	for k := 0; k < 6; k++ {
		// each derived-for-sweep strand: core + full shared block
		swDerived = append(swDerived, append([]int{1000}, block...))
	}
	// neutral derived strands: core + one distinct private marker each
	neuDerived := [][]int{}
	for k := 0; k < 6; k++ {
		priv := fid
		fid++
		mutPos[priv] = 160 + k // distinct positions
		nflank[k] = []int{priv}
		neuDerived = append(neuDerived, []int{2000, priv})
	}
	// 6 individuals: strand0 = sweep-derived, strand1 = neutral-derived => both
	// cores at freq 6/12 on their chromosomes.
	for k := 0; k < 6; k++ {
		inds = append(inds, indSpec{0, swDerived[k], neuDerived[k]})
	}
	// 6 ancestral individuals: carry neither core (empty strands).
	for k := 0; k < 6; k++ {
		inds = append(inds, indSpec{0, []int{}, []int{}})
	}
	pop := buildHaploPop(inds, mutPos, nil)
	ids := make([]int, len(inds))
	for i := range ids {
		ids[i] = i + 1
	}
	cfg := DefaultHaploConfig()
	cfg.NumFreqBin = 1
	res := ComputeHaploStats(m, pop, ids, cfg)
	if res == nil {
		t.Fatal("nil result")
	}
	sw := findCore(res.Cores, 1000)
	neu := findCore(res.Cores, 2000)
	if sw == nil || neu == nil {
		t.Fatalf("missing cores: sweep=%v neutral=%v (have %d)", sw != nil, neu != nil, len(res.Cores))
	}
	t.Logf("sweep  iHH_D=%.3f iHH_A=%.3f iHS=%.3f", sw.IHHDerived, sw.IHHAncestral, sw.UnstdIHS)
	t.Logf("neutral iHH_D=%.3f iHH_A=%.3f iHS=%.3f", neu.IHHDerived, neu.IHHAncestral, neu.UnstdIHS)
	if sw.IHHDerived <= neu.IHHDerived {
		t.Errorf("sweep iHH_derived (%.3f) should exceed neutral (%.3f)", sw.IHHDerived, neu.IHHDerived)
	}
	if !(sw.UnstdIHS < neu.UnstdIHS) {
		t.Errorf("sweep iHS (%.3f) should be more negative than neutral (%.3f)", sw.UnstdIHS, neu.UnstdIHS)
	}
}

func TestPairR2(t *testing.T) {
	// 4 strands: sites X,Y. X carried by strands {0,1}; Y by {0,1} => perfect LD.
	c := &chromHaps{nhap: 4, member: []map[int]bool{
		{10: true, 20: true}, {10: true, 20: true}, {}, {},
	}}
	approxHS(t, "r2-perfect", pairR2(c, 10, 20, 4), 1.0, 1e-9)
	// Y instead carried by {2,3} => opposite: still |r|=1 (r^2=1).
	c2 := &chromHaps{nhap: 4, member: []map[int]bool{
		{10: true}, {10: true}, {20: true}, {20: true},
	}}
	approxHS(t, "r2-opposite", pairR2(c2, 10, 20, 4), 1.0, 1e-9)
	// independent: X {0,1}, Y {0,2} => pA=.5 pB=.5 pAB=.25 D=0 => r2=0.
	c3 := &chromHaps{nhap: 4, member: []map[int]bool{
		{10: true, 20: true}, {10: true}, {20: true}, {},
	}}
	approxHS(t, "r2-indep", pairR2(c3, 10, 20, 4), 0.0, 1e-9)
}

func TestParseDemePair(t *testing.T) {
	cases := []struct {
		in   string
		a, b int
	}{
		{"", -1, -1},
		{"0,1", 0, 1},
		{"2;3", 2, 3},
		{"1 4", 1, 4},
		{"x", -1, -1},
		{"5", -1, -1},
	}
	for _, tc := range cases {
		a, b := ParseDemePair(tc.in)
		if a != tc.a || b != tc.b {
			t.Errorf("ParseDemePair(%q) = (%d,%d), want (%d,%d)", tc.in, a, b, tc.a, tc.b)
		}
	}
}

func TestStandardizeIHS(t *testing.T) {
	// three cores in the same freq bin with UnstdIHS 1,2,3 => mean 2, sd sqrt(2/3).
	cores := []CoreEHH{
		{DerivedFreq: 0.5, UnstdIHS: 1},
		{DerivedFreq: 0.5, UnstdIHS: 2},
		{DerivedFreq: 0.5, UnstdIHS: 3},
	}
	standardizeIHS(cores, 1)
	var sum float64
	for _, c := range cores {
		if math.IsNaN(c.StdIHS) {
			t.Fatalf("unexpected NaN StdIHS")
		}
		sum += c.StdIHS
	}
	approxHS(t, "std-mean0", sum/3, 0.0, 1e-9)
	approxHS(t, "std-mid", cores[1].StdIHS, 0.0, 1e-9)
}

// XP-EHH: deme 0 carriers share a long block (high iES); deme 1 diverse (low
// iES) => UnstdXPEHH = ln(iESA/iESB) > 0 at the shared site.
func TestXPEHHSign(t *testing.T) {
	m := oneChromModel()
	mutPos := map[int]int{500: 55, 501: 60, 502: 65, 503: 70}
	var inds []indSpec
	// deme 0: 4 individuals, both strands carry the shared block => high iES.
	for k := 0; k < 4; k++ {
		inds = append(inds, indSpec{0, []int{500, 501, 502, 503}, []int{500, 501, 502, 503}})
	}
	// deme 1: 4 individuals; only strand0 carries 500 (+ a distinct private flank),
	// strand1 empty => site 500 stays segregating overall AND deme-1 EHHS is low.
	fid := 600
	for k := 0; k < 4; k++ {
		p := fid
		fid++
		mutPos[p] = 60 + k
		inds = append(inds, indSpec{1, []int{500, p}, []int{}})
	}
	pop := buildHaploPop(inds, mutPos, nil)
	ids := make([]int, len(inds))
	for i := range ids {
		ids[i] = i + 1
	}
	cfg := DefaultHaploConfig()
	cfg.DemeA, cfg.DemeB = 0, 1
	res := ComputeHaploStats(m, pop, ids, cfg)
	if res == nil || res.DemeA != 0 || res.DemeB != 1 {
		t.Fatalf("XP-EHH demes not set: %+v", res)
	}
	var site *XPEHHSite
	for i := range res.XPEHH {
		if res.XPEHH[i].MutID == 500 {
			site = &res.XPEHH[i]
		}
	}
	if site == nil {
		t.Fatalf("XP-EHH site 500 not found (have %d)", len(res.XPEHH))
	}
	t.Logf("iES_A=%.3f iES_B=%.3f XP-EHH=%.3f", site.IESA, site.IESB, site.UnstdXPEHH)
	if !(site.UnstdXPEHH > 0) {
		t.Errorf("XP-EHH (%.3f) should be > 0 (deme A more extended)", site.UnstdXPEHH)
	}
}

func TestComputeHaploStatsGuards(t *testing.T) {
	m := oneChromModel()
	pop := buildHaploPop([]indSpec{{0, []int{1}, []int{1}}}, map[int]int{1: 50}, nil)
	if ComputeHaploStats(m, pop, []int{1}, DefaultHaploConfig()) != nil {
		t.Error("expected nil for <2 samples")
	}
}
