package core

import (
	"math"
	"testing"
)

// A two-arm chromosome map with hand-picked cM so cumulative-cM and its inverse have
// exact, round-numbered expected values:
//
//	arm 0 (p): bits [0,100),   cM = 50   -> 0.50 cM/bit
//	arm 1 (q): bits [100,200), cM = 150  -> 1.50 cM/bit
//
// totalCM = 200. cM/bit differs 3x between arms, so BitAtCM must concentrate
// crossovers in the q arm.
func handMap() *GeneticMap {
	arms := map[int]map[int][]int{1: {0: {0, 100}, 1: {100, 100}}}
	armCM := map[int]map[int]float64{1: {0: 50, 1: 150}}
	return BuildGeneticMap(arms, armCM, 1.0)
}

func approx(a, b float64) bool { return math.Abs(a-b) < 1e-9 }

func TestGeneticMap_TotalAndSpan(t *testing.T) {
	g := handMap()
	if !g.Has(1) {
		t.Fatal("expected chromosome 1 present")
	}
	if g.Has(2) {
		t.Error("chromosome 2 should be absent")
	}
	if !approx(g.TotalCM(1), 200) {
		t.Errorf("TotalCM = %v, want 200", g.TotalCM(1))
	}
	if s, e := g.Span(1); s != 0 || e != 200 {
		t.Errorf("Span = [%d,%d), want [0,200)", s, e)
	}
}

func TestGeneticMap_CumCM(t *testing.T) {
	g := handMap()
	cases := []struct {
		pos  int
		want float64
	}{
		{0, 0},     // chromosome start
		{50, 25},   // halfway through p arm: 50*0.5
		{100, 50},  // p/q boundary: whole p arm
		{150, 125}, // halfway through q arm: 50 + 50*1.5
		{200, 200}, // clamps to totalCM past the end
		{999, 200}, // far past the end
	}
	for _, c := range cases {
		if got := g.CumCM(1, c.pos); !approx(got, c.want) {
			t.Errorf("CumCM(pos=%d) = %v, want %v", c.pos, got, c.want)
		}
	}
	// Monotone non-decreasing across the whole span.
	prev := -1.0
	for p := 0; p <= 200; p++ {
		v := g.CumCM(1, p)
		if v < prev-1e-12 {
			t.Fatalf("CumCM not monotone at pos %d: %v < %v", p, v, prev)
		}
		prev = v
	}
}

func TestGeneticMap_BitAtCM_InvertsCumCM(t *testing.T) {
	g := handMap()
	// Points chosen to land strictly inside a segment so the linear inverse is exact.
	cases := []struct {
		cm   float64
		want int
	}{
		{25, 50},   // p arm midpoint
		{125, 150}, // q arm midpoint
	}
	for _, c := range cases {
		if got := g.BitAtCM(1, c.cm); got != c.want {
			t.Errorf("BitAtCM(cm=%v) = %d, want %d", c.cm, got, c.want)
		}
		// round-trip: the bit's CumCM returns the cM we asked for.
		if back := g.CumCM(1, c.want); !approx(back, c.cm) {
			t.Errorf("round-trip CumCM(BitAtCM(%v)) = %v, want %v", c.cm, back, c.cm)
		}
	}
	// Concentration check: uniform cM draws should place most crossovers in the q arm
	// (3x the rate). Enumerate cM in fine steps and count which arm they land in.
	pArm, qArm := 0, 0
	for cm := 0.0; cm < 200; cm += 0.1 {
		if g.BitAtCM(1, cm) < 100 {
			pArm++
		} else {
			qArm++
		}
	}
	if !(qArm > 2*pArm) {
		t.Errorf("expected q arm to attract ~3x crossovers, got p=%d q=%d", pArm, qArm)
	}
}

func TestBuildGeneticMap_DerivedFallback(t *testing.T) {
	arms := map[int]map[int][]int{1: {0: {0, 100}, 1: {100, 100}}}
	// No cM table at all -> every arm uses length*defaultCMPerBit.
	g := BuildGeneticMap(arms, nil, 2.0)
	if !approx(g.TotalCM(1), 400) {
		t.Errorf("derived TotalCM = %v, want 400 (200 bits * 2 cM/bit)", g.TotalCM(1))
	}
	// Partial table: q arm specified, p arm falls back to the derived rate.
	armCM := map[int]map[int]float64{1: {1: 300}}
	g2 := BuildGeneticMap(arms, armCM, 2.0)
	if !approx(g2.TotalCM(1), 200+300) {
		t.Errorf("partial TotalCM = %v, want 500 (p:100*2 derived + q:300)", g2.TotalCM(1))
	}
}
