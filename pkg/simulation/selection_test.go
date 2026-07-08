package simulation

import (
	"math"
	"testing"

	"drift/pkg/core"
)

func modelWithSelection(trackMutations int, mode string) *core.Model {
	m := &core.Model{
		Parameters:   map[string]float64{"track_mutations": float64(trackMutations)},
		StringParams: map[string]string{},
	}
	if mode != "" {
		m.StringParams["selection_mode"] = mode
	}
	return m
}

func TestSelectionMode(t *testing.T) {
	cases := []struct {
		name          string
		track         int
		mode          string
		wantMode      string
		wantFecundity bool
		wantViability bool
	}{
		{"mutations off collapses to none", 0, "both", selNone, false, false},
		{"default is fecundity", 1, "", selFecundity, true, false},
		{"explicit fecundity", 1, "fecundity", selFecundity, true, false},
		{"viability only", 1, "viability", selViability, false, true},
		{"both", 1, "both", selBoth, true, true},
		{"none disables selection", 1, "none", selNone, false, false},
		{"unknown falls back to fecundity", 1, "bogus", selFecundity, true, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := modelWithSelection(tc.track, tc.mode)
			if got := selectionMode(m); got != tc.wantMode {
				t.Errorf("selectionMode = %q, want %q", got, tc.wantMode)
			}
			if got := fecunditySelection(m); got != tc.wantFecundity {
				t.Errorf("fecunditySelection = %v, want %v", got, tc.wantFecundity)
			}
			if got := viabilitySelection(m); got != tc.wantViability {
				t.Errorf("viabilitySelection = %v, want %v", got, tc.wantViability)
			}
		})
	}
}

func TestViabilityHazardFactor(t *testing.T) {
	cases := []struct {
		fitness float64
		want    float64
	}{
		{1.0, 1.0}, // neutral: no effect on hazard
		{0.9, 1.1}, // deleterious load: higher hazard
		{0.5, 1.5}, // strongly deleterious
		{1.2, 0.8}, // beneficial: lower hazard
		{3.0, 0.0}, // pathological super-fit: floored at 0, never negative
	}
	for _, tc := range cases {
		if got := viabilityHazardFactor(tc.fitness); math.Abs(got-tc.want) > 1e-12 {
			t.Errorf("viabilityHazardFactor(%.2f) = %.4f, want %.4f", tc.fitness, got, tc.want)
		}
	}
}
