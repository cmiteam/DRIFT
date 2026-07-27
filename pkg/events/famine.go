package events

import (
	"drift/pkg/core"
	"fmt"
)

// famine is a windowed environmental-mortality event (roadmap §2): during the
// inclusive year span [start, end] it multiplies every individual's actuarial
// death risk by `mortality`. This raises the number of deaths each year WITHOUT
// changing the RNG stream — deathStandard (Step 1) still draws exactly one death
// roll per individual, in the same sorted-id order, and only the death threshold
// moves. So a scheduled famine perturbs the trajectory while an unscheduled run is
// byte-identical, and even a scheduled run stays deterministic.
//
// It composes with the rest of death.go: the bottleneck window (Step 2's hard
// bottleneck_size), the logistic carrying_capacity recovery (Step 3), and the
// per-deme cull (§6c) all act elsewhere, so a famine can thin a population that a
// logistic K then regrows, or deepen a bottleneck.
//
// Spec: famine start=<year> end=<year> mortality=<factor>
//   - start   (required) first year the famine acts
//   - end     (optional, default = start) last year it acts; must be >= start
//   - mortality (required, > 0) death-risk multiplier; > 1 is a famine, a value in
//     (0,1) is a boon (lowered mortality). The per-individual adjusted risk is
//     clamped by the usual roll, so an extreme factor simply kills everyone eligible.
type famine struct {
	start, end int
	mortality  float64
}

func newFamine(fields map[string]string) (Event, error) {
	start, err := reqInt(fields, "start")
	if err != nil {
		return nil, fmt.Errorf("famine: %w", err)
	}
	end, err := optInt(fields, "end", start)
	if err != nil {
		return nil, fmt.Errorf("famine: %w", err)
	}
	if end < start {
		return nil, fmt.Errorf("famine: end (%d) is before start (%d)", end, start)
	}
	mortality, err := reqFloat(fields, "mortality")
	if err != nil {
		return nil, fmt.Errorf("famine: %w", err)
	}
	if mortality <= 0 {
		return nil, fmt.Errorf("famine: mortality must be > 0, got %g", mortality)
	}
	return &famine{start: start, end: end, mortality: mortality}, nil
}

// Apply raises the model's transient mortality factor during the famine window.
// The Schedule resets the factor to 1.0 each year before events run, so this is a
// clean per-year multiply (overlapping famines compound multiplicatively).
func (f *famine) Apply(model *core.Model, pop *core.Pop, year int) {
	if year >= f.start && year <= f.end {
		model.EventMortalityFactor *= f.mortality
	}
}

func init() { RegisterEvent("famine", newFamine) }
