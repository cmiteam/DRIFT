package simulation

import (
	"testing"

	"drift/pkg/core"
)

// TestEventMortalityFactor covers the death-path helper that a famine event
// (roadmap §2) drives. The zero-value (no event scheduler wrote a factor) must
// read as 1.0 so the death threshold is unchanged — the byte-identical no-op — and
// a positive value passes through so a famine actually raises mortality.
func TestEventMortalityFactor(t *testing.T) {
	cases := []struct {
		name   string
		factor float64
		want   float64
	}{
		{"zero-value reads neutral", 0, 1.0},
		{"negative reads neutral", -3, 1.0},
		{"famine multiplier passes through", 2.5, 2.5},
		{"boon (sub-1) passes through", 0.5, 0.5},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := &core.Model{EventMortalityFactor: tc.factor}
			if got := eventMortalityFactor(m); got != tc.want {
				t.Errorf("eventMortalityFactor(%g) = %g, want %g", tc.factor, got, tc.want)
			}
		})
	}
}

// TestEventMortalityByteIdenticalWhenNeutral asserts that applying the neutral
// factor to a death risk is an exact identity: deathrisk * 1.0 == deathrisk for
// every risk value, so an unscheduled run's death threshold is bit-for-bit
// unchanged. (The engine additionally guards on EventScheduler != nil, so it never
// even performs this multiply without events — this is the belt-and-braces check
// that the multiply itself is a no-op.)
func TestEventMortalityByteIdenticalWhenNeutral(t *testing.T) {
	m := &core.Model{EventMortalityFactor: 0} // no event wrote a factor
	factor := eventMortalityFactor(m)
	for _, risk := range []float64{0, 0.001, 0.0123456789, 0.5, 0.9999999, 1.0} {
		if risk*factor != risk {
			t.Errorf("deathrisk %g * neutral factor %g = %g, want exactly %g", risk, factor, risk*factor, risk)
		}
	}
}
