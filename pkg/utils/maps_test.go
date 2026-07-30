package utils

import (
	"bytes"
	"drift/pkg/core"
	"testing"
)

// straitModel is a single-row map Land|Water|Land|Land. An individual at column 0
// can, with wander=2, draw a destination at column 2 (Land) whose straight path
// crosses the Water barrier at column 1 — the canonical "jumpable strait" that the
// crossing test should close.
func straitModel(barriers string) *core.Model {
	return &core.Model{
		Parameters:   map[string]float64{"wander": 2},
		StringParams: map[string]string{"movement_barriers": barriers},
		Map:          [][]int{{1, 3, 1, 1}},
	}
}

// TestWanderUnsetByteIdentical proves that with no movement_barriers table Wander's
// move-accept predicate is exactly the pre-barrier one (IsHabitable only): over a
// seeded sequence of calls it produces identical positions AND leaves the RNG in the
// identical state as an inline legacy oracle that draws the same two offsets.
func TestWanderUnsetByteIdentical(t *testing.T) {
	m := straitModel("") // barriers inactive
	const n = 300

	// Oracle: the original Wander body (two offset draws, accept iff habitable).
	SeedRNG(12345)
	oracleLat, oracleLon := 0, 0
	maxDistance := int(m.Parameters["wander"])
	for i := 0; i < n; i++ {
		offLat := RandIntn(2*maxDistance+1) - maxDistance
		offLon := RandIntn(2*maxDistance+1) - maxDistance
		nl, no := oracleLat+offLat, oracleLon+offLon
		if m.IsHabitable(nl, no) {
			oracleLat, oracleLon = nl, no
		}
	}
	oracleState := RNGState()

	// Real Wander over the same seed must match position and RNG state exactly.
	SeedRNG(12345)
	lat, lon := 0, 0
	for i := 0; i < n; i++ {
		lat, lon = Wander(m, lat, lon)
	}
	if lat != oracleLat || lon != oracleLon {
		t.Errorf("Wander unset diverged: got (%d,%d), oracle (%d,%d)", lat, lon, oracleLat, oracleLon)
	}
	if !bytes.Equal(RNGState(), oracleState) {
		t.Error("Wander unset consumed a different RNG stream than the legacy oracle")
	}
}

// TestWanderBarrierBlocksStrait is the distinguishable-outcome case: the same seed
// and the same map, with vs without a Water barrier. Without the barrier the
// individual jumps the strait and reaches the far shore (column >= 2); with the
// barrier it is pinned behind the strait at its origin. Crucially the barrier adds no
// RNG draws, so both runs end on the identical RNG state.
func TestWanderBarrierBlocksStrait(t *testing.T) {
	const n = 500

	run := func(barriers string) (maxLon, finalLat, finalLon int, state []byte) {
		m := straitModel(barriers)
		SeedRNG(99)
		lat, lon := 0, 0
		for i := 0; i < n; i++ {
			lat, lon = Wander(m, lat, lon)
			if lon > maxLon {
				maxLon = lon
			}
		}
		return maxLon, lat, lon, RNGState()
	}

	offMax, _, _, offState := run("")
	onMax, onLat, onLon, onState := run("3")

	if offMax < 2 {
		t.Fatalf("without a barrier the individual should cross the strait (maxLon>=2), got %d", offMax)
	}
	if onMax != 0 {
		t.Errorf("with a Water barrier the individual should stay pinned at origin (maxLon=0), got %d", onMax)
	}
	if onLat != 0 || onLon != 0 {
		t.Errorf("barrier-pinned individual ended at (%d,%d), want (0,0)", onLat, onLon)
	}
	if !bytes.Equal(offState, onState) {
		t.Error("barrier changed the RNG stream; it must only change accept/reject, not draws")
	}
}
