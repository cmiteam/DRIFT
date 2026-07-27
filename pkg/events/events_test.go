package events

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"testing"
)

// modelWith builds a bare model whose environmental_events param is set to spec.
func modelWith(spec string) *core.Model {
	return &core.Model{
		Parameters:     map[string]float64{},
		StringParams:   map[string]string{"environmental_events": spec},
		FreeParameters: map[string]int{},
	}
}

// makePop builds a population of n individuals all in deme 0.
func makePop(n int) *core.Pop {
	pop := &core.Pop{IndData: make(map[int][]int)}
	for i := 0; i < n; i++ {
		d := individual.MakeIndData()
		d[individual.Deme] = 0
		pop.IndData[i] = d
	}
	return pop
}

func demeCount(pop *core.Pop, deme int) int {
	n := 0
	for _, d := range pop.IndData {
		if d[individual.Deme] == deme {
			n++
		}
	}
	return n
}

// --- Loader / parser ---------------------------------------------------------

// TestLoadScheduleEmptyIsNoOp is the strict no-op guard: an unset (or all-blank)
// environmental_events param yields no scheduler, so the run takes the byte-
// identical existing path.
func TestLoadScheduleEmptyIsNoOp(t *testing.T) {
	for _, spec := range []string{"", "   ", " ; ;; "} {
		m := modelWith(spec)
		s, err := LoadSchedule(m)
		if err != nil {
			t.Fatalf("LoadSchedule(%q) unexpected error: %v", spec, err)
		}
		if s != nil {
			t.Errorf("LoadSchedule(%q) = non-nil schedule, want nil (no-op)", spec)
		}
	}
	// A model with no StringParams map at all must also be a no-op.
	m := &core.Model{Parameters: map[string]float64{}, FreeParameters: map[string]int{}}
	if s, err := LoadSchedule(m); err != nil || s != nil {
		t.Errorf("LoadSchedule(no StringParams) = (%v, %v), want (nil, nil)", s, err)
	}
}

func TestLoadScheduleParsesBothTypes(t *testing.T) {
	m := modelWith("famine start=200 end=205 mortality=2.0; migration_pulse year=300 from=0 to=1 fraction=0.1")
	s, err := LoadSchedule(m)
	if err != nil {
		t.Fatalf("LoadSchedule error: %v", err)
	}
	if s.Len() != 2 {
		t.Fatalf("parsed %d events, want 2", s.Len())
	}
	fm, ok := s.events[0].(*famine)
	if !ok {
		t.Fatalf("event 0 is %T, want *famine", s.events[0])
	}
	if fm.start != 200 || fm.end != 205 || fm.mortality != 2.0 {
		t.Errorf("famine parsed as %+v, want {200 205 2}", *fm)
	}
	mp, ok := s.events[1].(*migrationPulse)
	if !ok {
		t.Fatalf("event 1 is %T, want *migrationPulse", s.events[1])
	}
	if mp.year != 300 || mp.from != 0 || mp.to != 1 || mp.fraction != 0.1 {
		t.Errorf("migration_pulse parsed as %+v, want {300 0 1 0.1}", *mp)
	}
}

// TestFamineEndDefaultsToStart checks the optional-field inheritance: omitting
// `end` makes the famine a single-year event.
func TestFamineEndDefaultsToStart(t *testing.T) {
	s, err := LoadSchedule(modelWith("famine start=50 mortality=3"))
	if err != nil {
		t.Fatalf("LoadSchedule error: %v", err)
	}
	fm := s.events[0].(*famine)
	if fm.start != 50 || fm.end != 50 {
		t.Errorf("famine end = %d, want 50 (defaults to start)", fm.end)
	}
}

func TestLoadScheduleErrors(t *testing.T) {
	cases := []struct{ name, spec string }{
		{"unknown type", "plague year=10 severity=2"},
		{"famine missing start", "famine mortality=2"},
		{"famine missing mortality", "famine start=10"},
		{"famine non-numeric", "famine start=ten mortality=2"},
		{"famine end before start", "famine start=100 end=50 mortality=2"},
		{"famine mortality zero", "famine start=10 mortality=0"},
		{"malformed field", "famine start=10 mortality"},
		{"pulse same deme", "migration_pulse year=10 from=1 to=1 fraction=0.5"},
		{"pulse fraction too big", "migration_pulse year=10 from=0 to=1 fraction=2"},
		{"pulse fraction zero", "migration_pulse year=10 from=0 to=1 fraction=0"},
		{"pulse missing to", "migration_pulse year=10 from=0 fraction=0.5"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := LoadSchedule(modelWith(tc.spec)); err == nil {
				t.Errorf("LoadSchedule(%q) = nil error, want a parse error", tc.spec)
			}
		})
	}
}

// --- Famine behavior ---------------------------------------------------------

// TestFamineMortalityFactorWindow verifies Apply sets the transient mortality
// factor to the famine multiplier inside the window and resets it to 1.0 outside —
// the distinguishable-outcome vs. no-op behavior the death path reads.
func TestFamineMortalityFactorWindow(t *testing.T) {
	m := modelWith("famine start=10 end=12 mortality=2.5")
	s, err := LoadSchedule(m)
	if err != nil {
		t.Fatalf("LoadSchedule error: %v", err)
	}
	pop := makePop(1)
	cases := []struct {
		year int
		want float64
	}{
		{9, 1.0},  // before: neutral
		{10, 2.5}, // window start
		{12, 2.5}, // window end (inclusive)
		{13, 1.0}, // after: reset to neutral
	}
	for _, tc := range cases {
		m.FreeParameters["year"] = tc.year
		s.Apply(m, pop)
		if m.EventMortalityFactor != tc.want {
			t.Errorf("year %d: EventMortalityFactor = %g, want %g", tc.year, m.EventMortalityFactor, tc.want)
		}
	}
}

// TestOverlappingFaminesCompound checks that two famines active in the same year
// multiply, and that the factor is still reset cleanly the year they both end.
func TestOverlappingFaminesCompound(t *testing.T) {
	m := modelWith("famine start=10 end=20 mortality=2; famine start=15 end=25 mortality=3")
	s, _ := LoadSchedule(m)
	pop := makePop(1)

	m.FreeParameters["year"] = 17 // both active
	s.Apply(m, pop)
	if m.EventMortalityFactor != 6 {
		t.Errorf("both famines active: factor = %g, want 6 (2*3)", m.EventMortalityFactor)
	}
	m.FreeParameters["year"] = 30 // neither active
	s.Apply(m, pop)
	if m.EventMortalityFactor != 1.0 {
		t.Errorf("after both windows: factor = %g, want 1.0 (reset)", m.EventMortalityFactor)
	}
}

// --- Migration pulse behavior ------------------------------------------------

// TestMigrationPulseRelabels checks a pulse moves round(fraction*|from|) members
// only in its firing year, and is a no-op in every other year.
func TestMigrationPulseRelabels(t *testing.T) {
	utils.SeedRNG(1)
	m := modelWith("migration_pulse year=5 from=0 to=1 fraction=0.25")
	s, _ := LoadSchedule(m)
	pop := makePop(100)

	// A non-firing year moves nobody.
	m.FreeParameters["year"] = 4
	s.Apply(m, pop)
	if demeCount(pop, 1) != 0 {
		t.Fatalf("year 4 (not firing): deme 1 has %d, want 0", demeCount(pop, 1))
	}

	// The firing year moves 25 of 100.
	m.FreeParameters["year"] = 5
	s.Apply(m, pop)
	if got := demeCount(pop, 1); got != 25 {
		t.Errorf("pulse moved %d to deme 1, want 25 (0.25*100)", got)
	}
	if got := demeCount(pop, 0); got != 75 {
		t.Errorf("deme 0 has %d, want 75", got)
	}
}

// TestMigrationPulseDeterministic confirms the pulse consumes RNG deterministically
// (sorted-id partial Fisher-Yates), so a fixed seed reproduces the exact set moved.
func TestMigrationPulseDeterministic(t *testing.T) {
	pick := func() map[int]bool {
		utils.SeedRNG(777)
		m := modelWith("migration_pulse year=5 from=0 to=1 fraction=0.3")
		s, _ := LoadSchedule(m)
		pop := makePop(200)
		m.FreeParameters["year"] = 5
		s.Apply(m, pop)
		moved := map[int]bool{}
		for id, d := range pop.IndData {
			if d[individual.Deme] == 1 {
				moved[id] = true
			}
		}
		return moved
	}
	a, b := pick(), pick()
	if len(a) != 60 {
		t.Fatalf("moved %d, want 60 (0.3*200)", len(a))
	}
	for id := range a {
		if !b[id] {
			t.Fatalf("migration pulse not deterministic under a fixed seed")
		}
	}
}

// TestScheduleNoRNGWhenIdle is the core reproducibility guarantee: in a year no
// event acts, Apply must draw NO randomness, so the RNG stream (and hence the whole
// run) is identical to a no-event baseline up to the first firing event. We probe
// this by checking the RNG state is untouched across an idle Apply.
func TestScheduleNoRNGWhenIdle(t *testing.T) {
	utils.SeedRNG(42)
	m := modelWith("famine start=100 end=110 mortality=2; migration_pulse year=200 from=0 to=1 fraction=0.5")
	s, _ := LoadSchedule(m)
	pop := makePop(50)

	before := utils.RNGState()
	m.FreeParameters["year"] = 5 // neither event active
	s.Apply(m, pop)
	after := utils.RNGState()

	if len(before) != len(after) {
		t.Fatalf("RNG state length changed: %d -> %d", len(before), len(after))
	}
	for i := range before {
		if before[i] != after[i] {
			t.Fatalf("idle Apply consumed RNG (state byte %d changed) — breaks byte-identity", i)
		}
	}
	// Famine still sets the neutral factor even when idle.
	if m.EventMortalityFactor != 1.0 {
		t.Errorf("idle year factor = %g, want 1.0", m.EventMortalityFactor)
	}
}

// TestAvailableEvents documents the registered set (guards against an accidental
// unregister).
func TestAvailableEvents(t *testing.T) {
	got := AvailableEvents()
	want := map[string]bool{"famine": true, "migration_pulse": true}
	if len(got) != len(want) {
		t.Fatalf("AvailableEvents() = %v, want keys %v", got, want)
	}
	for _, n := range got {
		if !want[n] {
			t.Errorf("unexpected registered event %q", n)
		}
	}
}
