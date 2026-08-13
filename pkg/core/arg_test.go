package core

import "testing"

func eqInts(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func TestParseARGLoci(t *testing.T) {
	cases := []struct {
		spec       string
		genomeBits int
		want       []int
	}{
		{"", 1000, []int{}},
		{"100", 1000, []int{100}},
		{"300,100,200", 1000, []int{100, 200, 300}},
		{"100 200\t300", 1000, []int{100, 200, 300}},
		{"100,100,100", 1000, []int{100}},          // deduplicated
		{"10-14", 1000, []int{10, 11, 12, 13, 14}}, // inclusive range
		{"0-20:5", 1000, []int{0, 5, 10, 15, 20}},  // stepped range
		{"0-9:4,100", 1000, []int{0, 4, 8, 100}},   // mixed
		{"50,5000", 1000, []int{50}},               // out of range dropped
		{"-7,50", 1000, []int{50}},                 // negative dropped
		{"900-1100:50", 1000, []int{900, 950}},     // range clipped at genome_bits
		{"abc,50", 1000, []int{50}},                // unparseable item skipped
		{"20-10", 1000, []int{}},                   // inverted range skipped
		{"0-20:0", 1000, []int{}},                  // zero step skipped
		{"7", 0, []int{7}},                         // genome_bits unknown ⇒ no upper clip
		{"1000-1004", 0, []int{1000, 1001, 1002, 1003, 1004}},
	}
	for _, c := range cases {
		got := ParseARGLoci(c.spec, c.genomeBits)
		if got == nil {
			t.Errorf("ParseARGLoci(%q) returned nil; want a non-nil slice so the cache is populated", c.spec)
			continue
		}
		if !eqInts(got, c.want) {
			t.Errorf("ParseARGLoci(%q, %d) = %v, want %v", c.spec, c.genomeBits, got, c.want)
		}
	}
}

func TestEnsureARGLociCaches(t *testing.T) {
	m := &Model{
		Parameters:     map[string]float64{"track_arg": 1},
		StringParams:   map[string]string{"arg_loci": "10,20,30"},
		FreeParameters: map[string]int{"genome_bits": 100},
	}
	first := EnsureARGLoci(m)
	if !eqInts(first, []int{10, 20, 30}) {
		t.Fatalf("got %v", first)
	}
	// Changing the param after the cache is warm must not silently re-parse: the loci
	// set is fixed for the run, which is what makes a locus sweep reproducible.
	m.StringParams["arg_loci"] = "99"
	if second := EnsureARGLoci(m); !eqInts(second, first) {
		t.Errorf("cache was bypassed: got %v, want %v", second, first)
	}
	if !ARGFocalOn(m) {
		t.Error("ARGFocalOn should be true with track_arg=1 and a non-empty locus set")
	}
}

func TestARGFocalOffCases(t *testing.T) {
	mk := func(track float64, spec string) *Model {
		return &Model{
			Parameters:     map[string]float64{"track_arg": track},
			StringParams:   map[string]string{"arg_loci": spec},
			FreeParameters: map[string]int{"genome_bits": 100},
		}
	}
	if ARGFocalOn(mk(0, "10,20")) {
		t.Error("track_arg=0 must leave capture off even with loci configured")
	}
	if ARGFocalOn(mk(1, "")) {
		t.Error("an empty arg_loci must leave capture off even with track_arg=1")
	}
	if ARGFocalOn(nil) {
		t.Error("nil model must be off")
	}
}

// The stored bit must be the inverse of the mask bit, because meiosis() takes parental
// copy 0 where the mask is SET. Rather than restate that, build a mask, record it, and
// check the recorded strand against the position's actual origin.
func TestARGCaptureMatchesMaskConvention(t *testing.T) {
	loci := []int{0, 5, 63, 64, 130}
	// Mask with bits 5 and 64 set: those two positions come from parental copy 0.
	mask := make([]uint64, 3)
	mask[0] |= 1 << 5
	mask[1] |= 1 << 0 // genome bit 64
	// A distinct mask for the maternal gamete: bits 0 and 130.
	maskMom := make([]uint64, 3)
	maskMom[0] |= 1 << 0
	maskMom[2] |= 1 << 2 // genome bit 130

	p := &Pop{}
	p.ARGCapture(7, loci, mask, maskMom)

	wantDad := []int{1, 0, 1, 0, 1} // strand 0 exactly where the mask bit is set
	wantMom := []int{0, 1, 1, 1, 0}
	for i := range loci {
		gotDad, ok := p.ARGStrand(7, 0, i, len(loci))
		if !ok {
			t.Fatalf("locus %d: no paternal record", loci[i])
		}
		if gotDad != wantDad[i] {
			t.Errorf("locus %d paternal: got strand %d, want %d", loci[i], gotDad, wantDad[i])
		}
		gotMom, ok := p.ARGStrand(7, 1, i, len(loci))
		if !ok {
			t.Fatalf("locus %d: no maternal record", loci[i])
		}
		if gotMom != wantMom[i] {
			t.Errorf("locus %d maternal: got strand %d, want %d", loci[i], gotMom, wantMom[i])
		}
	}
}

// A locus set wider than one word must pack correctly — the index-space/position-space
// distinction is the easiest thing to get wrong here.
func TestARGCaptureAcrossWordBoundary(t *testing.T) {
	const n = 200
	loci := make([]int, n)
	for i := range loci {
		loci[i] = i * 3 // positions 0..597, sparse; indices 0..199, dense
	}
	mask := make([]uint64, (600+63)/64)
	// Every third focal locus comes from copy 0.
	for i := 0; i < n; i += 3 {
		pos := loci[i]
		mask[pos/64] |= 1 << (uint(pos) % 64)
	}
	p := &Pop{}
	p.ARGCapture(1, loci, mask, make([]uint64, len(mask)))
	for i := 0; i < n; i++ {
		got, ok := p.ARGStrand(1, 0, i, n)
		if !ok {
			t.Fatalf("index %d: no record", i)
		}
		want := 1
		if i%3 == 0 {
			want = 0
		}
		if got != want {
			t.Fatalf("index %d (pos %d): got strand %d, want %d", i, loci[i], got, want)
		}
	}
	// The all-clear maternal mask means every maternal locus came from copy 1.
	for i := 0; i < n; i++ {
		if got, _ := p.ARGStrand(1, 1, i, n); got != 1 {
			t.Fatalf("index %d maternal: got strand %d, want 1", i, got)
		}
	}
}

// A missing record is reported as such, never as a default strand — the walker relies on
// this to terminate at founders instead of inventing an ancestor.
func TestARGStrandMissingRecord(t *testing.T) {
	p := &Pop{}
	if _, ok := p.ARGStrand(1, 0, 0, 4); ok {
		t.Error("empty ARGDB must report no record")
	}
	p.ARGCapture(1, []int{0, 1}, make([]uint64, 1), make([]uint64, 1))
	if _, ok := p.ARGStrand(2, 0, 0, 2); ok {
		t.Error("an unrecorded individual must report no record")
	}
	if _, ok := p.ARGStrand(1, 0, 5, 2); ok {
		t.Error("an out-of-range locus index must report no record")
	}
	if _, ok := p.ARGStrand(1, 2, 0, 2); ok {
		t.Error("an out-of-range gamete must report no record")
	}
}

// Capture is inert with no loci: no map is allocated on the hot path.
func TestARGCaptureNoLociIsInert(t *testing.T) {
	p := &Pop{}
	p.ARGCapture(1, nil, make([]uint64, 1), make([]uint64, 1))
	if p.ARGDB != nil {
		t.Error("ARGCapture allocated storage with no focal loci configured")
	}
}

// Provenance must be pruned by the pedigree cascade, not by a second refcount.
func TestARGPrunedWithPedigree(t *testing.T) {
	p := &Pop{}
	p.PedEnsure(1, 0, 0)
	p.PedEnsure(2, 0, 1)
	p.PedRecordBirth(3, 1, 2, 10, 0)
	p.ARGCapture(3, []int{0}, make([]uint64, 1), make([]uint64, 1))

	p.PedRecordDeath(1)
	p.PedRecordDeath(2)
	if _, ok := p.ARGDB[3]; !ok {
		t.Fatal("3 is alive; its provenance must be retained")
	}
	p.PedRecordDeath(3)
	if _, ok := p.ARGDB[3]; ok {
		t.Error("provenance outlived the pedigree node it indexes")
	}
	if len(p.Pedigree) != 0 {
		t.Errorf("pedigree not fully pruned: %d nodes", len(p.Pedigree))
	}
}
