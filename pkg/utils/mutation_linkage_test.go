package utils

import (
	"reflect"
	"testing"

	"drift/pkg/core"
)

// tiledArmsModel builds a model whose ChromosomeArms tile [0, genomeBits)
// contiguously, mirroring the real chromosome_data.csv layout (p arm then q arm,
// no gaps). Two chromosomes: chr1 p[0,10) q[10,25); chr2 p[25,35) q[35,50).
func tiledArmsModel(linkage string) *core.Model {
	m := &core.Model{
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"genome_bits": 50},
		ChromosomeArms: map[int]map[int][]int{
			1: {0: {0, 10}, 1: {10, 15}},
			2: {0: {25, 10}, 1: {35, 15}},
		},
	}
	if linkage != "" {
		m.StringParams["linkage_model"] = linkage
	}
	return m
}

// TestPositionArm checks that a flat genome bit index maps back to the unique
// (chrom, arm) whose tiled range contains it, and that out-of-range positions are
// reported as unmapped.
func TestPositionArm(t *testing.T) {
	m := tiledArmsModel("")
	cases := []struct {
		pos       int
		chrom     int
		arm       int
		inRange   bool
	}{
		{0, 1, 0, true},    // chr1 p arm start
		{9, 1, 0, true},    // chr1 p arm end (inclusive)
		{10, 1, 1, true},   // chr1 q arm start (boundary belongs to q)
		{24, 1, 1, true},   // chr1 q arm end
		{25, 2, 0, true},   // chr2 p arm start
		{49, 2, 1, true},   // chr2 q arm end (genome_bits-1)
		{50, 0, 0, false},  // one past the genome
		{-1, 0, 0, false},  // negative
	}
	for _, tc := range cases {
		chrom, arm, ok := PositionArm(m, tc.pos)
		if ok != tc.inRange || (ok && (chrom != tc.chrom || arm != tc.arm)) {
			t.Errorf("PositionArm(%d) = (%d,%d,%v), want (%d,%d,%v)",
				tc.pos, chrom, arm, ok, tc.chrom, tc.arm, tc.inRange)
		}
	}
}

// inheritOracleLegacy is a verbatim copy of the pre-linkage-model InheritMutations
// body (inverted polarity: strand0 kept where mask==0, strand1 where mask==1). The
// default linkage_model="legacy" path must reproduce it exactly (byte-identity).
func inheritOracleLegacy(pop *core.Pop, genomemask []uint64, parent, child, copy int) {
	if _, exists := pop.IndMutations[child]; !exists {
		pop.IndMutations[child] = map[int][]int{0: {}, 1: {}}
	}
	for _, mutationID := range pop.IndMutations[parent][0] {
		mutation := pop.MutationPool[mutationID]
		wordIdx := mutation.Position / 64
		bitIdx := mutation.Position % 64
		inheritedStrand := (genomemask[wordIdx] >> bitIdx) & 1
		if int(inheritedStrand) == 0 {
			pop.IndMutations[child][copy] = append(pop.IndMutations[child][copy], mutationID)
			pop.MutationCount++
			pop.MutationPool[mutationID] = mutation
		}
	}
	for _, mutationID := range pop.IndMutations[parent][1] {
		mutation := pop.MutationPool[mutationID]
		wordIdx := mutation.Position / 64
		bitIdx := mutation.Position % 64
		inheritedStrand := (genomemask[wordIdx] >> bitIdx) & 1
		if int(inheritedStrand) == 1 {
			pop.IndMutations[child][copy] = append(pop.IndMutations[child][copy], mutationID)
			mutation.Count++
			pop.MutationPool[mutationID] = mutation
		}
	}
}

// inheritPop builds a parent (id 1) carrying the given strand-0 / strand-1
// mutation ids, with a pool placing each id at position id (so id doubles as its
// genome position for these small fixtures).
func inheritPop(strand0, strand1 []int) *core.Pop {
	pool := map[int]core.Mutation{}
	for _, id := range append(append([]int{}, strand0...), strand1...) {
		pool[id] = core.Mutation{Id: id, Position: id, Count: 1}
	}
	return &core.Pop{
		IndMutations: map[int]map[int][]int{1: {0: strand0, 1: strand1}},
		MutationPool: pool,
	}
}

// TestInheritMutationsLegacyByteIdentical asserts the default path reproduces the
// original inverted-polarity logic field-for-field over a range of masks.
func TestInheritMutationsLegacyByteIdentical(t *testing.T) {
	strand0 := []int{3, 40, 100}
	strand1 := []int{5, 41, 101}
	// A few masks with bits set at various positions (need >=2 words for pos 100+).
	masks := [][]uint64{
		{0, 0},                               // all zero
		{^uint64(0), ^uint64(0)},             // all one
		{1 << 3, 1 << (100 % 64)},            // bit 3 and bit 100
		{1<<5 | 1<<40, 1<<(41%64) | 1<<(101%64)},
	}
	for mi, mask := range masks {
		want := inheritPop(strand0, strand1)
		inheritOracleLegacy(want, mask, 1, 2, 0)

		got := inheritPop(strand0, strand1)
		InheritMutations(tiledArmsModel("legacy"), got, mask, 1, 2, 0)
		// also verify the unset default (no linkage_model key) behaves as legacy
		gotDefault := inheritPop(strand0, strand1)
		InheritMutations(tiledArmsModel(""), gotDefault, mask, 1, 2, 0)

		if !reflect.DeepEqual(got.IndMutations, want.IndMutations) {
			t.Errorf("mask %d: IndMutations = %v, want %v", mi, got.IndMutations, want.IndMutations)
		}
		if !reflect.DeepEqual(gotDefault.IndMutations, want.IndMutations) {
			t.Errorf("mask %d: default(unset) IndMutations = %v, want %v", mi, gotDefault.IndMutations, want.IndMutations)
		}
		if got.MutationCount != want.MutationCount {
			t.Errorf("mask %d: MutationCount = %d, want %d", mi, got.MutationCount, want.MutationCount)
		}
		if !reflect.DeepEqual(got.MutationPool, want.MutationPool) {
			t.Errorf("mask %d: MutationPool = %v, want %v", mi, got.MutationPool, want.MutationPool)
		}
	}
}

// TestInheritMutationsArmPolarity asserts that linkage_model="arm" flips the mask
// polarity to match meiosis: a strand-0 mutation is inherited where mask==1 (the
// same rule meiosis uses to take the copy-0 founder allele), the OPPOSITE of the
// legacy path — so the pool now co-segregates with the founder bitfield.
func TestInheritMutationsArmPolarity(t *testing.T) {
	// Mutation id 3 on strand0 at position 3; id 5 on strand1 at position 5.
	// Mask sets bit 3 (==1) and leaves bit 5 (==0).
	mask := []uint64{1 << 3}

	// Legacy: strand0 kept where mask==0 -> bit3 set -> id3 NOT inherited.
	//         strand1 kept where mask==1 -> bit5 clear -> id5 NOT inherited.
	legacy := inheritPop([]int{3}, []int{5})
	InheritMutations(tiledArmsModel("legacy"), legacy, mask, 1, 2, 0)
	if got := legacy.IndMutations[2][0]; len(got) != 0 {
		t.Errorf("legacy inherited %v, want none", got)
	}

	// Arm: strand0 kept where mask==1 -> bit3 set -> id3 inherited.
	//      strand1 kept where mask==0 -> bit5 clear -> id5 inherited.
	arm := inheritPop([]int{3}, []int{5})
	InheritMutations(tiledArmsModel("arm"), arm, mask, 1, 2, 0)
	gotArm := arm.IndMutations[2][0]
	if len(gotArm) != 2 || !((gotArm[0] == 3 && gotArm[1] == 5) || (gotArm[0] == 5 && gotArm[1] == 3)) {
		t.Errorf("arm inherited %v, want both id 3 and id 5", gotArm)
	}
}
