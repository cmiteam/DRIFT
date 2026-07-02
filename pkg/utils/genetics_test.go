package utils

import (
	"drift/pkg/core"
	"testing"
)

// setBits builds a 2-word (128-bit) strand copy with the given positions set,
// using the canonical convention: position p -> word p/64, bit p%64.
func setBits(positions ...int) []uint64 {
	words := make([]uint64, 2)
	for _, p := range positions {
		words[p/64] |= 1 << uint(p%64)
	}
	return words
}

func newBlockModel(arms map[int]map[int][]int) *core.Model {
	return &core.Model{ChromosomeArms: arms}
}

func blockPop(copy0 []uint64) *core.Pop {
	return &core.Pop{
		Chromosomes: map[int][][]uint64{
			1: {copy0, {0, 0}},
		},
	}
}

// Regression: a single contiguous run spanning a 64-bit word boundary must
// count as ONE block. The old %064b string approach split it into two.
func TestCountContiguousBlocks_SpansWordBoundary(t *testing.T) {
	model := newBlockModel(map[int]map[int][]int{
		0: {0: {0, 128}, 1: {128, 0}},
	})
	// positions 60..67 are contiguous and straddle the word0/word1 boundary.
	pop := blockPop(setBits(60, 61, 62, 63, 64, 65, 66, 67))

	if got := CountContiguousBlocks(model, pop, 1, 0); got != 1 {
		t.Errorf("contiguous block across word boundary counted as %d, want 1", got)
	}
}

func TestCountContiguousBlocks_SeparateRuns(t *testing.T) {
	model := newBlockModel(map[int]map[int][]int{
		0: {0: {0, 128}, 1: {128, 0}},
	})
	// Two runs separated by a gap: {3,4,5} and {70,71}.
	pop := blockPop(setBits(3, 4, 5, 70, 71))

	if got := CountContiguousBlocks(model, pop, 1, 0); got != 2 {
		t.Errorf("two separated runs counted as %d, want 2", got)
	}
}

// Runs are counted per arm, so a block touching an arm boundary is not merged
// across it: a single arm fully set is one block, both arms set is two.
func TestCountContiguousBlocks_BoundedPerArm(t *testing.T) {
	model := newBlockModel(map[int]map[int][]int{
		0: {0: {0, 64}, 1: {64, 64}},
	})

	pArmOnly := blockPop(setBits(rangeBits(0, 64)...))
	if got := CountContiguousBlocks(model, pArmOnly, 1, 0); got != 1 {
		t.Errorf("full p arm counted as %d, want 1", got)
	}

	bothArms := blockPop(setBits(rangeBits(0, 128)...))
	if got := CountContiguousBlocks(model, bothArms, 1, 0); got != 2 {
		t.Errorf("both arms fully set counted as %d, want 2 (one per arm)", got)
	}
}

func TestCountContiguousBlocks_Empty(t *testing.T) {
	model := newBlockModel(map[int]map[int][]int{
		0: {0: {0, 64}, 1: {64, 64}},
	})
	if got := CountContiguousBlocks(model, blockPop([]uint64{0, 0}), 1, 0); got != 0 {
		t.Errorf("empty genome counted as %d, want 0", got)
	}
}

// rangeBits returns the positions [start, end).
func rangeBits(start, end int) []int {
	out := make([]int, 0, end-start)
	for p := start; p < end; p++ {
		out = append(out, p)
	}
	return out
}
