package utils

import (
	"drift/pkg/core"
	"math/bits"
)

func CountSetBits(words []uint64) int {
	count := 0
	for _, word := range words {
		count += bits.OnesCount64(word)
	}
	return count
}

func CountSetBitsSingleVar(n uint64) int {
	count := 0
	for n > 0 {
		count += int(n & 1)
		n >>= 1
	}
	return count
}

// CountContiguousBlocks counts maximal runs of set bits (haplotype blocks
// descended from the seed) within each chromosome arm of the given strand copy.
//
// Bits are read with the same convention they are SET with elsewhere
// (meiosis, createMask, the SFS counter, IBD extraction): genome position p
// lives at (words[p/64] >> (p%64)) & 1, LSB = position 0. Counting is done
// per arm so a run cannot be merged across an arm boundary, and runs are never
// split at 64-bit word boundaries (the previous %064b string approach reversed
// within-word bit order and scrambled adjacency at word boundaries, inflating
// the count for any block straddling a boundary).
func CountContiguousBlocks(model *core.Model, pop *core.Pop, ind int, copy int) int {
	words := pop.Chromosomes[ind][copy]
	blockCount := 0
	for chrom := range model.ChromosomeArms {
		for arm := 0; arm < 2; arm++ {
			armData := model.ChromosomeArms[chrom][arm]
			if len(armData) < 2 {
				continue
			}
			blockCount += countBlocksInBitRange(words, armData[0], armData[1])
		}
	}
	return blockCount
}

// countBlocksInBitRange counts maximal runs of set bits in positions
// [start, start+length) of the genome word array.
func countBlocksInBitRange(words []uint64, start, length int) int {
	blocks := 0
	inRun := false
	for p := start; p < start+length; p++ {
		wordIdx := p / 64
		if wordIdx >= len(words) {
			break
		}
		set := (words[wordIdx]>>(uint(p)%64))&1 == 1
		if set && !inRun {
			blocks++ // start of a new run
		}
		inRun = set
	}
	return blocks
}

func BitwiseOR(a, b []uint64) []uint64 {
	result := make([]uint64, len(a))
	for i := range a {
		result[i] = a[i] | b[i]
	}
	return result
}

// BitwiseAND returns the per-word AND of two genome bit arrays. Length is the
// shorter of the two inputs, so mismatched-length arrays are handled safely.
func BitwiseAND(a, b []uint64) []uint64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	result := make([]uint64, n)
	for i := 0; i < n; i++ {
		result[i] = a[i] & b[i]
	}
	return result
}
