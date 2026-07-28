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

// PositionArm maps a genome bit position back to the chromosome and arm whose
// [start, start+length) range contains it (roadmap §1, linkage-aware positions).
// The chromosome arms tile the genome coordinate space contiguously, so a de-novo
// mutation's flat bit index already lives on a specific arm; this makes that
// mapping explicit so downstream analysis can partition mutations by chromosome/
// arm, and so linkage-aware inheritance reasons in the same coordinate frame as
// the founder bitfield. Arms never overlap, so the containing (chrom, arm) is
// unique. Returns (chrom, arm, true) on a hit, or (0, 0, false) when the position
// falls outside every arm (only possible for a malformed/partial map).
func PositionArm(model *core.Model, pos int) (int, int, bool) {
	for chrom := range model.ChromosomeArms {
		for arm := 0; arm < 2; arm++ {
			armData := model.ChromosomeArms[chrom][arm]
			if len(armData) < 2 {
				continue
			}
			if pos >= armData[0] && pos < armData[0]+armData[1] {
				return chrom, arm, true
			}
		}
	}
	return 0, 0, false
}

// EnsureGeneticMap returns the model's recombination map (roadmap §6a), building and
// caching it on first use. It builds only when a map is actually wanted — either the
// chromosome file carried a cM column (len(ArmCM) > 0) or recombination_model="map"
// requests one — and returns nil otherwise, which is the signal every caller uses to
// take the legacy path (createMask) or the ne_cM_per_bit scalar (recombFraction).
//
// When ArmCM is empty but a map is requested, the whole map is derived from a uniform
// recomb_cM_per_bit rate (default 1.0 cM/bit ≈ human average at ~1 bit per Mb); when
// ArmCM is present, missing arms fall back to that same rate. Pure (no RNG), built at
// most once, so it does not perturb the RNG stream and is safe to call in the sim loop.
func EnsureGeneticMap(model *core.Model) *core.GeneticMap {
	if model.GeneticMap != nil {
		return model.GeneticMap
	}
	hasCM := len(model.ArmCM) > 0
	mapMode := model.StringParam("recombination_model", "legacy") == "map"
	if !hasCM && !mapMode {
		return nil
	}
	rate, ok := model.Parameters["recomb_cM_per_bit"]
	if !ok {
		rate = 1.0
	}
	model.GeneticMap = core.BuildGeneticMap(model.ChromosomeArms, model.ArmCM, rate)
	return model.GeneticMap
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
