package utils

import (
	"drift/pkg/core"
	"fmt"
	"math/bits"
	"strings"
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

func CountContiguousBlocks(model *core.Model, pop *core.Pop, ind int, copy int) int {
	blockCount := 0
	genomestring := uint64ArrayToBitString(pop.Chromosomes[ind][copy])
	for chrom, _ := range model.ChromosomeArms {
		pstart := model.ChromosomeArms[chrom][0][0]
		plen := model.ChromosomeArms[chrom][0][1]
		if pstart+plen <= len(genomestring) {
			psegment := genomestring[pstart : pstart+plen]
			blockCount += countBlocksInRange(psegment)
		}

		qstart := model.ChromosomeArms[chrom][1][0]
		qlen := model.ChromosomeArms[chrom][1][1]
		if qstart+qlen <= len(genomestring) {
			qsegment := genomestring[qstart : qstart+qlen]
			blockCount += countBlocksInRange(qsegment)
		}
	}
	return blockCount
}

func countBlocksInRange(genomesegment string) int {
	blocks := strings.Split(genomesegment, "0")
	blockCount := 0
	for _, block := range blocks {
		if len(block) > 0 {
			blockCount++
		}
	}
	return blockCount
}

func uint64ArrayToBitString(genomesegment []uint64) string {
	var bitString strings.Builder
	for _, value := range genomesegment {
		bitString.WriteString(fmt.Sprintf("%064b", value))
	}
	return bitString.String()
}

func BitwiseOR(a, b []uint64) []uint64 {
	result := make([]uint64, len(a))
	for i := range a {
		result[i] = a[i] | b[i]
	}
	return result
}
