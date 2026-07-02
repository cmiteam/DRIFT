package analysis

// Identity-by-descent (IBD) extraction — Level 1 prototype.
//
// This is the "co-inherited founder-derived block" approach described in the
// roadmap (§7, cutting-edge bet 2: long-IBD as a signature of recent common
// ancestry). It runs against the EXISTING Pop structures with no changes to
// meiosis or the genome representation.
//
// WHAT A SET BIT MEANS HERE: seeding (see pkg/methods/seeding_*.go) sets the
// founder's chromosomes to all-1s, so every set bit marks a locus carrying
// material descended from the seed founder. Recombination in meiosis breaks
// those 1-runs into contiguous blocks. The bitwise AND of two individuals is
// therefore their co-inherited founder-derived material, and maximal runs of
// set bits in that AND are candidate shared (IBD) segments.
//
// HONEST LIMITATIONS (see roadmap for the Level 2 fix):
//   1. Because the founder's two chromosome copies are BOTH all-1s, descendants
//      cannot distinguish them. So this measures IBD-to-founder and over-coarsens
//      true pairwise IBD — the segment-length DISTRIBUTION has the right shape,
//      but absolute sharing is inflated. True IBD needs unique per-founder-
//      haplotype labels threaded through meiosis (Level 2 / tree sequences).
//   2. Only seed-derived material is visible (the genome is a descent tracer,
//      not a sequence). Fine for a single-couple model.
//   3. O(N^2) over sampled individuals × genome scan — sampling is mandatory.
//   4. Lengths are in BITS. Converting to centimorgans (what the IBD literature
//      uses) needs the recombination map from roadmap §6a.
//
// Wired into drift.go behind the track_IBD parameter (default off).

import (
	"drift/pkg/core"
	"drift/pkg/utils"
	"encoding/csv"
	"fmt"
	"os"
	"sort"
)

// IBDSegment is one maximal run of co-inherited founder-derived material shared
// between two individuals, bounded within a single chromosome arm.
type IBDSegment struct {
	Ind1, Ind2 int
	Chrom      int
	Arm        int // 0 = p arm, 1 = q arm
	StartBit   int // genome-global bit coordinate, inclusive
	EndBit     int // genome-global bit coordinate, exclusive
	LengthBits int
}

// IBDResult is the output of an extraction pass.
type IBDResult struct {
	Segments    []IBDSegment
	NumPairs    int            // number of individual pairs compared
	LengthHist  map[int]int    // log2-binned segment-length counts (the headline product)
	TotalByPair map[[2]int]int // total shared length per pair (basis for relatedness)
	MinLenBits  int
}

// ExtractIBD runs the pairwise pass over the given sampled individual IDs and
// returns all shared segments at least minLenBits long. Pass a SAMPLE of IDs —
// this is O(N^2). Individuals without tracked chromosomes are skipped.
func ExtractIBD(model *core.Model, pop *core.Pop, ids []int, minLenBits int) *IBDResult {
	if minLenBits < 1 {
		minLenBits = 1
	}
	res := &IBDResult{
		LengthHist:  map[int]int{},
		TotalByPair: map[[2]int]int{},
		MinLenBits:  minLenBits,
	}

	// Precompute the "carries founder material at locus p" mask (OR of the two
	// strand copies) once per individual.
	merged := make(map[int][]uint64, len(ids))
	for _, id := range ids {
		if m := orCopies(pop.Chromosomes[id]); m != nil {
			merged[id] = m
		}
	}

	for i := 0; i < len(ids); i++ {
		ca, ok := merged[ids[i]]
		if !ok {
			continue
		}
		for j := i + 1; j < len(ids); j++ {
			cb, ok := merged[ids[j]]
			if !ok {
				continue
			}
			res.NumPairs++
			shared := utils.BitwiseAND(ca, cb)

			// IBD runs must be bounded WITHIN an arm (adjacent arms/chromosomes
			// in the linear bit layout are not physically contiguous), mirroring
			// CountContiguousBlocks in pkg/utils/genetics.go.
			for chrom := range model.ChromosomeArms {
				for arm := 0; arm < 2; arm++ {
					armData := model.ChromosomeArms[chrom][arm]
					if len(armData) < 2 {
						continue
					}
					appendRuns(res, shared, armData[0], armData[1], minLenBits,
						ids[i], ids[j], chrom, arm)
				}
			}
		}
	}
	return res
}

// appendRuns scans bits [start, start+length) of the shared genome and records
// every maximal run of set bits whose length >= minLen as an IBDSegment.
//
// IMPORTANT: bits are read with the SAME convention meiosis writes them with —
// bit p lives at (word >> (p%64)) & 1, LSB = position 0. This deliberately does
// NOT use the fmt.Sprintf("%064b", ...) string method in genetics.go, which
// orders bits MSB-first within each word; the two disagree, and only this
// convention yields correct segment coordinates. (The string-based block
// counter has a latent inconsistency worth revisiting.)
func appendRuns(res *IBDResult, shared []uint64, start, length, minLen int,
	ind1, ind2, chrom, arm int) {

	end := start + length
	runStart := -1

	flush := func(runEnd int) {
		if runStart < 0 {
			return
		}
		if runEnd-runStart >= minLen {
			seg := IBDSegment{
				Ind1: ind1, Ind2: ind2, Chrom: chrom, Arm: arm,
				StartBit: runStart, EndBit: runEnd, LengthBits: runEnd - runStart,
			}
			res.Segments = append(res.Segments, seg)
			res.LengthHist[logBin(seg.LengthBits)]++
			res.TotalByPair[pairKey(ind1, ind2)] += seg.LengthBits
		}
		runStart = -1
	}

	for p := start; p < end; p++ {
		wordIdx := p / 64
		if wordIdx >= len(shared) {
			break
		}
		set := (shared[wordIdx]>>(uint(p)%64))&1 == 1
		if set {
			if runStart < 0 {
				runStart = p
			}
		} else {
			flush(p)
		}
	}
	flush(end) // close any run reaching the arm boundary
}

// orCopies returns the OR of an individual's two strand copies, or nil if the
// individual tracks no chromosomes.
func orCopies(c [][]uint64) []uint64 {
	if len(c) < 2 || c[0] == nil || c[1] == nil {
		return nil
	}
	return utils.BitwiseOR(c[0], c[1])
}

// pairKey returns a canonical (low, high) key for an unordered pair.
func pairKey(a, b int) [2]int {
	if a <= b {
		return [2]int{a, b}
	}
	return [2]int{b, a}
}

// logBin returns floor(log2(n)); bin k covers segment lengths [2^k, 2^(k+1)).
func logBin(n int) int {
	bin := 0
	for n > 1 {
		n >>= 1
		bin++
	}
	return bin
}

// SampleIDs returns up to n individual IDs that have tracked chromosomes,
// chosen at random. If n <= 0 or n >= the number of eligible individuals, all
// eligible IDs are returned.
//
// Uses the shared per-run RNG (utils), so sampling is reproducible under a fixed seed.
func SampleIDs(pop *core.Pop, n int) []int {
	ids := make([]int, 0, len(pop.Chromosomes))
	for id := range pop.Chromosomes {
		ids = append(ids, id)
	}
	if n <= 0 || n >= len(ids) {
		sort.Ints(ids) // stable output when not subsampling
		return ids
	}
	utils.RandShuffle(len(ids), func(i, j int) { ids[i], ids[j] = ids[j], ids[i] })
	ids = ids[:n]
	sort.Ints(ids)
	return ids
}

// SaveIBD writes two CSVs: the full segment table and the log-binned
// segment-length distribution (the headline scientific product). It also prints
// a one-line summary, matching the other analyses' end-of-run reporting.
func SaveIBD(model *core.Model, result *IBDResult) error {
	if result == nil {
		return nil
	}

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]

	// --- Segment table ---
	segPath := fmt.Sprintf("%s/%s_ibd_segments_run%d_year%d.csv",
		model.ResultsDir, model.ModelName, run, year)
	if err := writeSegments(segPath, result); err != nil {
		return err
	}

	// --- Length distribution ---
	histPath := fmt.Sprintf("%s/%s_ibd_lengthdist_run%d_year%d.csv",
		model.ResultsDir, model.ModelName, run, year)
	if err := writeLengthDist(histPath, result); err != nil {
		return err
	}

	totalLen := 0
	for _, seg := range result.Segments {
		totalLen += seg.LengthBits
	}
	meanLen := 0.0
	if len(result.Segments) > 0 {
		meanLen = float64(totalLen) / float64(len(result.Segments))
	}
	fmt.Printf("IBD: %d segments >= %d bits across %d pairs (mean length %.0f bits)\n",
		len(result.Segments), result.MinLenBits, result.NumPairs, meanLen)

	return nil
}

func writeSegments(path string, result *IBDResult) error {
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create IBD segment file: %v", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	if err := writer.Write([]string{
		"Ind1", "Ind2", "Chrom", "Arm", "StartBit", "EndBit", "LengthBits",
	}); err != nil {
		return err
	}
	for _, s := range result.Segments {
		row := []string{
			fmt.Sprintf("%d", s.Ind1),
			fmt.Sprintf("%d", s.Ind2),
			fmt.Sprintf("%d", s.Chrom),
			fmt.Sprintf("%d", s.Arm),
			fmt.Sprintf("%d", s.StartBit),
			fmt.Sprintf("%d", s.EndBit),
			fmt.Sprintf("%d", s.LengthBits),
		}
		if err := writer.Write(row); err != nil {
			return err
		}
	}
	return nil
}

func writeLengthDist(path string, result *IBDResult) error {
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create IBD length-distribution file: %v", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	if err := writer.Write([]string{
		"Log2Bin", "MinLengthBits", "MaxLengthBits", "SegmentCount",
	}); err != nil {
		return err
	}

	bins := make([]int, 0, len(result.LengthHist))
	for bin := range result.LengthHist {
		bins = append(bins, bin)
	}
	sort.Ints(bins)

	for _, bin := range bins {
		low := 1 << uint(bin)
		high := (1 << uint(bin+1)) - 1
		row := []string{
			fmt.Sprintf("%d", bin),
			fmt.Sprintf("%d", low),
			fmt.Sprintf("%d", high),
			fmt.Sprintf("%d", result.LengthHist[bin]),
		}
		if err := writer.Write(row); err != nil {
			return err
		}
	}
	return nil
}
