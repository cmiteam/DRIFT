package analysis

// Runs of homozygosity (ROH) — roadmap §4 "Analysis, validation & reproducibility".
//
// A run of homozygosity is a contiguous stretch of an individual's diploid genome
// where the two inherited haplotypes are identical, i.e. the individual is
// homozygous at every position in the stretch. Long ROH are the classic signature
// of recent shared ancestry (consanguinity, a small-founder event, a bottleneck):
// a long tract survives only if it was inherited from a recent common ancestor
// before recombination had time to break it up. The genome-wide burden of ROH
// (F_ROH — the fraction of the autosome lying in runs) and its partition into
// short / medium / long classes is the standard, directly-1000G/HGDP-comparable
// read-out of inbreeding and demographic history.
//
// WHAT A SITE / HOMOZYGOTE MEANS HERE. DRIFT's genome is the founder-descent
// bitfield (pop.Chromosomes[id] = two uint64 strand copies). An individual is
// HOMOZYGOUS at genome bit p iff its two strand copies agree there
// (bitAt(c[0],p) == bitAt(c[1],p)) — whether both carry the ancestral (0) or both
// the derived (1) allele. A HETEROZYGOUS site (the two copies differ) breaks a run.
// This is exactly the substrate real ROH callers use (dense per-position diploid
// genotypes with heterozygotes as run-breakers), so the statistic transfers
// unchanged. Contrast the IBD module (ibd.go), which ORs the two strands together
// to measure BETWEEN-individual founder sharing; ROH compares an individual's two
// strands AGAINST EACH OTHER to measure WITHIN-individual auto-zygosity. The two
// are complementary and read the same bitfield with opposite strand operators.
//
// WHY THE FOUNDER BITFIELD, NOT THE MUTATION POOL. ROH is a dense-genotype concept
// (a contiguous tract of the physical genome). The de-novo MutationPool is sparse
// and carries no contiguous per-position genotype array, so it is the wrong
// substrate. The bitfield is dense, contiguous, and tiled by ChromosomeArms — the
// natural coordinate space, shared with the VCF export and Fst scans (buildContigs).
// Empirically (probe, seed_style=2 "max_het"): founders start heterozygous at every
// site (F_ROH = 0) and homozygosity accumulates via drift + recombination; total
// F_ROH rises toward a drift equilibrium while LONG-tract F_ROH peaks early then
// decays in an outbred population — and a bottleneck leaves a persistent excess of
// long ROH (≈4× at the long threshold) plus roughly halved heterozygosity long
// after demographic recovery. The short/medium/long partition captures exactly that
// split: total/short ≈ cumulative inbreeding, long ≈ recent bottleneck/consanguinity.
//
// UNITS. Run length is reported natively in genome BITS and, via the §6a genetic
// map (utils.EnsureGeneticMap, uniform recomb_cM_per_bit fallback), in centiMorgans.
// The short/medium/long class boundaries are expressed in cM (roh_short_max_cM /
// roh_long_min_cM) because "long" is physically a recombination-length notion — a
// long tract is one recombination has not yet broken, which is a cM statement, not a
// bit-count one — and cM is what makes the classes comparable to real ROH studies.
// A native-bit floor (roh_min_bits) filters single-site het-gap noise before
// classification, and works even with no configured map.
//
// RUN BOUNDARIES. Runs are bounded WITHIN A SINGLE ARM (chosen with Rob), matching
// the IBD module's per-arm bounding and DRIFT's legacy centromere-bracketing
// recombination — adjacent arms are not treated as one contiguous stretch.
//
// Wired into drift.go behind the track_ROH parameter (default off); requires
// track_DNA so the Chromosomes bitfield exists (individuals without two strand
// copies are skipped, so an empty bitfield yields an empty result). Read-only:
// with roh_sample_size = 0 it is a deterministic full-population scan that consumes
// no RNG, so a track_ROH run is byte-identical to a track_ROH-off run. Also runs on
// an imported real-data panel (runImportedAnalyses), where it measures ROH on actual
// 1000G/HGDP genotypes.

import (
	"drift/pkg/core"
	"drift/pkg/utils"
	"sort"
)

// ROH length classes. Boundaries are in centiMorgans (roh_short_max_cM,
// roh_long_min_cM); medium is the interval between them.
const (
	ROHShort = iota
	ROHMedium
	ROHLong
	numROHClass // = 3
)

// ROHRun is one maximal homozygous stretch within a single chromosome arm.
type ROHRun struct {
	ID         int
	Chrom, Arm int
	StartBit   int // global bit coordinate, inclusive
	EndBit     int // global bit coordinate, exclusive
	LengthBits int
	LengthCM   float64
	Class      int // ROHShort / ROHMedium / ROHLong (by cM)
}

// ROHIndividual is one sampled individual's ROH summary.
type ROHIndividual struct {
	ID        int
	HetSites  int     // heterozygous bits over the scanned genome (run-breakers)
	NumRuns   int     // called ROH (length >= roh_min_bits)
	TotalBits int     // SROH: total bits in called ROH
	TotalCM   float64 // SROH in cM
	FROHBits  float64 // TotalBits / scanned autosome bits
	FROHCM    float64 // TotalCM / total map cM
	// Per-class breakdown (index ROHShort/ROHMedium/ROHLong).
	ClassCount [numROHClass]int
	ClassBits  [numROHClass]int
	ClassCM    [numROHClass]float64
}

// ROHResult is the output of a ROH pass over a sample.
type ROHResult struct {
	Individuals []ROHIndividual // per sampled individual, sorted by ID
	NumSampled  int             // individuals with a full 2-strand bitfield
	ScannedBits int             // autosome bits scanned (sum of arm lengths)
	TotalMapCM  float64         // total genetic length scanned (sum of per-chrom cM)
	HasMap      bool            // a real cM map was present (vs the uniform fallback)

	// Class configuration echoed for the output header.
	MinBits    int
	ShortMaxCM float64
	LongMinCM  float64

	// Population aggregates (means over sampled individuals).
	MeanHetFrac       float64
	MeanNumRuns       float64
	MeanFROHBits      float64
	MeanFROHCM        float64
	MeanClassCount    [numROHClass]float64
	MeanClassFROHBits [numROHClass]float64 // mean per-class fraction of the scanned genome

	// Pooled run-length distribution (the headline product), log2-binned in bits
	// (reusing ibd.go's logBin), and the pooled per-class run totals.
	LengthHist  map[int]int
	ClassTotals [numROHClass]int
}

// ROHConfig carries the parsed ROH parameters.
type ROHConfig struct {
	MinBits    int     // native-bit floor to call a homozygous stretch an ROH
	ShortMaxCM float64 // upper bound of the short class (cM)
	LongMinCM  float64 // lower bound of the long class (cM)
}

// DefaultROHConfig returns the ROH configuration from the model parameters, with
// sane fallbacks so the pass is well-defined even on a minimally-configured model.
func DefaultROHConfig(model *core.Model) ROHConfig {
	cfg := ROHConfig{MinBits: 5, ShortMaxCM: 0.5, LongMinCM: 1.5}
	if v, ok := model.Parameters["roh_min_bits"]; ok && v > 0 {
		cfg.MinBits = int(v)
	}
	if v, ok := model.Parameters["roh_short_max_cM"]; ok && v > 0 {
		cfg.ShortMaxCM = v
	}
	if v, ok := model.Parameters["roh_long_min_cM"]; ok && v > 0 {
		cfg.LongMinCM = v
	}
	// Guard against an inverted configuration (long boundary below short).
	if cfg.LongMinCM < cfg.ShortMaxCM {
		cfg.LongMinCM = cfg.ShortMaxCM
	}
	return cfg
}

// classifyCM assigns a run to a length class by its cM length.
func (cfg ROHConfig) classifyCM(cm float64) int {
	switch {
	case cm < cfg.ShortMaxCM:
		return ROHShort
	case cm < cfg.LongMinCM:
		return ROHMedium
	default:
		return ROHLong
	}
}

// rohGeneticMap returns the model's genetic map, or builds a uniform one from
// recomb_cM_per_bit (default 0.05 cM/bit) when none is configured, so ROH always
// has a cM axis. The second return reports whether a real (non-fallback) map was
// present. Mirrors haploGeneticMap (haplostats.go); the locally built map is NOT
// stored on the model (no side effect on the recombination kernel).
func rohGeneticMap(model *core.Model) (*core.GeneticMap, bool) {
	if gm := utils.EnsureGeneticMap(model); gm != nil {
		return gm, true
	}
	rate := model.Parameters["recomb_cM_per_bit"]
	if rate <= 0 {
		rate = 0.05
	}
	return core.BuildGeneticMap(model.ChromosomeArms, model.ArmCM, rate), false
}

// ComputeROH runs the ROH pass over the given sampled individual IDs (see
// SampleIDs). Individuals without two tracked strand copies are skipped. Returns a
// result with no individuals (and zero aggregates) when nothing is eligible.
func ComputeROH(model *core.Model, pop *core.Pop, ids []int) *ROHResult {
	cfg := DefaultROHConfig(model)
	gm, hasMap := rohGeneticMap(model)

	// Deterministic chromosome order for scanning and total-cM accumulation.
	chroms := make([]int, 0, len(model.ChromosomeArms))
	for chrom := range model.ChromosomeArms {
		chroms = append(chroms, chrom)
	}
	sort.Ints(chroms)

	scannedBits, totalMapCM := 0, 0.0
	for _, chrom := range chroms {
		for arm := 0; arm < 2; arm++ {
			armData := model.ChromosomeArms[chrom][arm]
			if len(armData) >= 2 && armData[1] > 0 {
				scannedBits += armData[1]
			}
		}
		totalMapCM += gm.TotalCM(chrom)
	}

	res := &ROHResult{
		ScannedBits: scannedBits,
		TotalMapCM:  totalMapCM,
		HasMap:      hasMap,
		MinBits:     cfg.MinBits,
		ShortMaxCM:  cfg.ShortMaxCM,
		LongMinCM:   cfg.LongMinCM,
		LengthHist:  map[int]int{},
	}

	// Sorted, deduplicated id order so per-individual rows and the pooled histogram
	// are a deterministic function of the input.
	sorted := append([]int(nil), ids...)
	sort.Ints(sorted)

	var sumHetFrac, sumRuns, sumFROHBits, sumFROHCM float64
	var sumClassCount, sumClassFROHBits [numROHClass]float64

	for _, id := range sorted {
		c, ok := pop.Chromosomes[id]
		if !ok || len(c) < 2 || c[0] == nil || c[1] == nil {
			continue
		}
		ind := scanIndividual(id, c, model, chroms, gm, cfg, scannedBits, totalMapCM, res)
		res.Individuals = append(res.Individuals, ind)
		res.NumSampled++

		if scannedBits > 0 {
			sumHetFrac += float64(ind.HetSites) / float64(scannedBits)
		}
		sumRuns += float64(ind.NumRuns)
		sumFROHBits += ind.FROHBits
		sumFROHCM += ind.FROHCM
		for k := 0; k < numROHClass; k++ {
			sumClassCount[k] += float64(ind.ClassCount[k])
			if scannedBits > 0 {
				sumClassFROHBits[k] += float64(ind.ClassBits[k]) / float64(scannedBits)
			}
		}
	}

	if n := float64(res.NumSampled); n > 0 {
		res.MeanHetFrac = sumHetFrac / n
		res.MeanNumRuns = sumRuns / n
		res.MeanFROHBits = sumFROHBits / n
		res.MeanFROHCM = sumFROHCM / n
		for k := 0; k < numROHClass; k++ {
			res.MeanClassCount[k] = sumClassCount[k] / n
			res.MeanClassFROHBits[k] = sumClassFROHBits[k] / n
		}
	}
	return res
}

// scanIndividual scans one individual's bitfield arm-by-arm, records every called
// ROH (length >= cfg.MinBits) into res's pooled tallies, and returns the
// per-individual summary.
func scanIndividual(id int, c [][]uint64, model *core.Model, chroms []int,
	gm *core.GeneticMap, cfg ROHConfig, scannedBits int, totalMapCM float64, res *ROHResult) ROHIndividual {

	ind := ROHIndividual{ID: id}
	for _, chrom := range chroms {
		for arm := 0; arm < 2; arm++ {
			armData := model.ChromosomeArms[chrom][arm]
			if len(armData) < 2 || armData[1] <= 0 {
				continue
			}
			start, length := armData[0], armData[1]
			end := start + length

			runStart := -1
			flush := func(runEnd int) {
				if runStart < 0 {
					return
				}
				lenBits := runEnd - runStart
				if lenBits >= cfg.MinBits {
					cm := gm.CumCM(chrom, runEnd) - gm.CumCM(chrom, runStart)
					class := cfg.classifyCM(cm)
					ind.NumRuns++
					ind.TotalBits += lenBits
					ind.TotalCM += cm
					ind.ClassCount[class]++
					ind.ClassBits[class] += lenBits
					ind.ClassCM[class] += cm
					res.LengthHist[logBin(lenBits)]++
					res.ClassTotals[class]++
				}
				runStart = -1
			}

			for p := start; p < end; p++ {
				a0 := bitAt(c[0], p)
				a1 := bitAt(c[1], p)
				if a0 == a1 { // homozygous
					if runStart < 0 {
						runStart = p
					}
				} else { // heterozygous — breaks the run
					ind.HetSites++
					flush(p)
				}
			}
			flush(end) // close a run reaching the arm boundary
		}
	}

	if scannedBits > 0 {
		ind.FROHBits = float64(ind.TotalBits) / float64(scannedBits)
	}
	if totalMapCM > 0 {
		ind.FROHCM = ind.TotalCM / totalMapCM
	}
	return ind
}
