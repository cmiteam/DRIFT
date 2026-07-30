package analysis

// VCF import — roadmap §6a "Real-data interoperability". The inverse of ExportVCF
// (vcf.go): parse a real 1000 Genomes / HGDP (or DRIFT-exported) VCF into DRIFT's
// structures so the SAME analysis pipeline (Fst / f-stats / joint-SFS / single-pop
// SFS / §6g EHH-iHS-XP-EHH / §6h neutral validation / §6e dating) runs on simulated
// and on real data. A critic can then be answered in the field's own statistics,
// computed by the same code, on their own panel.
//
// WHY BOTH STRUCTURES. An empirical probe (documented in the §6a project notes)
// established that the pipelines split across DRIFT's two variant substrates:
//   * the founder BITFIELD (pop.Chromosomes) feeds Fst / f-stats / joint-SFS /
//     single-pop SFS / temporal-Ne / LD / het-distance (all via buildContigs +
//     demeMembers/countSite/bitAt);
//   * the de-novo POOL (pop.IndMutations + pop.MutationPool) feeds §6g haplotype
//     stats, §6h neutral validation, §6e dating, §6f load.
// A panel loaded into only one is invisible to the readers of the other. So an
// import that must serve ALL the named pipelines populates BOTH, CONSISTENTLY: the
// same VCF site becomes one genome bit AND one pool Mutation at that same Position
// (Effect 0 — a genotype panel carries no selection coefficient), and the carrying
// strands are recorded in the bitfield word and the IndMutations strand list in
// lockstep. This is the same bitfield<->pool reconciliation linkage_model="arm"
// reaches for in a live run (§1); here it is imposed directly at load.
//
// COORDINATES. Real bp coordinates cannot survive: DRIFT's genome is a few thousand
// bits, a real chromosome is hundreds of millions of bp with millions of variants.
// So import REMAPS each retained biallelic-SNP site to one genome bit, in file order
// grouped by contig, and synthesizes a ChromosomeArms layout (one arm per contig,
// [runningBitOffset, siteCount]) plus genome_bits = total sites. POS on a subsequent
// ExportVCF is therefore the per-contig site INDEX, not the original bp — physical
// distance is intentionally discarded (a cM map from bp is a documented follow-up).
//
// POLARIZATION. DRIFT bit 0 = ancestral, bit 1 = derived; a VCF's REF/ALT is NOT
// ancestral/derived. Policy "aa" (default) reads the ancestral allele from the AA
// INFO field (1000G/HGDP carry it) when it is a known base, and flips the derived
// code to REF when AA==ALT; it falls back to REF=ancestral (derived=ALT) when AA is
// absent/unknown. Policy "ref" always assumes REF=ancestral. Mis-polarization folds
// the SFS and flips iHS, so the AA path matters for those stats.
//
// SITE FILTERING. Strict biallelic SNPs only: multiallelic (ALT has a comma) and
// indel/MNP (REF or ALT not a single base) sites are skipped, matching DRIFT's
// biallelic single-bit model. Genotypes are phased (the two strand copies are
// tracked independently) — 1000G/HGDP are phased; an unphased "/" separator is
// accepted in column order and counted. Missing genotypes ("." ) read as ancestral;
// a site whose missing fraction exceeds MaxMissingFrac is dropped. Monomorphic sites
// (no variation in the panel) are dropped — they carry no information and would only
// inflate genome_bits.

import (
	"bufio"
	"compress/gzip"
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"fmt"
	"os"
	"sort"
	"strings"
)

// ImportPolicy configures how a VCF is mapped onto DRIFT's biallelic bit model.
type ImportPolicy struct {
	// Polarize selects ancestral/derived assignment: "aa" (default) reads the AA
	// INFO field and falls back to REF=ancestral; "ref" always assumes REF=ancestral.
	Polarize string
	// MaxMissingFrac drops any site whose fraction of missing genotype calls exceeds
	// this (over 2*NumSamples strands). Default 1.0 keeps every site.
	MaxMissingFrac float64
	// PopMap maps a sample name to a deme id (island-model subpopulation). Nil ⇒
	// every sample is deme 0 (single-population; Fst/f-stats then report no structure).
	PopMap map[string]int
	// MaxSites caps the number of retained sites (0 = no cap). Useful for tests and
	// for bounding a large panel; pre-filter with bcftools for real work.
	MaxSites int
}

// DefaultImportPolicy returns Rob's chosen defaults: AA polarization with REF
// fallback, keep all sites, single deme.
func DefaultImportPolicy() ImportPolicy {
	return ImportPolicy{Polarize: "aa", MaxMissingFrac: 1.0}
}

// ImportStats summarizes an import pass (returned for logging / tests).
type ImportStats struct {
	Path        string
	NumSamples  int
	NumSites    int // biallelic SNP sites kept = genome_bits
	GenomeBits  int
	NumContigs  int
	NumDemes    int
	ContigNames []string

	SkipMulti     int // multiallelic sites skipped
	SkipIndel     int // indel/MNP sites skipped
	SkipMissing   int // sites dropped for excess missingness
	SkipMonomorph int // sites with no variation in the panel
	PolarizedAA   int // sites polarized from a known AA INFO value
	PolarizedFlip int // of those, sites where AA==ALT (derived flipped to REF)
	PolarizedFall int // sites falling back to REF=ancestral (AA absent/unknown)
	UnphasedCalls int // genotype calls seen with a "/" (unphased) separator
	MissingCalls  int // missing genotype calls read as ancestral
}

// importSite is a parsed, retained site held until bit positions are assigned.
type importSite struct {
	contig  string
	derived []bool // length 2*NumSamples; strand s of sample i at index 2*i+s
}

// ImportVCF parses the VCF at path and populates model + pop with the panel in BOTH
// the founder bitfield and the de-novo mutation pool, consistently (see file doc).
// It OVERWRITES model.ChromosomeArms / model.FreeParameters["genome_bits"] and
// model.Parameters["num_demes"] to match the imported panel, and fills pop.IndData /
// Chromosomes / IndMutations / MutationPool. A path ending in ".gz" is decompressed.
func ImportVCF(path string, model *core.Model, pop *core.Pop, pol ImportPolicy) (*ImportStats, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open VCF: %w", err)
	}
	defer f.Close()

	var r *bufio.Reader
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("gzip: %w", err)
		}
		defer gz.Close()
		r = bufio.NewReader(gz)
	} else {
		r = bufio.NewReader(f)
	}

	if pol.Polarize == "" {
		pol.Polarize = "aa"
	}
	if pol.MaxMissingFrac <= 0 {
		pol.MaxMissingFrac = 1.0
	}

	stats := &ImportStats{Path: path}
	var sampleNames []string
	var sites []importSite
	contigOrder := []string{} // first-seen contig order
	contigSeen := map[string]bool{}

	// bufio.Reader.ReadString handles arbitrarily long variant lines (a panel with
	// thousands of samples has very wide rows).
	for {
		line, readErr := r.ReadString('\n')
		line = strings.TrimRight(line, "\r\n")
		if line != "" {
			if strings.HasPrefix(line, "##") {
				// meta line; ignore (we synthesize our own contigs)
			} else if strings.HasPrefix(line, "#CHROM") {
				cols := strings.Split(line, "\t")
				if len(cols) > 9 {
					sampleNames = append(sampleNames, cols[9:]...)
				}
			} else if len(sampleNames) > 0 {
				site, kind := parseVariantLine(line, len(sampleNames), pol, stats)
				switch kind {
				case siteKept:
					if pol.MaxSites == 0 || len(sites) < pol.MaxSites {
						sites = append(sites, site)
						if !contigSeen[site.contig] {
							contigSeen[site.contig] = true
							contigOrder = append(contigOrder, site.contig)
						}
					}
				}
			}
		}
		if readErr != nil {
			break
		}
	}

	if len(sampleNames) == 0 {
		return nil, fmt.Errorf("no #CHROM header / no samples in %s", path)
	}
	stats.NumSamples = len(sampleNames)
	stats.NumSites = len(sites)
	stats.GenomeBits = len(sites)
	stats.NumContigs = len(contigOrder)
	stats.ContigNames = contigOrder
	if len(sites) == 0 {
		return stats, fmt.Errorf("no biallelic SNP sites retained from %s", path)
	}

	// --- Assign bit positions: group by contig in first-seen order, file order
	// within a contig. Build a synthetic ChromosomeArms (one arm per contig). ---
	sort.SliceStable(sites, func(i, j int) bool {
		ci := indexOf(contigOrder, sites[i].contig)
		cj := indexOf(contigOrder, sites[j].contig)
		return ci < cj // stable: preserves file order within a contig
	})

	arms := map[int]map[int][]int{}
	contigStart := map[string]int{}
	contigCount := map[string]int{}
	for _, s := range sites {
		contigCount[s.contig]++
	}
	running := 0
	for chromIdx, name := range contigOrder {
		contigStart[name] = running
		arms[chromIdx+1] = map[int][]int{0: {running, contigCount[name]}}
		running += contigCount[name]
	}
	genomeBits := running

	// --- Deme assignment: distinct populations → 0..k-1 (sorted for determinism). ---
	demeOf := make([]int, len(sampleNames))
	numDemes := 1
	if pol.PopMap != nil {
		raw := map[int]bool{}
		for _, name := range sampleNames {
			if d, ok := pol.PopMap[name]; ok {
				raw[d] = true
			}
		}
		ids := make([]int, 0, len(raw))
		for d := range raw {
			ids = append(ids, d)
		}
		sort.Ints(ids)
		remap := map[int]int{}
		for i, d := range ids {
			remap[d] = i
		}
		for i, name := range sampleNames {
			demeOf[i] = remap[pol.PopMap[name]] // absent ⇒ 0 (map zero-value key -> remap[0]); default deme 0
		}
		if len(ids) > 0 {
			numDemes = len(ids)
		}
	}
	stats.NumDemes = numDemes

	// --- Populate model ---
	if model.FreeParameters == nil {
		model.FreeParameters = map[string]int{}
	}
	if model.Parameters == nil {
		model.Parameters = map[string]float64{}
	}
	model.ChromosomeArms = arms
	model.ArmCM = nil      // imported panels carry no cM map (physical distance discarded)
	model.GeneticMap = nil // rebuilt lazily if a recombination model needs it
	model.FreeParameters["genome_bits"] = genomeBits
	model.FreeParameters["mutID"] = genomeBits // pool mutation ids are 1..genomeBits
	model.Parameters["num_demes"] = float64(numDemes)

	// --- Populate pop: bitfield + pool, consistently ---
	words := (genomeBits + 63) / 64
	if pop.IndData == nil {
		pop.IndData = map[int][]int{}
	}
	if pop.Chromosomes == nil {
		pop.Chromosomes = map[int][][]uint64{}
	}
	if pop.IndMutations == nil {
		pop.IndMutations = map[int]map[int][]int{}
	}
	if pop.MutationPool == nil {
		pop.MutationPool = map[int]core.Mutation{}
	}
	if pop.MutationHist == nil {
		pop.MutationHist = map[int]int{}
	}
	for i := range sampleNames {
		row := individual.MakeIndData()
		row[individual.Deme] = demeOf[i]
		row[individual.BirthYear] = 0
		pop.IndData[i] = row
		pop.Chromosomes[i] = [][]uint64{make([]uint64, words), make([]uint64, words)}
		pop.IndMutations[i] = map[int][]int{0: {}, 1: {}}
	}
	// One mutation lineage per site (id = bit+1, Position = bit). Walk sites in the
	// assigned (contig-grouped) order so bit == global site index.
	for bit, s := range sites {
		mid := bit + 1
		count := 0
		for i := range sampleNames {
			for strand := 0; strand < 2; strand++ {
				if s.derived[2*i+strand] {
					pop.Chromosomes[i][strand][bit/64] |= 1 << uint(bit%64)
					pop.IndMutations[i][strand] = append(pop.IndMutations[i][strand], mid)
					count++
				}
			}
		}
		pop.MutationPool[mid] = core.Mutation{
			Id: mid, Position: bit, Effect: 0, Origin: -1, Count: count, Class: 0, Size: 1,
		}
	}
	// Allele counts (mirrors the seeders) so downstream reporting is consistent.
	for i := range sampleNames {
		pop.IndData[i][individual.AlleleCount] =
			utils.CountSetBits(pop.Chromosomes[i][0]) + utils.CountSetBits(pop.Chromosomes[i][1])
	}

	return stats, nil
}

type siteKind int

const (
	siteKept siteKind = iota
	siteSkipMulti
	siteSkipIndel
	siteSkipMissing
	siteSkipMono
)

// parseVariantLine parses one VCF data row into an importSite, applying the biallelic
// filter and polarization. It updates the running counters in stats and returns the
// disposition. The returned importSite is valid only when kind == siteKept.
func parseVariantLine(line string, numSamples int, pol ImportPolicy, stats *ImportStats) (importSite, siteKind) {
	cols := strings.Split(line, "\t")
	if len(cols) < 10 {
		return importSite{}, siteSkipIndel // malformed: treat as unusable
	}
	chrom := cols[0]
	ref := cols[3]
	alt := cols[4]
	info := cols[7]
	format := cols[8]

	if strings.Contains(alt, ",") {
		stats.SkipMulti++
		return importSite{}, siteSkipMulti
	}
	if len(ref) != 1 || len(alt) != 1 {
		stats.SkipIndel++
		return importSite{}, siteSkipIndel
	}

	// Polarization: derivedCode is the allele code (0=REF, 1=ALT) that is DERIVED.
	derivedCode := 1 // default REF=ancestral ⇒ ALT is derived
	if pol.Polarize == "aa" {
		if aa, ok := ancestralAllele(info); ok {
			stats.PolarizedAA++
			switch aa {
			case strings.ToUpper(ref):
				derivedCode = 1
			case strings.ToUpper(alt):
				derivedCode = 0
				stats.PolarizedFlip++
			default:
				derivedCode = 1
				stats.PolarizedFall++
				stats.PolarizedAA-- // AA present but neither REF nor ALT ⇒ not usable
			}
		} else {
			stats.PolarizedFall++
		}
	} else {
		stats.PolarizedFall++
	}

	gtIdx := formatGTIndex(format)

	site := importSite{contig: chrom, derived: make([]bool, 2*numSamples)}
	missing := 0
	derivedCount := 0
	for i := 0; i < numSamples; i++ {
		col := 9 + i
		if col >= len(cols) {
			// short row: treat remaining calls as missing
			missing += 2
			continue
		}
		a0, a1, miss, unphased := parseGenotype(cols[col], gtIdx)
		if unphased {
			stats.UnphasedCalls++
		}
		if miss {
			missing += 2
		}
		if a0 == derivedCode {
			site.derived[2*i] = true
			derivedCount++
		}
		if a1 == derivedCode {
			site.derived[2*i+1] = true
			derivedCount++
		}
	}
	stats.MissingCalls += missing

	total := 2 * numSamples
	if total > 0 && float64(missing)/float64(total) > pol.MaxMissingFrac {
		stats.SkipMissing++
		return importSite{}, siteSkipMissing
	}
	if derivedCount == 0 || derivedCount == total {
		stats.SkipMonomorph++
		return importSite{}, siteSkipMono
	}
	return site, siteKept
}

// ancestralAllele extracts a known ancestral base (A/C/G/T, uppercased) from the AA
// INFO field. 1000G writes e.g. "AA=g|||" (base before the first '|', lower-case =
// low confidence); "." / "-" / "N" / absent ⇒ unknown (ok=false).
func ancestralAllele(info string) (string, bool) {
	for _, kv := range strings.Split(info, ";") {
		if strings.HasPrefix(kv, "AA=") {
			v := kv[3:]
			if i := strings.IndexByte(v, '|'); i >= 0 {
				v = v[:i]
			}
			v = strings.ToUpper(strings.TrimSpace(v))
			switch v {
			case "A", "C", "G", "T":
				return v, true
			}
			return "", false
		}
	}
	return "", false
}

// formatGTIndex returns the position of the GT sub-field in a FORMAT string (0 when
// absent, since GT is conventionally first).
func formatGTIndex(format string) int {
	for i, k := range strings.Split(format, ":") {
		if k == "GT" {
			return i
		}
	}
	return 0
}

// parseGenotype parses a sample's GT into two allele codes. Returns (a0, a1, missing,
// unphased). A code is 0 (REF), 1 (ALT), or -1 (missing/other, read as ancestral by
// the caller). A haploid call (single allele) is treated as homozygous.
func parseGenotype(field string, gtIdx int) (int, int, bool, bool) {
	sub := field
	if gtIdx > 0 || strings.IndexByte(field, ':') >= 0 {
		parts := strings.Split(field, ":")
		if gtIdx < len(parts) {
			sub = parts[gtIdx]
		} else {
			sub = parts[0]
		}
	}
	unphased := strings.IndexByte(sub, '/') >= 0
	sep := "|"
	if unphased {
		sep = "/"
	}
	alleles := strings.Split(sub, sep)
	a0 := alleleCode(alleles[0])
	a1 := a0
	if len(alleles) > 1 {
		a1 = alleleCode(alleles[len(alleles)-1])
	}
	miss := a0 < 0 && a1 < 0
	return a0, a1, miss, unphased
}

// alleleCode maps a VCF allele token to 0 (REF), 1 (ALT), or -1 (missing/other).
func alleleCode(a string) int {
	switch a {
	case "0":
		return 0
	case "1":
		return 1
	default:
		return -1
	}
}

// LoadPanelFile reads a sample→population panel (e.g. the 1000 Genomes
// integrated_call_samples...panel) and returns a sample-name→deme-id map suitable
// for ImportPolicy.PopMap, plus the population labels in deme-id order. The file is
// whitespace/tab-separated with the sample id in column 0 and the population in
// column `popCol` (1-based data columns; the 1000G panel uses col 1 = "pop", col 2 =
// "super_pop"). A header row whose first token is "sample" is skipped. Populations
// are assigned deme ids 0..k-1 in first-seen order.
func LoadPanelFile(path string, popCol int) (map[string]int, []string, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, nil, fmt.Errorf("open panel file: %w", err)
	}
	if popCol < 1 {
		popCol = 1
	}
	popID := map[string]int{}
	var popNames []string
	out := map[string]int{}
	for _, line := range strings.Split(string(b), "\n") {
		line = strings.TrimRight(line, "\r")
		if strings.TrimSpace(line) == "" {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) <= popCol {
			continue
		}
		if fields[0] == "sample" {
			continue // header
		}
		pop := fields[popCol]
		if _, ok := popID[pop]; !ok {
			popID[pop] = len(popNames)
			popNames = append(popNames, pop)
		}
		out[fields[0]] = popID[pop]
	}
	if len(out) == 0 {
		return nil, nil, fmt.Errorf("no sample→population rows parsed from %s", path)
	}
	return out, popNames, nil
}

func indexOf(s []string, v string) int {
	for i, x := range s {
		if x == v {
			return i
		}
	}
	return -1
}
