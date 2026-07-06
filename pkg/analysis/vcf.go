package analysis

// VCF export — roadmap §6a "Real-data interoperability".
//
// Writes the sampled population's genomes as a standard VCF v4.2 file so the
// SAME external pipeline (PLINK, vcftools, ADMIXTURE, EIGENSOFT, scikit-allel)
// can run on DRIFT output and on real 1000 Genomes / HGDP data. The whole point
// is that a critic can't dismiss a result produced with their own tools on a
// standard interchange format. This is the highest-leverage single feature in
// §6a.
//
// WHAT A SITE / ALLELE MEANS HERE. DRIFT's genome is a bit array that traces
// descent, not a nucleotide sequence: bit = 0 is the ancestral state, bit = 1 is
// a founder-derived (or, with mutation tracking, mutation-carrying) allele. Each
// bit position is therefore emitted as one biallelic SNP with a fixed
// REF=A (ancestral) / ALT=T (derived) convention — the letters are placeholders;
// only the 0/1 dosage carries information. Downstream tools operate on the
// genotype dosages, so this is sufficient for PCA, ADMIXTURE, Fst, f-stats, etc.
//
// PHASING. The two strand copies (chromosomes[id][0] and [1]) are tracked
// independently through meiosis, so the genotypes are genuinely phased — emitted
// with the '|' separator, which unlocks haplotype-based tools directly.
//
// COORDINATES. Global bit positions are mapped to per-chromosome 1-based
// coordinates using model.ChromosomeArms: each chromosome's base is the minimum
// arm start, so POS = bit - base + 1 and each contig's length is its arm span.
//
// By default only sites segregating WITHIN THE SAMPLE are written (monomorphic
// sites carry no information and would bloat the file); set vcf_include_fixed to
// emit every position. Sampling (vcf_sample_size) bounds the column count.
//
// Wired into drift.go behind the export_VCF parameter (default off); requires
// track_DNA so chromosomes exist.

import (
	"bufio"
	"drift/pkg/core"
	"fmt"
	"os"
	"sort"
)

// vcfContig describes one chromosome's placement in the linear bit genome.
type vcfContig struct {
	chrom int
	base  int // global bit coordinate of this chromosome's first arm (POS origin)
	span  int // length in bits (max arm end - base); the ##contig length
	arms  [][2]int
}

// VCFStats summarizes an export pass (returned for logging / tests).
type VCFStats struct {
	NumSamples int
	NumSites   int // records written
	Path       string
}

// buildContigs turns model.ChromosomeArms into a chromosome-sorted list with the
// per-chromosome base offset and span needed for VCF coordinates.
func buildContigs(model *core.Model) []vcfContig {
	chroms := make([]int, 0, len(model.ChromosomeArms))
	for chrom := range model.ChromosomeArms {
		chroms = append(chroms, chrom)
	}
	sort.Ints(chroms)

	contigs := make([]vcfContig, 0, len(chroms))
	for _, chrom := range chroms {
		base, maxEnd := -1, 0
		arms := make([][2]int, 0, 2)
		for arm := 0; arm < 2; arm++ {
			a, ok := model.ChromosomeArms[chrom][arm]
			if !ok || len(a) < 2 {
				continue
			}
			start, length := a[0], a[1]
			arms = append(arms, [2]int{start, length})
			if base < 0 || start < base {
				base = start
			}
			if start+length > maxEnd {
				maxEnd = start + length
			}
		}
		if len(arms) == 0 {
			continue
		}
		// Emit arms in ascending genome order for deterministic POS order.
		sort.Slice(arms, func(i, j int) bool { return arms[i][0] < arms[j][0] })
		contigs = append(contigs, vcfContig{chrom: chrom, base: base, span: maxEnd - base, arms: arms})
	}
	return contigs
}

// bitAt reports whether bit p (LSB-first, matching meiosis and ibd.go) is set in
// the given strand words; out-of-range bits read as 0.
func bitAt(words []uint64, p int) bool {
	w := p / 64
	if w < 0 || w >= len(words) {
		return false
	}
	return (words[w]>>(uint(p)%64))&1 == 1
}

// ExportVCF writes the sampled individuals' genomes to a VCF v4.2 file in the
// model's results directory and returns a summary. Pass a SAMPLE of ids (see
// SampleIDs); individuals without tracked chromosomes are skipped. When
// includeFixed is false only sites polymorphic within the sample are written.
func ExportVCF(model *core.Model, pop *core.Pop, ids []int, includeFixed bool) (*VCFStats, error) {
	// Keep only individuals that actually carry two strand copies, in the given
	// (already sorted) id order, so the sample columns are stable.
	samples := make([]int, 0, len(ids))
	for _, id := range ids {
		if c, ok := pop.Chromosomes[id]; ok && len(c) >= 2 && c[0] != nil && c[1] != nil {
			samples = append(samples, id)
		}
	}

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	path := fmt.Sprintf("%s/%s_run%d_year%d.vcf", model.ResultsDir, model.ModelName, run, year)

	file, err := os.Create(path)
	if err != nil {
		return nil, fmt.Errorf("failed to create VCF file: %v", err)
	}
	defer file.Close()

	w := bufio.NewWriter(file)
	defer w.Flush()

	contigs := buildContigs(model)

	// --- Header ---
	fmt.Fprintln(w, "##fileformat=VCFv4.2")
	fmt.Fprintf(w, "##source=DRIFT model=%s run=%d year=%d\n", model.ModelName, run, year)
	fmt.Fprintln(w, "##FILTER=<ID=PASS,Description=\"All filters passed\">")
	fmt.Fprintln(w, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for the ALT allele\">")
	fmt.Fprintln(w, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">")
	fmt.Fprintln(w, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	fmt.Fprintln(w, "##ALT allele encodes the founder-derived (bit=1) state; REF the ancestral (bit=0) state; letters are placeholders.")
	for _, ct := range contigs {
		fmt.Fprintf(w, "##contig=<ID=%d,length=%d>\n", ct.chrom, ct.span)
	}

	// #CHROM header line (tab-separated) with one column per sample.
	fmt.Fprint(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
	for _, id := range samples {
		fmt.Fprintf(w, "\tIND%d", id)
	}
	fmt.Fprintln(w)

	an := 2 * len(samples)

	// --- Records ---
	// Reused genotype-string buffer to avoid per-site allocation.
	gt := make([]byte, 0, 3*len(samples)+3*len(samples))
	numSites := 0
	for _, ct := range contigs {
		for _, arm := range ct.arms {
			start, length := arm[0], arm[1]
			for bit := start; bit < start+length; bit++ {
				ac := 0
				gt = gt[:0]
				for _, id := range samples {
					c := pop.Chromosomes[id]
					a0, a1 := bitAt(c[0], bit), bitAt(c[1], bit)
					gt = append(gt, '\t')
					gt = appendAllele(gt, a0)
					gt = append(gt, '|')
					gt = appendAllele(gt, a1)
					if a0 {
						ac++
					}
					if a1 {
						ac++
					}
				}
				if !includeFixed && (ac == 0 || ac == an) {
					continue
				}
				pos := bit - ct.base + 1
				fmt.Fprintf(w, "%d\t%d\t.\tA\tT\t.\tPASS\tAC=%d;AN=%d\tGT",
					ct.chrom, pos, ac, an)
				w.Write(gt)
				fmt.Fprintln(w)
				numSites++
			}
		}
	}

	fmt.Printf("VCF: %d sites × %d samples → %s\n", numSites, len(samples), path)
	return &VCFStats{NumSamples: len(samples), NumSites: numSites, Path: path}, nil
}

func appendAllele(b []byte, set bool) []byte {
	if set {
		return append(b, '1')
	}
	return append(b, '0')
}
