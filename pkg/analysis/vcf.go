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
// THE DE-NOVO POOL (TMR4A.md W6) is a SECOND variant substrate this file wrote
// nothing about until vcf_include_mutations existed; see vcf_merged.go for what
// merging it in means, why it is one record per mutation rather than an OR onto
// the bit, and why linkage_model="arm" is required for the result to be coherent.
//
// Wired into drift.go behind the export_VCF parameter (default off); requires
// track_DNA so chromosomes exist.

import (
	"bufio"
	"drift/pkg/core"
	"fmt"
	"os"
	"sort"
	"strconv"
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
	NumSites   int // records written (founder + de-novo)
	Path       string

	// Merged-export breakdown (TMR4A.md W6); all zero unless
	// VCFOptions.IncludeMutations is set.
	NumFounderSites  int // records from the founder bitfield
	NumMutationSites int // records from the de-novo mutation pool
	MutationsIndexed int // distinct pool mutations carried by the sample
	// Carried mutation ids with no MutationPool entry, so no position: dropped, not
	// placed at position 0. Non-zero means mutation_count_model="legacy" GC'd a
	// still-carried lineage (utils/mutation.go).
	MutationsUnpooled int
	// Indexed mutations whose Position fell outside every chromosome arm, so no VCF
	// coordinate exists for them. Dropped and counted rather than silently absent.
	MutationsOutOfRange int
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

// hasPoolMutations reports whether an individual carries any de-novo mutation that
// can actually be placed — i.e. one whose id resolves to a Mutation in MutationPool.
// An id carried but unpooled has no Position (see buildMutationIndex), so it cannot
// contribute a site and must not make an otherwise-empty individual look exportable.
func hasPoolMutations(pop *core.Pop, id int) bool {
	strands := pop.IndMutations[id]
	for s := 0; s < 2; s++ {
		for _, mid := range strands[s] {
			if _, ok := pop.MutationPool[mid]; ok {
				return true
			}
		}
	}
	return false
}

// hasBitfield reports whether an individual has two allocated strand bitfields.
// Absence is not an error: pop.Chromosomes[child] is only allocated when a parent
// had set bits or presented a created germ cell (pkg/simulation/birth.go:123), so on
// a POOL-ONLY run (init_heterozygosity = 0, no created alleles) NOBODY has one while
// the mutation pool is full. Such an individual reads as ancestral at every founder
// site, which is exactly right — it inherited no founder allele because there were
// none to inherit.
func hasBitfield(pop *core.Pop, id int) bool {
	c, ok := pop.Chromosomes[id]
	return ok && len(c) >= 2 && c[0] != nil && c[1] != nil
}

// exportSamples picks the individuals a genotype export can describe, preserving the
// given (already sorted) id order so sample columns are stable. ONE definition shared
// by ExportVCF and ExportARGweaver, because drift.go draws a single sample for both
// and the two files are only joinable if they agree on who is in it.
//
// An individual qualifies on EITHER substrate: an allocated bitfield, or at least one
// placeable pool mutation. Requiring the bitfield — as this did before — made a
// pool-only run unexportable (ARGweaver errored with "no sampled individual has
// tracked chromosomes"; the VCF wrote a sample-less header) even though every site it
// needed was sitting in the pool. Worse, in the patchy case it silently DROPPED the
// individuals with no founder alleles, biasing the sample while still producing a
// well-formed file.
//
// withBitfield is how many of the returned samples have one, so a caller can skip the
// founder-site walk entirely when the answer is none.
func exportSamples(pop *core.Pop, ids []int) (samples []int, withBitfield int) {
	samples = make([]int, 0, len(ids))
	for _, id := range ids {
		bf := hasBitfield(pop, id)
		if !bf && !hasPoolMutations(pop, id) {
			continue
		}
		if bf {
			withBitfield++
		}
		samples = append(samples, id)
	}
	return samples, withBitfield
}

// ExportVCF writes the sampled individuals' genomes to a VCF v4.2 file in the
// model's results directory and returns a summary. Pass a SAMPLE of ids (see
// SampleIDs); individuals without tracked chromosomes are skipped. When
// opts.IncludeFixed is false only sites polymorphic within the sample are written;
// see VCFOptions (vcf_merged.go) for the W6 merged-export and truth-annotation
// switches. The zero VCFOptions reproduces the pre-W6 output byte for byte.
func ExportVCF(model *core.Model, pop *core.Pop, ids []int, opts VCFOptions) (*VCFStats, error) {
	// Individuals this export can describe, on either substrate (see exportSamples).
	samples, withBitfield := exportSamples(pop, ids)

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

	// De-novo pool overlay (TMR4A.md W6) and the W3 truth annotation, both built once
	// over the sample before any record is written. Nil when not requested, which is
	// what keeps the record loop below on its pre-W6 path.
	var mutIdx *mutationIndex
	if opts.IncludeMutations {
		mutIdx = buildMutationIndex(pop, samples)
		warnMergedLinkage(model, mutIdx)
	}
	var labels map[int][]string
	if opts.AlleleLabels {
		labels = buildLabelFields(model, pop, samples)
	}

	// --- Header ---
	fmt.Fprintln(w, "##fileformat=VCFv4.2")
	fmt.Fprintf(w, "##source=DRIFT model=%s run=%d year=%d\n", model.ModelName, run, year)
	fmt.Fprintln(w, "##FILTER=<ID=PASS,Description=\"All filters passed\">")
	fmt.Fprintln(w, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for the ALT allele\">")
	fmt.Fprintln(w, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">")
	if mutIdx != nil {
		fmt.Fprintln(w, "##INFO=<ID=SRC,Number=1,Type=String,Description=\"Variant substrate: founder (seeded bitfield) or denovo (mutation pool)\">")
	}
	fmt.Fprintln(w, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	if labels != nil {
		fmt.Fprintln(w, "##FORMAT=<ID=FL,Number=2,Type=Integer,Description=\"Ground-truth created-allele label (TMR4A W3) of each phased strand at this focal locus; . = not recoverable\">")
	}
	fmt.Fprintln(w, "##ALT allele encodes the founder-derived (bit=1) state; REF the ancestral (bit=0) state; letters are placeholders.")
	if mutIdx != nil {
		fmt.Fprintln(w, "##Records are keyed by IDENTITY, not coordinate: F<bit> is a founder bitfield site, M<id> a de-novo mutation. Positions recur, so several records may share one POS — that is the infinite-sites model, not a duplicate.")
	}
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
	// INFO suffix and record ids only exist in merged mode; hoisted so the default
	// path formats exactly the string it always did.
	srcFounder, srcDenovo := "", ""
	if mutIdx != nil {
		srcFounder, srcDenovo = ";SRC=founder", ";SRC=denovo"
	}

	// --- Records ---
	// Reused genotype-string buffer to avoid per-site allocation.
	gt := make([]byte, 0, 8*len(samples))
	stats := &VCFStats{NumSamples: len(samples), Path: path}
	mutSeen := 0
	for _, ct := range contigs {
		for _, arm := range ct.arms {
			start, length := arm[0], arm[1]
			for bit := start; bit < start+length; bit++ {
				pos := bit - ct.base + 1
				// The FL annotation is a property of the LOCUS, so both the founder
				// record and any de-novo record at this position carry it.
				fl := labels[bit]
				format := "GT"
				if fl != nil {
					format = "GT:FL"
				}

				// --- founder bitfield site ---
				// Skipped wholesale on a pool-only run: with no bitfield anywhere there
				// are no founder sites, and emitting one all-reference record per genome
				// bit under IncludeFixed would bury the pool records under millions of
				// monomorphic rows.
				ac := 0
				gt = gt[:0]
				if withBitfield > 0 {
					for si, id := range samples {
						// A sample with no allocated bitfield reads as ancestral here; it
						// inherited no founder allele. bitAt already tolerates a short
						// slice, so only the c[0]/c[1] indexing needs the length guard.
						a0, a1 := false, false
						if c := pop.Chromosomes[id]; len(c) >= 2 {
							a0, a1 = bitAt(c[0], bit), bitAt(c[1], bit)
						}
						gt = appendSample(gt, a0, a1, fl, si)
						if a0 {
							ac++
						}
						if a1 {
							ac++
						}
					}
				}
				if withBitfield > 0 && (opts.IncludeFixed || (ac != 0 && ac != an)) {
					recID := "."
					if mutIdx != nil {
						recID = "F" + strconv.Itoa(bit)
					}
					writeVCFRecord(w, ct.chrom, pos, recID, srcFounder, ac, an, format, gt)
					stats.NumSites++
					stats.NumFounderSites++
				}

				// --- de-novo pool sites at this position (W6) ---
				if mutIdx == nil {
					continue
				}
				for _, mid := range mutIdx.byPos[bit] {
					mutSeen++
					gt, ac = mutIdx.appendGenotypes(gt[:0], mid, len(samples), fl)
					if !opts.IncludeFixed && (ac == 0 || ac == an) {
						continue
					}
					writeVCFRecord(w, ct.chrom, pos, "M"+strconv.Itoa(mid), srcDenovo, ac, an, format, gt)
					stats.NumSites++
					stats.NumMutationSites++
				}
			}
		}
	}

	if mutIdx != nil {
		stats.MutationsIndexed = mutIdx.total
		stats.MutationsUnpooled = mutIdx.unpooled
		stats.MutationsOutOfRange = mutIdx.total - mutSeen
		fmt.Printf("VCF: %d sites × %d samples (%d founder + %d de-novo of %d carried) → %s\n",
			stats.NumSites, len(samples), stats.NumFounderSites, stats.NumMutationSites,
			mutIdx.total, path)
		if stats.MutationsUnpooled > 0 || stats.MutationsOutOfRange > 0 {
			fmt.Printf("VCF: dropped %d unpooled and %d out-of-range mutations (no position to write them at)\n",
				stats.MutationsUnpooled, stats.MutationsOutOfRange)
		}
	} else {
		fmt.Printf("VCF: %d sites × %d samples → %s\n", stats.NumSites, len(samples), path)
	}
	return stats, nil
}

// writeVCFRecord emits one record line. With recID "." and src "" this is exactly the
// line the pre-W6 exporter wrote.
func writeVCFRecord(w *bufio.Writer, chrom, pos int, recID, src string, ac, an int, format string, gt []byte) {
	fmt.Fprintf(w, "%d\t%d\t%s\tA\tT\t.\tPASS\tAC=%d;AN=%d%s\t%s", chrom, pos, recID, ac, an, src, format)
	w.Write(gt)
	fmt.Fprintln(w)
}

// warnMergedLinkage says loudly when a merged export is about to describe haplotypes
// no meiosis produced. Under linkage_model="legacy" (the DEFAULT) the de-novo pool
// uses the opposite meiosis-mask polarity from the founder bitfield, so the two
// substrates ANTI-segregate and their combination is an artefact — see vcf_merged.go.
// A warning rather than a refusal: the merge is still the right thing to do on an
// arm-linkage run, and a legacy run's file is diagnostically useful as long as nobody
// mistakes it for a genealogy.
func warnMergedLinkage(model *core.Model, idx *mutationIndex) {
	if idx == nil || idx.total == 0 {
		fmt.Println("WARNING: vcf_include_mutations is set but the sample carries no pool mutations; " +
			"is track_mutations on?")
		return
	}
	if model.StringParam("linkage_model", "legacy") != "arm" {
		fmt.Println("WARNING: merged VCF export under linkage_model=\"legacy\": de-novo mutations " +
			"anti-segregate with the founder bitfield, so the exported haplotypes are an artefact " +
			"of the mixed polarity. Set linkage_model=\"arm\" for anything fed to an ARG inference tool.")
	}
}

func appendAllele(b []byte, set bool) []byte {
	if set {
		return append(b, '1')
	}
	return append(b, '0')
}
