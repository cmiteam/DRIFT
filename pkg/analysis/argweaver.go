package analysis

// ARGweaver export — TMR4A.md work item W7, the half that does not need ARGweaver
// installed: the converter from DRIFT's sampled genomes to ARGweaver's native `.sites`
// format, plus the region BED, the scaling rationale (argweaver_scaling.go) and a ready
// arg-sample command line carrying the derived rates.
//
// NOT REIMPLEMENTING THE MCMC. README 5b and W7 are explicit: run the real tool. This
// module's entire job is to hand it correctly-scaled input and to record, in the output
// directory, exactly what scaling was claimed. The parser for ARGweaver's output — and
// the join against the W5 ground truth — is the other half of W7 and is not here.
//
// THE FORMAT (verified against the ARGweaver manual, not assumed):
//
//	NAMES   <hap1>  <hap2>  ...              one name per HAPLOID genome, two per individual
//	REGION  <chrom> <start> <end>            1-based, END-INCLUSIVE
//	<pos>   <alleles>                        1-based absolute; one base per haplotype
//
// Tab-delimited; sites must be sorted by coordinate; any position not listed is taken to
// be invariant across the samples. So this writes SEGREGATING SITES ONLY regardless of
// vcf_include_fixed — a site fixed within the sample is invariant, which the REGION line
// already accounts for, and listing it would only inflate the file.
//
// ONE FILE PER CHROMOSOME, because arg-sample runs one region at a time and because the
// realized recombination rate is chromosome-specific (see argweaver_scaling.go: the
// obligate crossover makes a short chromosome recombine far more per bp than a long one).
// A genome-wide average rho would be wrong for every region it was applied to.
//
// COLLISIONS ARE THE REASON bp-PER-BIT MUST EXCEED 1. DRIFT positions recur — the founder
// bit at position P and any number of de-novo mutations can all sit at P — and §8c
// measured 213 of 223 de-novo records sharing a coordinate with a founder record in a
// routine run. The `.sites` format keys on position and requires sorted coordinates, so
// those sites need DISTINCT bp. Each genome bit owns a bp window of width bp_per_bit and
// colliding sites are laid out consecutively inside it. Over ≤100 bp the recombination
// this implies between them is ~1e-6 per generation, i.e. nothing; the alternative —
// dropping all but one site per bit — would discard most of the mutational history the
// export exists to carry. A bit whose window overflows drops the excess and counts it.
//
// THE SAME SITES AS THE W6 VCF, BY CONSTRUCTION: same sample, same contig walk, same
// segregating filter, same buildMutationIndex. TestARGweaverSitesMatchVCFRecords asserts
// the two files list the same sites in the same order, because the whole study is a join
// between what ARGweaver infers here and what the W5 walker knows, and a silent
// divergence between the two renderings would corrupt it invisibly.

import (
	"bufio"
	"drift/pkg/core"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// Allele characters. DRIFT bit 0 is ancestral and bit 1 derived; the letters are
// placeholders exactly as in the VCF's REF=A/ALT=T convention, so the two exports agree
// base for base and a truth join can be made on either.
const (
	argAncestral = 'A'
	argDerived   = 'T'
)

// ARGweaverOptions configures an export pass.
type ARGweaverOptions struct {
	// BpPerBit is the declared bp-per-genome-bit convention (see argweaver_scaling.go).
	// It sets resolution and the absolute per-bp rates, and cancels out of theta/rho.
	BpPerBit int
	// Chroms restricts the export to these chromosomes; empty exports all. A real study
	// scopes to one arm (§1A: one arm is larger than ARGweaver's native ~2 Mb block).
	Chroms []int
}

// ARGweaverOptionsFromModel reads the export settings off the model parameters.
func ARGweaverOptionsFromModel(model *core.Model) ARGweaverOptions {
	opts := ARGweaverOptions{BpPerBit: 100}
	if model == nil {
		return opts
	}
	if v := int(model.Parameters["argweaver_bp_per_bit"]); v > 0 {
		opts.BpPerBit = v
	}
	for _, f := range strings.FieldsFunc(model.StringParam("argweaver_chroms", ""), func(r rune) bool {
		return r == ',' || r == ' ' || r == '\t' || r == ';'
	}) {
		if c, err := strconv.Atoi(strings.TrimSpace(f)); err == nil {
			opts.Chroms = append(opts.Chroms, c)
		} else {
			fmt.Printf("WARNING: argweaver_chroms: bad chromosome %q, skipped\n", f)
		}
	}
	return opts
}

// ARGweaverRegion is one exported chromosome.
type ARGweaverRegion struct {
	Chrom   int
	Path    string
	StartBp int // 1-based inclusive, as written on the REGION line
	EndBp   int
	Sites   int
	Dropped int // sites whose bit's bp window overflowed
}

// ARGweaverStats summarizes an export pass.
type ARGweaverStats struct {
	Regions       []ARGweaverRegion
	NumHaplotypes int
	TotalSites    int
	TotalDropped  int
	BedPath       string
	RationalePath string
	CommandPath   string
	Anchors       *RateAnchors
}

// argSite is one biallelic site awaiting a bp coordinate: either a founder bitfield bit
// or a pool mutation, distinguished exactly as in the W6 VCF.
type argSite struct {
	bit   int
	mutID int // -1 for a founder bitfield site
}

// ExportARGweaver writes one `.sites` file per exported chromosome plus the region BED,
// the scaling rationale, and a runnable arg-sample command line. Pass the SAME sample the
// W6 VCF and the W5 walker were given.
func ExportARGweaver(model *core.Model, pop *core.Pop, ids []int, opts ARGweaverOptions) (*ARGweaverStats, error) {
	if opts.BpPerBit < 1 {
		opts.BpPerBit = 1
	}
	samples := make([]int, 0, len(ids))
	for _, id := range ids {
		if c, ok := pop.Chromosomes[id]; ok && len(c) >= 2 && c[0] != nil && c[1] != nil {
			samples = append(samples, id)
		}
	}
	if len(samples) == 0 {
		return nil, fmt.Errorf("ARGweaver export: no sampled individual has tracked chromosomes")
	}

	contigs := buildContigs(model)
	if len(opts.Chroms) > 0 {
		want := make(map[int]bool, len(opts.Chroms))
		for _, c := range opts.Chroms {
			want[c] = true
		}
		kept := contigs[:0:0]
		for _, ct := range contigs {
			if want[ct.chrom] {
				kept = append(kept, ct)
			}
		}
		contigs = kept
	}
	if len(contigs) == 0 {
		return nil, fmt.Errorf("ARGweaver export: argweaver_chroms selected no chromosome that exists")
	}

	// The pool is always merged: ARGweaver needs the segregating sites INCLUDING the
	// de-novo mutations, which is the whole reason W6 exists (TMR4A.md W7).
	mutIdx := buildMutationIndex(pop, samples)
	if mutIdx.total == 0 {
		fmt.Println("WARNING: ARGweaver export: the sample carries no de-novo mutations. Only founder " +
			"bitfield sites will be written, so the inferred genealogy will rest on seeded variation " +
			"alone — is track_mutations on?")
	}

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	prefix := fmt.Sprintf("%s/%s_run%d_year%d", model.ResultsDir, model.ModelName, run, year)

	stats := &ARGweaverStats{
		NumHaplotypes: 2 * len(samples),
		Anchors:       ComputeRateAnchors(model, opts.BpPerBit, contigs),
	}

	for _, ct := range contigs {
		region, err := writeSitesFile(prefix, ct, model, pop, samples, mutIdx, opts.BpPerBit)
		if err != nil {
			return nil, err
		}
		stats.Regions = append(stats.Regions, *region)
		stats.TotalSites += region.Sites
		stats.TotalDropped += region.Dropped
	}

	var err error
	if stats.BedPath, err = writeRegionBed(prefix, stats.Regions); err != nil {
		return nil, err
	}
	if stats.RationalePath, err = writeRationaleFile(prefix, model, stats.Anchors); err != nil {
		return nil, err
	}
	if stats.CommandPath, err = writeCommandFile(prefix, stats, len(pop.IndData)); err != nil {
		return nil, err
	}

	fmt.Printf("ARGweaver: %d regions, %d sites x %d haplotypes at %d bp/bit -> %s_chr*.sites\n",
		len(stats.Regions), stats.TotalSites, stats.NumHaplotypes, opts.BpPerBit, prefix)
	if stats.TotalDropped > 0 {
		fmt.Printf("ARGweaver: dropped %d sites whose genome bit had more colliding sites than its "+
			"%d-bp window could hold; raise argweaver_bp_per_bit\n", stats.TotalDropped, opts.BpPerBit)
	}
	for _, msg := range stats.Anchors.Warnings {
		fmt.Printf("WARNING: ARGweaver scaling: %s\n", msg)
	}
	return stats, nil
}

// writeSitesFile emits one chromosome's `.sites`. The site walk mirrors ExportVCF's
// exactly — contig arms in order, founder bit then its de-novo mutations by ascending id,
// segregating within the sample only — so the two files describe the same site set.
func writeSitesFile(prefix string, ct vcfContig, model *core.Model, pop *core.Pop,
	samples []int, mutIdx *mutationIndex, bpPerBit int) (*ARGweaverRegion, error) {

	path := fmt.Sprintf("%s_chr%d.sites", prefix, ct.chrom)
	file, err := os.Create(path)
	if err != nil {
		return nil, fmt.Errorf("create sites file: %w", err)
	}
	defer file.Close()
	w := bufio.NewWriter(file)
	defer w.Flush()

	region := &ARGweaverRegion{
		Chrom:   ct.chrom,
		Path:    path,
		StartBp: 1,
		EndBp:   ct.span * bpPerBit, // 1-based end-INCLUSIVE, per the manual
	}

	// NAMES: one per haploid genome. IND<id>_1 / _2 line up with the VCF's IND<id>
	// column and its two phased strands, which is what makes the truth join mechanical.
	fmt.Fprint(w, "NAMES")
	for _, id := range samples {
		fmt.Fprintf(w, "\tIND%d_1\tIND%d_2", id, id)
	}
	fmt.Fprintln(w)
	fmt.Fprintf(w, "REGION\t%d\t%d\t%d\n", ct.chrom, region.StartBp, region.EndBp)

	an := 2 * len(samples)
	alleles := make([]byte, an)
	for _, arm := range ct.arms {
		start, length := arm[0], arm[1]
		for bit := start; bit < start+length; bit++ {
			// Every site at this bit, in the VCF's order, then laid out consecutively
			// inside the bit's bp window.
			sites := collectSitesAtBit(bit, mutIdx)
			if len(sites) == 0 {
				continue
			}
			bitBp := (bit-ct.base)*bpPerBit + 1
			offset := 0
			for _, s := range sites {
				ac := fillAlleles(alleles, s, pop, samples, mutIdx)
				if ac == 0 || ac == an {
					continue // invariant in the sample; the REGION line covers it
				}
				if offset >= bpPerBit {
					region.Dropped++
					continue
				}
				fmt.Fprintf(w, "%d\t%s\n", bitBp+offset, alleles)
				offset++
				region.Sites++
			}
		}
	}
	return region, nil
}

// collectSitesAtBit lists the sites sharing one genome bit: the founder bitfield site
// first, then each de-novo mutation by ascending id. Same order as the VCF's records, and
// deterministic because buildMutationIndex sorted the ids.
func collectSitesAtBit(bit int, mutIdx *mutationIndex) []argSite {
	muts := mutIdx.byPos[bit]
	sites := make([]argSite, 0, len(muts)+1)
	sites = append(sites, argSite{bit: bit, mutID: -1})
	for _, mid := range muts {
		sites = append(sites, argSite{bit: bit, mutID: mid})
	}
	return sites
}

// fillAlleles writes one site's per-haplotype bases into buf (haplotype 2i+s is strand s
// of sample i, matching the NAMES order) and returns the derived-allele count.
func fillAlleles(buf []byte, s argSite, pop *core.Pop, samples []int, mutIdx *mutationIndex) int {
	ac := 0
	if s.mutID < 0 {
		for i, id := range samples {
			c := pop.Chromosomes[id]
			for strand := 0; strand < 2; strand++ {
				if bitAt(c[strand], s.bit) {
					buf[2*i+strand] = argDerived
					ac++
				} else {
					buf[2*i+strand] = argAncestral
				}
			}
		}
		return ac
	}
	for i := range buf {
		buf[i] = argAncestral
	}
	for _, carrier := range mutIdx.carriers[s.mutID] {
		for strand := 0; strand < 2; strand++ {
			if carrier.strands&(1<<uint(strand)) != 0 {
				buf[2*carrier.sample+strand] = argDerived
				ac++
			}
		}
	}
	return ac
}

// writeRegionBed emits the exported spans as a 3-column BED, the form ARGweaver's
// --maskmap and the usual downstream tooling read. BED is 0-based half-open, so it says
// [0, len) where the REGION line says [1, len].
func writeRegionBed(prefix string, regions []ARGweaverRegion) (string, error) {
	path := prefix + "_regions.bed"
	file, err := os.Create(path)
	if err != nil {
		return "", fmt.Errorf("create region BED: %w", err)
	}
	defer file.Close()
	w := bufio.NewWriter(file)
	defer w.Flush()
	for _, r := range regions {
		fmt.Fprintf(w, "%d\t%d\t%d\n", r.Chrom, r.StartBp-1, r.EndBp)
	}
	return path, nil
}

func writeRationaleFile(prefix string, model *core.Model, a *RateAnchors) (string, error) {
	path := prefix + "_argweaver_scaling.txt"
	file, err := os.Create(path)
	if err != nil {
		return "", fmt.Errorf("create scaling rationale: %w", err)
	}
	defer file.Close()
	w := bufio.NewWriter(file)
	defer w.Flush()
	a.WriteRationale(w, model)
	return path, nil
}

// writeCommandFile emits a runnable arg-sample invocation per region with the derived
// rates already filled in. The rates are the single easiest thing to get wrong at the
// command line — they are two per-generation counts divided by a declared convention, and
// nothing about the .sites file itself records them — so they ship next to the data.
func writeCommandFile(prefix string, stats *ARGweaverStats, censusN int) (string, error) {
	path := prefix + "_argweaver_cmd.sh"
	file, err := os.Create(path)
	if err != nil {
		return "", fmt.Errorf("create command file: %w", err)
	}
	defer file.Close()
	w := bufio.NewWriter(file)
	defer w.Flush()

	a := stats.Anchors
	fmt.Fprintln(w, "#!/bin/sh")
	fmt.Fprintln(w, "# Generated by DRIFT (TMR4A.md W7). Rates are derived in the companion")
	fmt.Fprintf(w, "# %s_argweaver_scaling.txt -- read it before trusting any output.\n", basename(prefix))
	fmt.Fprintln(w, "#")
	fmt.Fprintf(w, "# --mutrate is genome-wide (%.6g per site per generation); --recombrate is\n", a.MuPerBp)
	fmt.Fprintln(w, "# PER REGION, because DRIFT's obligate crossover makes it chromosome-specific.")
	fmt.Fprintln(w, "#")
	fmt.Fprintln(w, "# --popsize is the ASSUMED effective size and is the primary object of study")
	fmt.Fprintln(w, "# (TMR4A.md §1B / W9): the constant-N coalescent predicts TMR4A ~ N generations")
	fmt.Fprintln(w, "# before any data is consulted, so sweep it rather than trusting one value.")
	fmt.Fprintf(w, "# DRIFT's own CENSUS size at export was %d diploids -- a starting point, not a\n", censusN)
	fmt.Fprintln(w, "# recommendation, and NOT the same quantity as a coalescent Ne.")
	fmt.Fprintln(w)
	fmt.Fprintf(w, "POPSIZE=${POPSIZE:-%d}\n", censusN)
	fmt.Fprintln(w)
	for _, r := range stats.Regions {
		rho := 0.0
		for _, rr := range a.Regions {
			if rr.Chrom == r.Chrom {
				rho = rr.RhoPerBp
			}
		}
		fmt.Fprintf(w, "arg-sample --sites %s \\\n", basename(r.Path))
		fmt.Fprintf(w, "  --popsize \"$POPSIZE\" --mutrate %.6g --recombrate %.6g \\\n", a.MuPerBp, rho)
		fmt.Fprintf(w, "  --output %s_chr%d\n\n", basename(prefix), r.Chrom)
	}
	return path, nil
}

func basename(p string) string {
	if i := strings.LastIndexAny(p, "/\\"); i >= 0 {
		return p[i+1:]
	}
	return p
}
