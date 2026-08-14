package analysis

// Merged VCF export — TMR4A.md work item W6, the prerequisite for the ARGweaver
// bridge (W7).
//
// WHAT WAS MISSING. ExportVCF (vcf.go) writes the founder BITFIELD only. DRIFT
// carries variation in two substrates — the per-individual bitfield and the
// de-novo mutation POOL (pop.IndMutations + pop.MutationPool) — and an external
// ARG inference tool needs the segregating sites from BOTH. A VCF holding only
// the founder bits describes a sample whose entire mutational history since
// founding is invisible, which is precisely the history TMR4A is claimed to date.
//
// WHY THE POOL IS NOT "OVERLAID" ONTO THE BITFIELD. TMR4A.md W6 says "overlay the
// pool onto the bitfield per sampled strand", and read as a literal OR that is
// wrong three separate ways, each of which destroys identity — §6 risk 4, "identity
// must key on mutation ID and founder label throughout, never on coordinate":
//
//  1. GenerateNewMutations draws RandIntn(genome_bits), so positions RECUR. Two
//     distinct mutations at the same position would merge into one site and their
//     two independent genealogies would be reported as one.
//  2. A mutation landing where the founder bit is already 1 would vanish entirely:
//     1 OR 1 = 1, no new site, one lost segregating site.
//  3. A founder-derived allele and a de-novo allele at the same coordinate are
//     different alleles with different histories. Collapsing them asserts a
//     coalescence that never happened.
//
// So each source site gets its OWN record, keyed by identity: the founder bit at
// position P is record F<P>, and de-novo mutation M is record M<id>, both mapped to
// the same POS. Multiple records at one POS are legal VCF and are exactly what an
// infinite-sites model means. Downstream joins key on the ID column, never on POS.
//
// LINKAGE — THE ONE THING THAT MUST BE SET RIGHT. Under the DEFAULT
// linkage_model="legacy" the pool uses the OPPOSITE meiosis-mask polarity from the
// bitfield (see utils/mutation.go), so a de-novo mutation at position P follows the
// opposite parental strand from the founder allele at P: the two substrates
// ANTI-segregate. A merged export from a legacy run therefore describes haplotypes
// that no meiosis ever produced, and handing that to an ARG inference tool would
// measure the artefact rather than the model. linkage_model="arm" aligns the
// polarities and is REQUIRED for a coherent merged export; ExportVCF warns when it
// is not set. This is the W6 counterpart of §8a's finding that anything fed to
// ARGweaver must use recombination_model="map".
//
// IMPORTED PANELS ARE ALREADY MERGED. ImportVCF deliberately mirrors every site
// into BOTH substrates (one genome bit AND one pool Mutation at the same Position)
// so that bitfield-reading and pool-reading analyses both see the panel. Merging
// again on the way back out would emit every site twice, so drift.go's
// re-export-after-import path forces IncludeMutations off and says so.
//
// TRUTH ANNOTATION. With vcf_allele_labels on, records at a focal locus carry a
// FORMAT FL field: the W3 created-allele label of each phased strand, resolved by
// walking the recorded genealogy (W1+W2+W4) back to its terminus. That is the join
// key for the truth-vs-inference plot — it puts DRIFT's ground-truth allele identity
// in the same file as the sequence an inference tool will read, per sample, per
// strand, so no separate id mapping has to be maintained. An unresolvable strand
// emits '.', never a shared default: assuming shared identity in the absence of
// evidence is the exact error the labels exist to prevent.

import (
	"drift/pkg/core"
	"fmt"
	"sort"
	"strconv"
)

// VCFOptions selects what an export pass writes. Zero value = the pre-W6 behaviour
// (founder bitfield only, segregating sites only, no annotation), so a default run
// is byte-identical.
type VCFOptions struct {
	// IncludeFixed emits every position, not only those segregating in the sample.
	IncludeFixed bool
	// IncludeMutations adds one record per de-novo pool mutation carried by the
	// sample (W6). Requires track_mutations to have run, and linkage_model="arm"
	// for the result to be internally coherent.
	IncludeMutations bool
	// AlleleLabels adds the FORMAT FL created-allele truth annotation at focal loci
	// (W3 labels resolved through the W2 provenance records).
	AlleleLabels bool
}

// VCFOptionsFromModel reads the export flags off the model parameters.
func VCFOptionsFromModel(model *core.Model) VCFOptions {
	if model == nil {
		return VCFOptions{}
	}
	return VCFOptions{
		IncludeFixed:     model.Parameters["vcf_include_fixed"] == 1,
		IncludeMutations: model.Parameters["vcf_include_mutations"] == 1,
		AlleleLabels:     model.Parameters["vcf_allele_labels"] == 1,
	}
}

// mutCarrier is one sampled individual carrying a mutation, and on which strands.
// `sample` is an index into the export's SAMPLE COLUMN ORDER, not an individual id,
// so emission is a single forward merge against the column loop.
type mutCarrier struct {
	sample  int
	strands uint8 // bit s set = carried on strand s
}

// mutationIndex is the de-novo pool restricted to the exported sample, keyed for
// position-ordered emission.
type mutationIndex struct {
	byPos    map[int][]int        // genome-bit position -> mutation ids, ascending
	carriers map[int][]mutCarrier // mutation id -> carriers, in column order
	total    int                  // distinct mutations indexed
	unpooled int                  // ids carried by the sample but absent from MutationPool
}

// buildMutationIndex collects the de-novo mutations carried by the sample.
//
// UNPOOLED IDS ARE COUNTED, NOT GUESSED. An id present in IndMutations but missing
// from MutationPool has no Position, so it cannot be placed. Under the DEFAULT
// mutation_count_model="refcount" this is zero; under the opt-in "legacy"
// asymmetric refcount a still-carried lineage can be garbage-collected while alive
// (utils/mutation.go), and the old code path silently resolved it to the zero
// Mutation at position 0. Dropping it and REPORTING the count is the honest
// alternative — a merged export that quietly parked mutations at position 0 would
// manufacture linkage exactly where §6f found it before.
func buildMutationIndex(pop *core.Pop, samples []int) *mutationIndex {
	idx := &mutationIndex{byPos: map[int][]int{}, carriers: map[int][]mutCarrier{}}
	if pop == nil || pop.IndMutations == nil {
		return idx
	}
	unpooled := map[int]bool{}
	for si, id := range samples {
		strands := pop.IndMutations[id]
		if len(strands) == 0 {
			continue
		}
		for s := 0; s < 2; s++ {
			for _, mid := range strands[s] {
				mut, ok := pop.MutationPool[mid]
				if !ok {
					unpooled[mid] = true
					continue
				}
				list, seen := idx.carriers[mid]
				if !seen {
					idx.byPos[mut.Position] = append(idx.byPos[mut.Position], mid)
					idx.total++
				}
				// Both strands of one individual are visited consecutively, so a
				// homozygous carrier extends the entry just appended rather than
				// adding a second one for the same column.
				if n := len(list); n > 0 && list[n-1].sample == si {
					list[n-1].strands |= 1 << uint(s)
				} else {
					list = append(list, mutCarrier{sample: si, strands: 1 << uint(s)})
				}
				idx.carriers[mid] = list
			}
		}
	}
	// Ascending mutation id within a position: map order is randomized and the
	// record order must be a property of the data, not of the run.
	for pos := range idx.byPos {
		sort.Ints(idx.byPos[pos])
	}
	idx.unpooled = len(unpooled)
	return idx
}

// appendGenotypes writes one mutation's phased genotype column for every sample and
// returns the buffer and the ALT allele count. The carrier list is in column order,
// so this is a linear merge rather than a per-sample lookup.
func (idx *mutationIndex) appendGenotypes(gt []byte, mid int, nSamples int, fl []string) ([]byte, int) {
	list := idx.carriers[mid]
	ac, c := 0, 0
	for si := 0; si < nSamples; si++ {
		var m uint8
		if c < len(list) && list[c].sample == si {
			m = list[c].strands
			c++
		}
		a0, a1 := m&1 != 0, m&2 != 0
		gt = appendSample(gt, a0, a1, fl, si)
		if a0 {
			ac++
		}
		if a1 {
			ac++
		}
	}
	return gt, ac
}

// appendSample writes one sample's phased genotype, plus the FL truth annotation
// when this position is a focal locus.
func appendSample(gt []byte, a0, a1 bool, fl []string, si int) []byte {
	gt = append(gt, '\t')
	gt = appendAllele(gt, a0)
	gt = append(gt, '|')
	gt = appendAllele(gt, a1)
	if fl != nil {
		gt = append(gt, ':')
		gt = append(gt, fl[si]...)
	}
	return gt
}

// buildLabelFields precomputes the FORMAT FL value for every sample at every focal
// locus: the created-allele label each of the sample's two strands carries there,
// as a Number=2 Integer field ("0,3"), with '.' where no label is recoverable.
//
// Keyed by GENOME-BIT position, because that is what the record loop iterates.
// Returns nil (no annotation, and a warning) when the substrate is absent — a file
// full of ".,." would look like a measurement of nothing rather than like a missing
// prerequisite.
func buildLabelFields(model *core.Model, pop *core.Pop, samples []int) map[int][]string {
	loci := core.EnsureARGLoci(model)
	if len(loci) == 0 || pop == nil || pop.Pedigree == nil || !pop.FounderLabelsOn() {
		fmt.Println("WARNING: vcf_allele_labels needs track_pedigree + track_arg/arg_loci + " +
			"track_allele_labels (or the created founder model); no FL annotation written")
		return nil
	}
	out := make(map[int][]string, len(loci))
	for i, pos := range loci {
		byInd := FounderAlleleLabelsAt(pop, samples, i, len(loci))
		fields := make([]string, len(samples))
		for si, id := range samples {
			pair, ok := byInd[id]
			if !ok {
				fields[si] = ".,."
				continue
			}
			fields[si] = labelField(pair[0]) + "," + labelField(pair[1])
		}
		out[pos] = fields
	}
	return out
}

func labelField(label int) string {
	if label < 0 {
		return "."
	}
	return strconv.Itoa(label)
}
