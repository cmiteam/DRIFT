package analysis

import (
	"drift/pkg/core"
	"os"
	"strings"
	"testing"
)

// Merged VCF export (TMR4A.md W6). The tests that matter here are the three ways a
// literal "OR the pool onto the bitfield" overlay loses a site — recurrent positions,
// a mutation under an already-set founder bit, and the founder/de-novo distinction —
// because those are silent losses: the file still parses, still looks complete, and
// simply describes fewer genealogies than the run contained.

// mergedFixture extends the shared IBD fixture (one 128-bit chromosome, p arm [0,64),
// q arm [64,128); inds 1 & 2 carry the whole p arm on copy 0, ind 3 its low 32 bits,
// ind 4 nothing) with a de-novo pool chosen to exercise every drop path at once.
func mergedFixture() (*core.Model, *core.Pop) {
	model := makeIBDTestModel()
	model.StringParams = map[string]string{"linkage_model": "arm"}
	pop := makeIBDTestPop()

	pop.MutationPool = map[int]core.Mutation{
		// 100 and 101 share position 70. GenerateNewMutations draws RandIntn(genome_bits)
		// so positions RECUR; an overlay keyed on coordinate would merge these two
		// independent genealogies into one site.
		100: {Id: 100, Position: 70},
		101: {Id: 101, Position: 70},
		// 102 sits at position 5, where every sampled copy-0 founder bit is ALREADY 1.
		// 1 OR 1 = 1, so an overlay would erase this site entirely.
		102: {Id: 102, Position: 5},
		// 103 is outside every chromosome arm: no VCF coordinate exists for it.
		103: {Id: 103, Position: 200},
		// 105 is carried by every sampled strand — fixed in the sample, so it is dropped
		// by the same rule that drops a monomorphic founder bit.
		105: {Id: 105, Position: 90},
		// 104 is deliberately ABSENT from the pool while still being carried (see below).
	}
	pop.IndMutations = map[int]map[int][]int{
		1: {0: {100, 102, 103, 105}, 1: {105}},
		2: {0: {104, 105}, 1: {101, 105}},
		3: {0: {105}, 1: {100, 105}},
	}
	return model, pop
}

// The core W6 behaviour: pool mutations become their own records, keyed by mutation id,
// at the same POS as the founder bit they sit under.
func TestExportVCFMergedEmitsPoolSites(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportVCF(model, pop, []int{1, 2, 3, 4},
		VCFOptions{IncludeMutations: true})
	if err != nil {
		t.Fatalf("ExportVCF failed: %v", err)
	}

	// 64 segregating founder bits (as in the founder-only test) + 3 emitted mutations:
	// 100 (2 carriers), 101 (1), 102 (1). 105 is fixed in the sample, 103 is out of
	// range, 104 has no pool entry.
	if stats.NumFounderSites != 64 || stats.NumMutationSites != 3 || stats.NumSites != 67 {
		t.Errorf("sites = %d (%d founder + %d de-novo), want 67 (64 + 3)",
			stats.NumSites, stats.NumFounderSites, stats.NumMutationSites)
	}
	if stats.MutationsIndexed != 5 {
		t.Errorf("MutationsIndexed = %d, want 5 (100,101,102,103,105)", stats.MutationsIndexed)
	}
	if stats.MutationsUnpooled != 1 {
		t.Errorf("MutationsUnpooled = %d, want 1 (id 104 is carried but has no pool entry, "+
			"so it has no position and must be dropped, not written at position 0)",
			stats.MutationsUnpooled)
	}
	if stats.MutationsOutOfRange != 1 {
		t.Errorf("MutationsOutOfRange = %d, want 1 (id 103 at position 200)", stats.MutationsOutOfRange)
	}

	header, records := readVCF(t, stats.Path)
	if !containsLine(header, "##INFO=<ID=SRC,Number=1,Type=String,Description=\"Variant substrate: founder (seeded bitfield) or denovo (mutation pool)\">") {
		t.Errorf("merged export must declare the SRC INFO field:\n%s", strings.Join(header, "\n"))
	}
	byID := indexRecords(t, records)

	// The already-set-bit case. Founder bit 5 is 1 on copy 0 for all three samples
	// (AC=3); mutation 102 sits at the same position on ind 1's strand 0 only (AC=1).
	// BOTH records exist, at the same POS, distinguishable only by their ids.
	f5, ok := byID["F5"]
	if !ok {
		t.Fatalf("no founder record F5")
	}
	m102, ok := byID["M102"]
	if !ok {
		t.Fatalf("no de-novo record M102 — a mutation under an already-set founder bit was lost")
	}
	if f5[1] != "6" || m102[1] != "6" {
		t.Errorf("POS of F5/M102 = %s/%s, want 6/6 (same coordinate, different alleles)", f5[1], m102[1])
	}
	if f5[7] != "AC=3;AN=6;SRC=founder" {
		t.Errorf("F5 INFO = %s, want AC=3;AN=6;SRC=founder", f5[7])
	}
	if m102[7] != "AC=1;AN=6;SRC=denovo" {
		t.Errorf("M102 INFO = %s, want AC=1;AN=6;SRC=denovo", m102[7])
	}
	if m102[9] != "1|0" || m102[10] != "0|0" || m102[11] != "0|0" {
		t.Errorf("M102 genotypes = %v, want 1|0 0|0 0|0", m102[9:12])
	}

	// The recurrence case: two distinct mutations at position 70, emitted as two records
	// at POS 71 with different carriers. Merging them on coordinate would have reported
	// one site with AC=3 and one genealogy where there are two.
	m100, ok100 := byID["M100"]
	m101, ok101 := byID["M101"]
	if !ok100 || !ok101 {
		t.Fatalf("recurrent position 70 lost a mutation: M100 present=%v, M101 present=%v", ok100, ok101)
	}
	if m100[1] != "71" || m101[1] != "71" {
		t.Errorf("POS of M100/M101 = %s/%s, want 71/71", m100[1], m101[1])
	}
	// 100 is on ind 1 strand 0 and ind 3 strand 1; 101 on ind 2 strand 1 only.
	if m100[9] != "1|0" || m100[10] != "0|0" || m100[11] != "0|1" {
		t.Errorf("M100 genotypes = %v, want 1|0 0|0 0|1", m100[9:12])
	}
	if m101[9] != "0|0" || m101[10] != "0|1" || m101[11] != "0|0" {
		t.Errorf("M101 genotypes = %v, want 0|0 0|1 0|0", m101[9:12])
	}

	// The q arm's founder bits are monomorphic and absent, so position 70 carries de-novo
	// records with no founder record — the merge is not an annotation of existing sites.
	if _, ok := byID["F70"]; ok {
		t.Errorf("founder bit 70 is monomorphic and should not have been written")
	}
	// Fixed-in-sample mutations obey the same rule as fixed founder bits.
	if _, ok := byID["M105"]; ok {
		t.Errorf("M105 is carried by every sampled strand and should be dropped when "+
			"IncludeFixed is false; records: %d", len(records))
	}
	// Nothing was parked at position 0 (the pre-refcount failure mode for an unpooled id).
	if rec, ok := byID["M104"]; ok {
		t.Errorf("unpooled mutation 104 was written anyway, at POS %s", rec[1])
	}
}

// IncludeFixed reaches the pool too: the fixed mutation and the monomorphic q arm both
// come back, and the out-of-range / unpooled drops are NOT overridden by it — they have
// no coordinate, which is a different thing from being uninformative.
func TestExportVCFMergedIncludeFixed(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportVCF(model, pop, []int{1, 2, 3},
		VCFOptions{IncludeFixed: true, IncludeMutations: true})
	if err != nil {
		t.Fatalf("ExportVCF failed: %v", err)
	}
	if stats.NumFounderSites != 128 {
		t.Errorf("founder sites = %d, want 128", stats.NumFounderSites)
	}
	if stats.NumMutationSites != 4 {
		t.Errorf("de-novo sites = %d, want 4 (100,101,102,105 — 103 is out of range and "+
			"104 unpooled, and neither has a position to write)", stats.NumMutationSites)
	}
}

// The default path cannot see the pool at all. Exporting a pop that carries mutations
// with the zero VCFOptions must produce the same bytes as exporting the same pop with
// the pool removed — which is the byte-identity discipline every TMR4A item is held to,
// applied to the one file this item changes.
func TestExportVCFDefaultPathIgnoresPool(t *testing.T) {
	withPool, popA := mergedFixture()
	withPool.ResultsDir = t.TempDir()
	statsA, err := ExportVCF(withPool, popA, []int{1, 2, 3, 4}, VCFOptions{})
	if err != nil {
		t.Fatalf("export with pool: %v", err)
	}

	noPool, popB := mergedFixture()
	noPool.ResultsDir = t.TempDir()
	popB.IndMutations, popB.MutationPool = nil, nil
	statsB, err := ExportVCF(noPool, popB, []int{1, 2, 3, 4}, VCFOptions{})
	if err != nil {
		t.Fatalf("export without pool: %v", err)
	}

	a, err := os.ReadFile(statsA.Path)
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	b, err := os.ReadFile(statsB.Path)
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	if string(a) != string(b) {
		t.Errorf("the de-novo pool leaked into a default (non-merged) export")
	}
	if statsA.NumMutationSites != 0 || statsA.MutationsIndexed != 0 {
		t.Errorf("merged counters are non-zero on the default path: %d sites, %d indexed",
			statsA.NumMutationSites, statsA.MutationsIndexed)
	}
	if statsA.NumSites != 64 || statsB.NumSites != 64 {
		t.Errorf("sites = %d / %d, want 64 / 64", statsA.NumSites, statsB.NumSites)
	}
}

// labelFixture is a two-generation pedigree with strand provenance and W3 labels, plus
// just enough bitfield to make the focal locus a segregating site worth writing.
//
//	founders 1 (labels 0,1) and 2 (labels 2,3) under A=4
//	child 3 takes dad strand 0 (label 0) and mom strand 0 (label 2)
//	child 4 takes dad strand 1 (label 1) and mom strand 0 (label 2)
func labelFixture(t *testing.T) (*core.Model, *core.Pop) {
	t.Helper()
	model := makeIBDTestModel()
	model.StringParams = map[string]string{"linkage_model": "arm"}
	model.ARGLoci = []int{10}
	model.ResultsDir = t.TempDir()

	b := newTMRKA([]int{10})
	b.founder(1, -30)
	b.founder(2, -30)
	b.pop.AssignFounderLabels([]int{1, 2}, 4)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 5, []int{1}, []int{0})
	b.alive(3, 4)

	// Only ind 3 carries the derived allele at bit 10, so the site segregates.
	b.pop.Chromosomes = map[int][][]uint64{
		3: chrom(1<<10, 0, 0, 0),
		4: chrom(0, 0, 0, 0),
	}
	return model, b.pop
}

// The truth annotation: FORMAT FL carries, per sample per phased strand, the created
// allele that strand actually descends from — resolved by walking the recorded
// genealogy, not inferred from the sequence. This is the join key W7 needs.
func TestExportVCFAlleleLabelsAnnotateFocalLoci(t *testing.T) {
	model, pop := labelFixture(t)

	stats, err := ExportVCF(model, pop, []int{3, 4}, VCFOptions{AlleleLabels: true})
	if err != nil {
		t.Fatalf("ExportVCF failed: %v", err)
	}
	header, records := readVCF(t, stats.Path)
	if !containsLine(header, "##FORMAT=<ID=FL,Number=2,Type=Integer,Description=\"Ground-truth created-allele label (TMR4A W3) of each phased strand at this focal locus; . = not recoverable\">") {
		t.Errorf("labelled export must declare the FL FORMAT field:\n%s", strings.Join(header, "\n"))
	}

	var focal []string
	for _, r := range records {
		f := strings.Split(r, "\t")
		if f[1] == "11" { // bit 10 -> POS 11
			focal = f
		} else if f[8] != "GT" {
			t.Errorf("non-focal record at POS %s carries FORMAT %q; FL belongs only at "+
				"focal loci", f[1], f[8])
		}
	}
	if focal == nil {
		t.Fatalf("no record at the focal locus (POS 11); got %d records", len(records))
	}
	if focal[8] != "GT:FL" {
		t.Fatalf("focal FORMAT = %q, want GT:FL", focal[8])
	}
	// ind 3: strand 0 -> founder 1 strand 0 -> label 0; strand 1 -> founder 2 strand 0 -> label 2.
	// ind 4: strand 0 -> founder 1 strand 1 -> label 1; strand 1 -> founder 2 strand 0 -> label 2.
	if focal[9] != "1|0:0,2" {
		t.Errorf("IND3 field = %q, want 1|0:0,2", focal[9])
	}
	if focal[10] != "0|0:1,2" {
		t.Errorf("IND4 field = %q, want 0|0:1,2", focal[10])
	}
}

// The annotation and the walker must agree, because they are two readings of the same
// genealogy: when no source is unlabelled, the distinct FL values across the sample are
// exactly the walker's CreatedAllelesSurviving. This is the cross-check that keeps the
// per-lineage resolver in tmrka_sources.go from drifting away from walkLocus.
func TestFounderAlleleLabelsAgreeWithWalker(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, -50)
	b.founder(2, -50)
	b.founder(3, -50)
	b.founder(4, -50)
	b.pop.AssignFounderLabels([]int{1, 2, 3, 4}, 2)
	b.child(5, 1, 2, 0, []int{0}, []int{0})
	b.child(6, 1, 2, 0, []int{0}, []int{0})
	b.child(7, 3, 4, 0, []int{1}, []int{0})
	b.alive(5, 6, 7)

	ids := []int{5, 6, 7}
	loc := ComputeTMRKA(b.model(4, 60), b.pop, ids).Loci[0]
	if loc.UnlabelledSources != 0 {
		t.Fatalf("fixture has %d unlabelled sources; the identity below only holds at 0",
			loc.UnlabelledSources)
	}

	got := FounderAlleleLabelsAt(b.pop, ids, 0, 1)
	distinct := map[int]bool{}
	for _, id := range ids {
		distinct[got[id][0]] = true
		distinct[got[id][1]] = true
	}
	if len(distinct) != loc.CreatedAllelesSurviving {
		t.Errorf("per-strand resolver sees %d distinct created alleles, walker reports %d: %v",
			len(distinct), loc.CreatedAllelesSurviving, got)
	}
}

// An unresolvable strand reports -1 (rendered '.' in the VCF), never a shared default.
// Two '.' strands may be two different unknown alleles, and the file must not claim
// otherwise — the same rule the walker's UnlabelledSources enforces.
func TestFounderAlleleLabelsReportUnknownAsUnknown(t *testing.T) {
	b := newTMRKA([]int{0})
	b.founder(1, -30)
	b.founder(2, -30)
	b.child(3, 1, 2, 0, []int{0}, []int{0})
	b.child(4, 1, 2, 5, []int{1}, []int{1})
	b.alive(3, 4)
	b.pop.AssignFounderLabels([]int{1}, 1) // founder 2 is left unlabelled

	got := FounderAlleleLabelsAt(b.pop, []int{3, 4}, 0, 1)
	if got[3][0] != 0 || got[4][0] != 0 {
		t.Errorf("paternal strands: got %d/%d, want 0/0 (both trace to founder 1, A=1)",
			got[3][0], got[4][0])
	}
	if got[3][1] != -1 || got[4][1] != -1 {
		t.Errorf("maternal strands: got %d/%d, want -1/-1 (founder 2 carries no label)",
			got[3][1], got[4][1])
	}
	if labelField(got[3][1]) != "." {
		t.Errorf("unknown label rendered as %q, want \".\"", labelField(got[3][1]))
	}
}

// indexRecords keys the emitted records by their VCF ID column, which in merged mode is
// the identity (F<bit> / M<mutid>) rather than the coordinate — several records share a
// POS by construction, so POS cannot be the key.
func indexRecords(t *testing.T, records []string) map[string][]string {
	t.Helper()
	out := make(map[string][]string, len(records))
	for _, r := range records {
		f := strings.Split(r, "\t")
		if len(f) < 10 {
			t.Fatalf("malformed record: %q", r)
		}
		if _, dup := out[f[2]]; dup {
			t.Errorf("duplicate record ID %q — ids must identify a site uniquely", f[2])
		}
		out[f[2]] = f
	}
	return out
}
