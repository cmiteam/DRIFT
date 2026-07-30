package analysis

import (
	"bytes"
	"compress/gzip"
	"drift/pkg/core"
	"drift/pkg/individual"
	"os"
	"path/filepath"
	"testing"
)

// A hand-built fixture VCF. Sites (tab-separated in the raw string below):
//
//	chrom1 POS100 A>G AA=A  GT 1|0 0|0 0|0 0|0  -> AA==REF -> derived=ALT(1); id0 s0
//	chrom1 POS200 C>T AA=T  GT 0|0 0|0 1|1 1|1  -> AA==ALT -> FLIP derived=REF(0); id0,id1 both strands
//	chrom1 POS300 A>C .     GT 0|1 0|0 0|0 0|0  -> no AA   -> fallback derived=ALT(1); id0 s1
//	chrom2 POS100 A>G,T .   GT ...               -> multiallelic: SKIP
//	chrom2 POS200 AC>G .    GT ...               -> indel: SKIP
//	chrom2 POS300 A>G AA=.  GT 0|0 0|0 0|0 0|0  -> monomorphic (0 derived): SKIP
//	chrom2 POS400 A>G AA=.  GT 1|1 1|1 1|1 1|1  -> monomorphic (all derived): SKIP
//	chrom2 POS500 G>A AA=G  GT 1|0 0|1 1|1 0|0  -> AA==REF -> derived=ALT(1); id0 s0, id1 s1, id2 both
//
// Kept sites, contig-grouped file order: bit0=c1/100, bit1=c1/200, bit2=c1/300, bit3=c2/500.
const fixtureVCF = "##fileformat=VCFv4.2\n" +
	"##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n" +
	"1\t100\t.\tA\tG\t.\tPASS\tAA=A\tGT\t1|0\t0|0\t0|0\t0|0\n" +
	"1\t200\t.\tC\tT\t.\tPASS\tAA=T\tGT\t0|0\t0|0\t1|1\t1|1\n" +
	"1\t300\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\t0|0\n" +
	"2\t100\t.\tA\tG,T\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\t0|0\n" +
	"2\t200\t.\tAC\tG\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\t0|0\n" +
	"2\t300\t.\tA\tG\t.\tPASS\tAA=.\tGT\t0|0\t0|0\t0|0\t0|0\n" +
	"2\t400\t.\tA\tG\t.\tPASS\tAA=.\tGT\t1|1\t1|1\t1|1\t1|1\n" +
	"2\t500\t.\tG\tA\t.\tPASS\tAA=G\tGT\t1|0\t0|1\t1|1\t0|0\n"

func writeTemp(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0644); err != nil {
		t.Fatalf("write temp: %v", err)
	}
	return p
}

func newImportModel(dir string) (*core.Model, *core.Pop) {
	m := &core.Model{
		Parameters:     map[string]float64{},
		StringParams:   map[string]string{},
		FreeParameters: map[string]int{"run": 0, "year": 0},
		ModelName:      "Imp",
		ResultsDir:     dir,
	}
	p := &core.Pop{
		IndData:      map[int][]int{},
		Chromosomes:  map[int][][]uint64{},
		IndMutations: map[int]map[int][]int{},
		MutationPool: map[int]core.Mutation{},
		MutationHist: map[int]int{},
	}
	return m, p
}

func hasBit(pop *core.Pop, id, strand, bit int) bool {
	return bitAt(pop.Chromosomes[id][strand], bit)
}

func strandHasMut(pop *core.Pop, id, strand, mid int) bool {
	for _, m := range pop.IndMutations[id][strand] {
		if m == mid {
			return true
		}
	}
	return false
}

func TestImportVCFParse(t *testing.T) {
	path := writeTemp(t, "fix.vcf", fixtureVCF)
	m, p := newImportModel(t.TempDir())
	st, err := ImportVCF(path, m, p, DefaultImportPolicy())
	if err != nil {
		t.Fatalf("ImportVCF: %v", err)
	}

	// --- stats ---
	if st.NumSamples != 4 || st.NumSites != 4 || st.GenomeBits != 4 || st.NumContigs != 2 {
		t.Fatalf("stats: samples=%d sites=%d bits=%d contigs=%d (want 4/4/4/2)",
			st.NumSamples, st.NumSites, st.GenomeBits, st.NumContigs)
	}
	if st.SkipMulti != 1 || st.SkipIndel != 1 || st.SkipMonomorph != 2 {
		t.Fatalf("skips: multi=%d indel=%d mono=%d (want 1/1/2)", st.SkipMulti, st.SkipIndel, st.SkipMonomorph)
	}
	// polarization: kept AA sites = POS100,POS200,POS500 (3); one flip (POS200);
	// fallbacks = chrom1 POS300 + the two AA=. mono sites (3).
	if st.PolarizedAA != 3 || st.PolarizedFlip != 1 || st.PolarizedFall != 3 {
		t.Fatalf("polarization: AA=%d flip=%d fall=%d (want 3/1/3)", st.PolarizedAA, st.PolarizedFlip, st.PolarizedFall)
	}

	// --- model patched ---
	if m.FreeParameters["genome_bits"] != 4 {
		t.Fatalf("genome_bits=%d want 4", m.FreeParameters["genome_bits"])
	}
	if got := m.ChromosomeArms[1][0]; got[0] != 0 || got[1] != 3 {
		t.Fatalf("chrom1 arm = %v want [0 3]", got)
	}
	if got := m.ChromosomeArms[2][0]; got[0] != 3 || got[1] != 1 {
		t.Fatalf("chrom2 arm = %v want [3 1]", got)
	}

	// --- bitfield contents (hand-computed) ---
	type bset struct{ id, strand, bit int }
	wantSet := []bset{
		{0, 0, 0},                                  // bit0: id0 s0
		{0, 0, 1}, {0, 1, 1}, {1, 0, 1}, {1, 1, 1}, // bit1: id0,id1 both strands
		{0, 1, 2},                                  // bit2: id0 s1
		{0, 0, 3}, {1, 1, 3}, {2, 0, 3}, {2, 1, 3}, // bit3: id0 s0, id1 s1, id2 both
	}
	setMap := map[bset]bool{}
	for _, b := range wantSet {
		setMap[b] = true
	}
	for id := 0; id < 4; id++ {
		for s := 0; s < 2; s++ {
			for bit := 0; bit < 4; bit++ {
				got := hasBit(p, id, s, bit)
				want := setMap[bset{id, s, bit}]
				if got != want {
					t.Errorf("bitfield id%d s%d bit%d = %v want %v", id, s, bit, got, want)
				}
			}
		}
	}

	// --- pool consistency: mid = bit+1, Position = bit; strand carries mid iff bit set ---
	for _, b := range wantSet {
		mid := b.bit + 1
		if !strandHasMut(p, b.id, b.strand, mid) {
			t.Errorf("pool: id%d s%d missing mut %d (bit %d)", b.id, b.strand, mid, b.bit)
		}
		if mut, ok := p.MutationPool[mid]; !ok || mut.Position != b.bit {
			t.Errorf("pool: mut %d Position=%d want %d", mid, mut.Position, b.bit)
		}
	}
	// and the pool must NOT carry a mutation where the bit is clear
	for id := 0; id < 4; id++ {
		for s := 0; s < 2; s++ {
			for bit := 0; bit < 4; bit++ {
				mid := bit + 1
				if strandHasMut(p, id, s, mid) != hasBit(p, id, s, bit) {
					t.Errorf("INCONSISTENT id%d s%d bit%d: pool=%v bitfield=%v",
						id, s, bit, strandHasMut(p, id, s, mid), hasBit(p, id, s, bit))
				}
			}
		}
	}
}

// TestImportVCFPolarizeRef: with policy "ref", AA is ignored so POS200 does NOT flip
// (derived=ALT), giving a different bitfield than the AA path — a distinguishable
// outcome driven purely by the polarization knob.
func TestImportVCFPolarizeRef(t *testing.T) {
	path := writeTemp(t, "fix.vcf", fixtureVCF)
	m, p := newImportModel(t.TempDir())
	pol := DefaultImportPolicy()
	pol.Polarize = "ref"
	if _, err := ImportVCF(path, m, p, pol); err != nil {
		t.Fatalf("ImportVCF: %v", err)
	}
	// bit1 = POS200 (C>T), REF-ancestral => derived=ALT(1). GT 0|0 0|0 1|1 1|1 =>
	// derived carried by id2,id3 (both strands), NOT id0/id1 (opposite of the AA/flip case).
	if hasBit(p, 0, 0, 1) || !hasBit(p, 2, 0, 1) {
		t.Fatalf("policy=ref bit1: id0=%v id2=%v (want id0 clear, id2 set — no flip)",
			hasBit(p, 0, 0, 1), hasBit(p, 2, 0, 1))
	}
}

// TestImportVCFDifferentiated builds a 2-deme, all-sites-fixed-different panel and
// checks the pipeline recovers Fst = 1 (bitfield reader) AND segregating sites from
// the pool (pool reader) — the both-structures payoff on a distinguishable outcome.
func TestImportVCFDifferentiated(t *testing.T) {
	// 4 samples; deme A = {S1,S2}, deme B = {S3,S4}. 3 SNP sites, each fixed derived
	// in A and fixed ancestral in B (with REF-ancestral, that's GT 1|1 in A, 0|0 in B).
	vcf := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n" +
		"1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\t1|1\t0|0\t0|0\n" +
		"1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\t1|1\t0|0\t0|0\n" +
		"1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\t1|1\t0|0\t0|0\n"
	path := writeTemp(t, "diff.vcf", vcf)
	m, p := newImportModel(t.TempDir())
	pol := DefaultImportPolicy()
	pol.PopMap = map[string]int{"S1": 5, "S2": 5, "S3": 9, "S4": 9} // arbitrary ids -> remapped 0,1
	st, err := ImportVCF(path, m, p, pol)
	if err != nil {
		t.Fatalf("ImportVCF: %v", err)
	}
	if st.NumDemes != 2 {
		t.Fatalf("NumDemes=%d want 2", st.NumDemes)
	}
	ids := []int{0, 1, 2, 3}
	fst := ComputeFst(m, p, ids)
	if len(fst.Pairs) != 1 {
		t.Fatalf("Fst pairs=%d want 1", len(fst.Pairs))
	}
	if fst.Pairs[0].Hudson < 0.999 {
		t.Fatalf("Hudson Fst=%.4f want ~1 (fully differentiated)", fst.Pairs[0].Hudson)
	}
	// pool reader sees all 3 sites segregating
	ns := ComputeNeutralStats(p, ids)
	if ns.SegSites != 3 {
		t.Fatalf("pool SegSites=%d want 3", ns.SegSites)
	}
}

// TestImportVCFRoundTrip: import -> ExportVCF -> re-import yields an identical
// bitfield and genome_bits, proving import is the faithful inverse of export.
func TestImportVCFRoundTrip(t *testing.T) {
	path := writeTemp(t, "fix.vcf", fixtureVCF)
	dir := t.TempDir()
	m, p := newImportModel(dir)
	if _, err := ImportVCF(path, m, p, DefaultImportPolicy()); err != nil {
		t.Fatalf("import: %v", err)
	}
	ids := []int{0, 1, 2, 3}
	vcfStats, err := ExportVCF(m, p, ids, false)
	if err != nil {
		t.Fatalf("export: %v", err)
	}
	// re-import the exported file (REF=A/ALT=T, no AA -> REF-ancestral fallback).
	m2, p2 := newImportModel(t.TempDir())
	if _, err := ImportVCF(vcfStats.Path, m2, p2, DefaultImportPolicy()); err != nil {
		t.Fatalf("re-import: %v", err)
	}
	if m2.FreeParameters["genome_bits"] != m.FreeParameters["genome_bits"] {
		t.Fatalf("genome_bits %d != %d", m2.FreeParameters["genome_bits"], m.FreeParameters["genome_bits"])
	}
	bits := m.FreeParameters["genome_bits"]
	for id := 0; id < 4; id++ {
		for s := 0; s < 2; s++ {
			for bit := 0; bit < bits; bit++ {
				if hasBit(p, id, s, bit) != hasBit(p2, id, s, bit) {
					t.Errorf("round-trip mismatch id%d s%d bit%d", id, s, bit)
				}
			}
		}
	}
}

func TestImportVCFGzip(t *testing.T) {
	var buf bytes.Buffer
	gz := gzip.NewWriter(&buf)
	if _, err := gz.Write([]byte(fixtureVCF)); err != nil {
		t.Fatal(err)
	}
	gz.Close()
	p := filepath.Join(t.TempDir(), "fix.vcf.gz")
	if err := os.WriteFile(p, buf.Bytes(), 0644); err != nil {
		t.Fatal(err)
	}
	m, pop := newImportModel(t.TempDir())
	st, err := ImportVCF(p, m, pop, DefaultImportPolicy())
	if err != nil {
		t.Fatalf("gzip import: %v", err)
	}
	if st.NumSites != 4 {
		t.Fatalf("gzip NumSites=%d want 4", st.NumSites)
	}
}

func TestLoadPanelFile(t *testing.T) {
	panel := "sample\tpop\tsuper_pop\tgender\n" +
		"S1\tYRI\tAFR\tmale\n" +
		"S2\tCEU\tEUR\tfemale\n" +
		"S3\tYRI\tAFR\tmale\n"
	path := writeTemp(t, "p.panel", panel)
	m, names, err := LoadPanelFile(path, 1) // col 1 = pop
	if err != nil {
		t.Fatalf("LoadPanelFile: %v", err)
	}
	if len(names) != 2 || names[0] != "YRI" || names[1] != "CEU" {
		t.Fatalf("populations = %v want [YRI CEU]", names)
	}
	if m["S1"] != 0 || m["S2"] != 1 || m["S3"] != 0 {
		t.Fatalf("deme map = %v want S1/S3=0 S2=1", m)
	}
	// super_pop column
	m2, names2, err := LoadPanelFile(path, 2)
	if err != nil {
		t.Fatalf("LoadPanelFile col2: %v", err)
	}
	if len(names2) != 2 || m2["S1"] != m2["S3"] {
		t.Fatalf("super_pop map = %v names=%v", m2, names2)
	}
}

// Guard: AlleleCount is populated so downstream reporting is consistent.
func TestImportVCFAlleleCount(t *testing.T) {
	path := writeTemp(t, "fix.vcf", fixtureVCF)
	m, p := newImportModel(t.TempDir())
	if _, err := ImportVCF(path, m, p, DefaultImportPolicy()); err != nil {
		t.Fatalf("import: %v", err)
	}
	// id0 set bits: bit0 s0, bit1 s0+s1, bit2 s1, bit3 s0 = 5 total.
	if got := p.IndData[0][individual.AlleleCount]; got != 5 {
		t.Fatalf("id0 AlleleCount=%d want 5", got)
	}
}
