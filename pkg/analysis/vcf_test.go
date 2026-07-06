package analysis

import (
	"bufio"
	"os"
	"strings"
	"testing"
)

// TestExportVCF_SegregatingOnly uses the shared IBD test fixtures (one 128-bit
// chromosome, p arm [0,64), q arm [64,128)). Samples 1 & 2 carry the whole p
// arm on copy 0; sample 3 carries only its low 32 bits; sample 4 has no
// chromosomes and must be dropped.
func TestExportVCF_SegregatingOnly(t *testing.T) {
	model := makeIBDTestModel()
	pop := makeIBDTestPop()
	model.ResultsDir = t.TempDir()

	stats, err := ExportVCF(model, pop, []int{1, 2, 3, 4}, false)
	if err != nil {
		t.Fatalf("ExportVCF failed: %v", err)
	}
	// Ind 4 (no chromosomes) is skipped -> 3 samples.
	if stats.NumSamples != 3 {
		t.Errorf("NumSamples = %d, want 3", stats.NumSamples)
	}
	// Every p-arm bit (0..63) is polymorphic within the sample; the q arm is all
	// zero (monomorphic) and excluded when includeFixed is false.
	if stats.NumSites != 64 {
		t.Errorf("NumSites = %d, want 64", stats.NumSites)
	}

	header, records := readVCF(t, stats.Path)

	// Contig length is the chromosome span (128 bits), not genome-global.
	if !containsLine(header, "##contig=<ID=0,length=128>") {
		t.Errorf("missing/incorrect ##contig header:\n%s", strings.Join(header, "\n"))
	}
	// #CHROM line carries exactly the retained samples, in id order.
	wantCols := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND1\tIND2\tIND3"
	if !containsLine(header, wantCols) {
		t.Errorf("column header mismatch; want %q in:\n%s", wantCols, strings.Join(header, "\n"))
	}

	// Bit 0: all three samples carry the derived allele on copy 0 -> 1|0 each,
	// AC=3. POS = bit - base(0) + 1 = 1.
	first := strings.Split(records[0], "\t")
	if first[0] != "0" || first[1] != "1" {
		t.Errorf("first record CHROM/POS = %s/%s, want 0/1", first[0], first[1])
	}
	if first[7] != "AC=3;AN=6" {
		t.Errorf("first record INFO = %s, want AC=3;AN=6", first[7])
	}
	if first[9] != "1|0" || first[10] != "1|0" || first[11] != "1|0" {
		t.Errorf("first record genotypes = %v, want 1|0 1|0 1|0", first[9:12])
	}

	// Bit 32: sample 3 lost its derived allele -> AC=2 (samples 1 & 2 only).
	rec32 := strings.Split(records[32], "\t")
	if rec32[1] != "33" {
		t.Errorf("record[32] POS = %s, want 33", rec32[1])
	}
	if rec32[7] != "AC=2;AN=6" {
		t.Errorf("record[32] INFO = %s, want AC=2;AN=6", rec32[7])
	}
	if rec32[11] != "0|0" {
		t.Errorf("record[32] IND3 genotype = %s, want 0|0", rec32[11])
	}
}

// TestExportVCF_IncludeFixed emits every position, so the monomorphic q arm is
// present too: 128 sites total.
func TestExportVCF_IncludeFixed(t *testing.T) {
	model := makeIBDTestModel()
	pop := makeIBDTestPop()
	model.ResultsDir = t.TempDir()

	stats, err := ExportVCF(model, pop, []int{1, 2, 3}, true)
	if err != nil {
		t.Fatalf("ExportVCF failed: %v", err)
	}
	if stats.NumSites != 128 {
		t.Errorf("NumSites = %d, want 128 (all positions)", stats.NumSites)
	}
}

func readVCF(t *testing.T, path string) (header, records []string) {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open VCF: %v", err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "#") {
			header = append(header, line)
		} else {
			records = append(records, line)
		}
	}
	if err := sc.Err(); err != nil {
		t.Fatalf("scan VCF: %v", err)
	}
	return header, records
}

func containsLine(lines []string, want string) bool {
	for _, l := range lines {
		if l == want {
			return true
		}
	}
	return false
}
