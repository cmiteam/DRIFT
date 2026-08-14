package analysis

import (
	"bufio"
	"os"
	"strconv"
	"strings"
	"testing"
)

// ARGweaver export (TMR4A.md W7). The load-bearing properties are format fidelity —
// ARGweaver requires sorted, distinct coordinates and one base per haploid genome — and
// agreement with the W6 VCF, since the study is a join between what a tool infers from
// this file and what the W5 walker knows about the same sample.

// readSites parses a .sites file back into its header lines and (position, alleles) rows.
func readSites(t *testing.T, path string) (names []string, region []string, pos []int, alleles []string) {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open sites: %v", err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	for sc.Scan() {
		fields := strings.Split(sc.Text(), "\t")
		switch fields[0] {
		case "NAMES":
			names = fields[1:]
		case "REGION":
			region = fields[1:]
		default:
			p, err := strconv.Atoi(fields[0])
			if err != nil {
				t.Fatalf("unparseable site line %q", sc.Text())
			}
			if len(fields) != 2 {
				t.Fatalf("site line must have exactly 2 tab-separated columns: %q", sc.Text())
			}
			pos = append(pos, p)
			alleles = append(alleles, fields[1])
		}
	}
	if err := sc.Err(); err != nil {
		t.Fatalf("scan sites: %v", err)
	}
	return names, region, pos, alleles
}

// The format contract, asserted directly against the manual's requirements: two names per
// individual, a 1-based end-inclusive REGION, sorted and DISTINCT positions, and one
// allele character per haploid genome on every row.
func TestARGweaverSitesFormat(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3, 4}, ARGweaverOptions{BpPerBit: 10})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	if stats.NumHaplotypes != 6 {
		t.Fatalf("haplotypes = %d, want 6 (3 individuals with chromosomes x 2 strands)",
			stats.NumHaplotypes)
	}
	if len(stats.Regions) != 1 {
		t.Fatalf("regions = %d, want 1", len(stats.Regions))
	}

	names, region, pos, alleles := readSites(t, stats.Regions[0].Path)

	want := []string{"IND1_1", "IND1_2", "IND2_1", "IND2_2", "IND3_1", "IND3_2"}
	if strings.Join(names, ",") != strings.Join(want, ",") {
		t.Errorf("NAMES = %v, want %v", names, want)
	}
	// The fixture chromosome spans 128 bits; at 10 bp/bit that is 1280 bp, written
	// 1-based end-inclusive.
	if len(region) != 3 || region[0] != "0" || region[1] != "1" || region[2] != "1280" {
		t.Errorf("REGION = %v, want [0 1 1280]", region)
	}

	if len(pos) == 0 {
		t.Fatal("no sites written")
	}
	for i := range pos {
		if i > 0 && pos[i] <= pos[i-1] {
			t.Fatalf("positions must be sorted and distinct: pos[%d]=%d follows %d",
				i, pos[i], pos[i-1])
		}
		if len(alleles[i]) != stats.NumHaplotypes {
			t.Fatalf("site at %d has %d allele chars, want %d",
				pos[i], len(alleles[i]), stats.NumHaplotypes)
		}
		if strings.Trim(alleles[i], "AT") != "" {
			t.Fatalf("site at %d has non-ACGT alleles %q", pos[i], alleles[i])
		}
		// A written site must be genuinely segregating: the format takes every
		// unlisted position as invariant, so an invariant listed row is dead weight.
		if !strings.Contains(alleles[i], "A") || !strings.Contains(alleles[i], "T") {
			t.Errorf("site at %d is invariant (%q) and should not have been written",
				pos[i], alleles[i])
		}
	}
}

// Colliding sites — the founder bit and the de-novo mutations sharing one genome bit —
// must land on DISTINCT bp inside that bit's window, because .sites keys on position.
// This is the case §8c measured at 213/223 in a real run, so getting it wrong would
// silently discard most of the mutational history.
func TestARGweaverSpreadsCollidingSitesWithinTheBitWindow(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3}, ARGweaverOptions{BpPerBit: 10})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	_, _, pos, _ := readSites(t, stats.Regions[0].Path)

	// Bit 5 carries the founder site AND mutation 102. At 10 bp/bit the bit's window is
	// [51, 60], so both must be written there, at 51 and 52.
	found := map[int]bool{}
	for _, p := range pos {
		found[p] = true
	}
	if !found[51] || !found[52] {
		t.Errorf("bit 5's two colliding sites should occupy 51 and 52; positions present near "+
			"the window: %v", nearby(pos, 45, 65))
	}
	// Bit 70 carries mutations 100 and 101 but no founder site (the q arm is
	// monomorphic), so its window [701, 710] holds exactly two sites.
	if !found[701] || !found[702] {
		t.Errorf("bit 70's two de-novo sites should occupy 701 and 702; got %v",
			nearby(pos, 695, 715))
	}
	if found[703] {
		t.Errorf("bit 70 should hold exactly 2 sites, but 703 is occupied too")
	}
	if stats.TotalDropped != 0 {
		t.Errorf("dropped %d sites at 10 bp/bit, want 0", stats.TotalDropped)
	}
}

// With no room in the window a site is dropped and COUNTED, never written at a coordinate
// that belongs to a different genome bit — which would move a site to another locus.
func TestARGweaverDropsAndCountsWindowOverflow(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	// bp_per_bit = 1 leaves each bit room for exactly one site.
	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3}, ARGweaverOptions{BpPerBit: 1})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	// Bit 5 (founder + M102) loses one; bit 70 (M100 + M101) loses one.
	if stats.TotalDropped != 2 {
		t.Errorf("TotalDropped = %d, want 2", stats.TotalDropped)
	}
	_, _, pos, _ := readSites(t, stats.Regions[0].Path)
	for i := 1; i < len(pos); i++ {
		if pos[i] <= pos[i-1] {
			t.Fatalf("positions still must be sorted and distinct after dropping: %d then %d",
				pos[i-1], pos[i])
		}
	}
}

// The .sites file and the W6 VCF must describe the SAME sites in the SAME order. They are
// two renderings of one site list and the study joins across them, so a divergence would
// corrupt the truth-vs-inference comparison invisibly.
func TestARGweaverSitesMatchVCFRecords(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()
	ids := []int{1, 2, 3}

	vcfStats, err := ExportVCF(model, pop, ids, VCFOptions{IncludeMutations: true})
	if err != nil {
		t.Fatalf("ExportVCF failed: %v", err)
	}
	argStats, err := ExportARGweaver(model, pop, ids, ARGweaverOptions{BpPerBit: 10})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}

	if argStats.TotalSites != vcfStats.NumSites {
		t.Errorf("sites disagree: .sites has %d, VCF has %d (%d founder + %d de-novo)",
			argStats.TotalSites, vcfStats.NumSites, vcfStats.NumFounderSites, vcfStats.NumMutationSites)
	}

	// Same genotypes, site for site: the VCF's phased "a|b" per sample must equal the
	// .sites row's two characters for that sample's two haplotypes.
	_, records := readVCF(t, vcfStats.Path)
	_, _, _, alleles := readSites(t, argStats.Regions[0].Path)
	if len(records) != len(alleles) {
		t.Fatalf("record counts differ: %d VCF vs %d sites", len(records), len(alleles))
	}
	for i, rec := range records {
		f := strings.Split(rec, "\t")
		var got strings.Builder
		for _, gt := range f[9:] {
			for _, half := range strings.Split(gt, "|") {
				if half == "1" {
					got.WriteByte('T')
				} else {
					got.WriteByte('A')
				}
			}
		}
		if got.String() != alleles[i] {
			t.Fatalf("site %d (%s): VCF genotypes render as %q, .sites says %q",
				i, f[2], got.String(), alleles[i])
		}
	}
}

// argweaver_chroms restricts the export, which is how a real study scopes to one arm
// (§1A: one chromosome arm is larger than ARGweaver's native ~2 Mb block).
func TestARGweaverChromSelection(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	if _, err := ExportARGweaver(model, pop, []int{1, 2, 3},
		ARGweaverOptions{BpPerBit: 10, Chroms: []int{99}}); err == nil {
		t.Error("selecting a chromosome that does not exist should be an error, not an empty export")
	}
	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3},
		ARGweaverOptions{BpPerBit: 10, Chroms: []int{0}})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	if len(stats.Regions) != 1 || stats.Regions[0].Chrom != 0 {
		t.Errorf("regions = %+v, want just chromosome 0", stats.Regions)
	}
}

// The region BED is 0-based half-open where the REGION line is 1-based end-inclusive.
// Both describe the same span; getting the convention wrong shifts every downstream mask.
func TestARGweaverRegionBedConvention(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3}, ARGweaverOptions{BpPerBit: 10})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	data, err := os.ReadFile(stats.BedPath)
	if err != nil {
		t.Fatalf("read BED: %v", err)
	}
	line := strings.TrimSpace(string(data))
	if line != "0\t0\t1280" {
		t.Errorf("BED = %q, want \"0\\t0\\t1280\" (0-based half-open for the 1..1280 REGION)", line)
	}
}

func nearby(pos []int, lo, hi int) []int {
	var out []int
	for _, p := range pos {
		if p >= lo && p <= hi {
			out = append(out, p)
		}
	}
	return out
}
