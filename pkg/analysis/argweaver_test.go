package analysis

import (
	"bufio"
	"drift/pkg/core"
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

// argweaver_emit can only ever REDUCE output. Every caller that predates the parameter
// builds ARGweaverOptions with a zero Emit, and must still get all four artifacts.
func TestARGweaverEmitZeroValueWritesEverything(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3}, ARGweaverOptions{BpPerBit: 10})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	for _, p := range []struct{ name, path string }{
		{"sites", stats.Regions[0].Path},
		{"bed", stats.BedPath},
		{"scaling", stats.RationalePath},
		{"cmd", stats.CommandPath},
	} {
		if p.path == "" {
			t.Errorf("%s path is empty; the zero Emit must resolve to all four", p.name)
			continue
		}
		if _, err := os.Stat(p.path); err != nil {
			t.Errorf("%s not written: %v", p.name, err)
		}
	}
}

// Suppressing an artifact must remove the FILE without disturbing the statistics the
// remaining artifacts report on — the site walk still has to run.
func TestARGweaverEmitSuppressesArtifacts(t *testing.T) {
	model, pop := mergedFixture()
	dir := t.TempDir()
	model.ResultsDir = dir

	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3},
		ARGweaverOptions{BpPerBit: 10, Emit: ARGweaverEmit{Sites: true}})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	if stats.BedPath != "" || stats.RationalePath != "" || stats.CommandPath != "" {
		t.Errorf("bed=%q scaling=%q cmd=%q, want all empty", stats.BedPath, stats.RationalePath, stats.CommandPath)
	}
	if _, err := os.Stat(stats.Regions[0].Path); err != nil {
		t.Errorf("sites was selected but not written: %v", err)
	}
	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatalf("readdir: %v", err)
	}
	for _, e := range entries {
		if !strings.HasSuffix(e.Name(), ".sites") {
			t.Errorf("unexpected file %q; only .sites was selected", e.Name())
		}
	}
}

// Dropping the .sites files must not cost the caller the region spans or site counts,
// because the scaling rationale is computed from them.
func TestARGweaverEmitScalingOnlyStillCountsSites(t *testing.T) {
	model, pop := mergedFixture()
	dir := t.TempDir()
	model.ResultsDir = dir

	full, err := ExportARGweaver(model, pop, []int{1, 2, 3}, ARGweaverOptions{BpPerBit: 10})
	if err != nil {
		t.Fatalf("full export failed: %v", err)
	}

	model.ResultsDir = t.TempDir()
	lean, err := ExportARGweaver(model, pop, []int{1, 2, 3},
		ARGweaverOptions{BpPerBit: 10, Emit: ARGweaverEmit{Scaling: true}})
	if err != nil {
		t.Fatalf("scaling-only export failed: %v", err)
	}
	if lean.TotalSites != full.TotalSites || lean.TotalDropped != full.TotalDropped {
		t.Errorf("sites=%d dropped=%d, want %d/%d — the walk must run even when discarded",
			lean.TotalSites, lean.TotalDropped, full.TotalSites, full.TotalDropped)
	}
	if len(lean.Regions) != len(full.Regions) || lean.Regions[0].EndBp != full.Regions[0].EndBp {
		t.Errorf("region spans differ between full and scaling-only exports")
	}
	if _, err := os.Stat(lean.Regions[0].Path); err == nil {
		t.Errorf("a .sites file was written despite sites being unselected")
	}
}

func TestParseEmit(t *testing.T) {
	for _, tc := range []struct {
		spec string
		want ARGweaverEmit
	}{
		{"", allEmit()},
		{"   ", allEmit()},
		{"all", allEmit()},
		{"sites", ARGweaverEmit{Sites: true}},
		{"sites,cmd", ARGweaverEmit{Sites: true, Cmd: true}},
		{"SITES; Bed", ARGweaverEmit{Sites: true, Bed: true}},
		{"sites bed cmd scaling", allEmit()},
		{"sites,nonsense", ARGweaverEmit{Sites: true}},
	} {
		if got := parseEmit(tc.spec); got != tc.want {
			t.Errorf("parseEmit(%q) = %+v, want %+v", tc.spec, got, tc.want)
		}
	}
}

// A spec of nothing but garbage resolves to everything rather than silently exporting an
// empty run — the warning already told the user, and losing the data is the worse failure.
func TestARGweaverEmitAllUnknownResolvesToEverything(t *testing.T) {
	if got := parseEmit("nonsense").resolved(); got != allEmit() {
		t.Errorf("resolved = %+v, want everything", got)
	}
}

// The sampling knobs are omitted when unset so the command file never claims a setting
// the user did not choose, and rendered verbatim when they are.
func TestARGweaverSamplingFlags(t *testing.T) {
	if got := samplingFlags(ARGweaverOptions{}); got != "" {
		t.Errorf("unset flags rendered %q, want empty", got)
	}
	got := samplingFlags(ARGweaverOptions{MaxTime: 2000, NTimes: 40, Seed: 42})
	want := "--maxtime 2000 --ntimes 40 --randseed 42"
	if got != want {
		t.Errorf("samplingFlags = %q, want %q", got, want)
	}
}

func TestARGweaverCommandFileCarriesSamplingFlags(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3},
		ARGweaverOptions{BpPerBit: 10, MaxTime: 2000, NTimes: 40, Seed: 42})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	data, err := os.ReadFile(stats.CommandPath)
	if err != nil {
		t.Fatalf("read command file: %v", err)
	}
	cmd := string(data)
	for _, want := range []string{"--maxtime 2000", "--ntimes 40", "--randseed 42"} {
		if !strings.Contains(cmd, want) {
			t.Errorf("command file missing %q", want)
		}
	}
	if strings.Contains(cmd, "argweaver_maxtime is UNSET") {
		t.Errorf("command file warns maxtime is unset when it was set")
	}
}

// Unset maxtime/seed must produce the warnings, because arg-sample's own defaults are a
// 200000-generation grid and a clock seed — both wrong for a DRIFT run.
func TestARGweaverCommandFileWarnsOnUnsetGrid(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3}, ARGweaverOptions{BpPerBit: 10})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	data, err := os.ReadFile(stats.CommandPath)
	if err != nil {
		t.Fatalf("read command file: %v", err)
	}
	cmd := string(data)
	for _, want := range []string{"argweaver_maxtime is UNSET", "argweaver_seed is UNSET"} {
		if !strings.Contains(cmd, want) {
			t.Errorf("command file missing warning %q", want)
		}
	}
	// The explanatory header names --maxtime, so only the runnable lines can be checked.
	for _, line := range strings.Split(cmd, "\n") {
		if strings.HasPrefix(strings.TrimSpace(line), "#") {
			continue
		}
		for _, flag := range []string{"--maxtime", "--ntimes", "--randseed"} {
			if strings.Contains(line, flag) {
				t.Errorf("command line %q carries %s, which was never set", line, flag)
			}
		}
	}
}

// The command file must not send the reader to a scaling rationale argweaver_emit
// suppressed — those warnings are what decide whether the rates are trustworthy.
func TestARGweaverCommandFileDoesNotCiteSuppressedScaling(t *testing.T) {
	model, pop := mergedFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportARGweaver(model, pop, []int{1, 2, 3},
		ARGweaverOptions{BpPerBit: 10, Emit: ARGweaverEmit{Sites: true, Cmd: true}})
	if err != nil {
		t.Fatalf("ExportARGweaver failed: %v", err)
	}
	data, err := os.ReadFile(stats.CommandPath)
	if err != nil {
		t.Fatalf("read command file: %v", err)
	}
	cmd := string(data)
	if strings.Contains(cmd, "_argweaver_scaling.txt") {
		t.Errorf("command file cites a scaling file that was never written")
	}
	if !strings.Contains(cmd, "argweaver_emit suppressed") {
		t.Errorf("command file does not say the rationale was suppressed")
	}
}

func TestHistoryGenerations(t *testing.T) {
	model, _ := mergedFixture()
	model.Parameters["generation_time"] = 25
	model.Parameters["seed_year"] = 0
	model.FreeParameters["year"] = 500
	if got := historyGenerations(model); got != 20 {
		t.Errorf("historyGenerations = %v, want 20 (500 years / 25)", got)
	}

	model.Parameters["generation_time"] = 0
	if got := historyGenerations(model); got != 0 {
		t.Errorf("generation_time 0 must yield 0, got %v", got)
	}
	if got := historyGenerations(nil); got != 0 {
		t.Errorf("nil model must yield 0, got %v", got)
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

// poolOnlyFixture is the merged fixture with the BITFIELD REMOVED: a full de-novo
// mutation pool and no pop.Chromosomes entry for anyone. This is not a contrived
// state — it is what a mainstream equilibrium run produces. SeedingPopulation returns
// without allocating anything at init_heterozygosity = 0, and meiosis only allocates a
// child's chromosomes when a parent had set bits or presented a created germ cell
// (pkg/simulation/birth.go:123), so on a pool-only model NOBODY has a bitfield while
// the pool fills up normally.
func poolOnlyFixture() (*core.Model, *core.Pop) {
	model, pop := mergedFixture()
	pop.Chromosomes = map[int][][]uint64{}
	return model, pop
}

// A pool-only sample must export, and the two artifacts must agree on the site count.
//
// Before this was fixed, both exports required an allocated bitfield to admit an
// individual at all: ExportARGweaver failed outright with "no sampled individual has
// tracked chromosomes" and ExportVCF wrote a sample-less header — on a run whose pool
// held every site they needed. That made the mainstream arm of the ARGweaver study
// unexportable.
func TestPoolOnlyExportAgreesWithVCF(t *testing.T) {
	ids := []int{1, 2, 3, 4}

	// --- VCF ---
	model, pop := poolOnlyFixture()
	model.ResultsDir = t.TempDir()
	vcfStats, err := ExportVCF(model, pop, ids, VCFOptions{IncludeMutations: true})
	if err != nil {
		t.Fatalf("ExportVCF on a pool-only sample failed: %v", err)
	}
	// Individual 4 carries neither substrate and stays out; 1, 2 and 3 qualify on the
	// pool alone.
	if vcfStats.NumSamples != 3 {
		t.Fatalf("VCF samples = %d, want 3 (1,2,3 qualify on pool mutations; 4 carries nothing)",
			vcfStats.NumSamples)
	}
	// No bitfield anywhere means no founder sites at all — not a monomorphic row per bit.
	// The emitted mutations are the same three the merged test expects: 100, 101, 102.
	// 105 is fixed in the sample, 103 is out of arm range, 104 has no pool entry.
	if vcfStats.NumFounderSites != 0 {
		t.Errorf("VCF founder sites = %d, want 0: there is no bitfield to segregate",
			vcfStats.NumFounderSites)
	}
	if vcfStats.NumMutationSites != 3 || vcfStats.NumSites != 3 {
		t.Errorf("VCF sites = %d (%d de-novo), want 3 (3)",
			vcfStats.NumSites, vcfStats.NumMutationSites)
	}

	// --- ARGweaver, same sample ---
	model2, pop2 := poolOnlyFixture()
	model2.ResultsDir = t.TempDir()
	argStats, err := ExportARGweaver(model2, pop2, ids, ARGweaverOptions{BpPerBit: 10})
	if err != nil {
		t.Fatalf("ExportARGweaver on a pool-only sample failed: %v", err)
	}
	if argStats.NumHaplotypes != 2*vcfStats.NumSamples {
		t.Errorf("haplotypes = %d, want %d (the VCF's %d samples x 2 strands)",
			argStats.NumHaplotypes, 2*vcfStats.NumSamples, vcfStats.NumSamples)
	}

	// THE POINT: the two files describe the same site set, which is what makes the
	// truth-vs-inference join across them valid.
	if argStats.TotalSites != vcfStats.NumSites {
		t.Errorf(".sites wrote %d sites but the VCF reports %d — the two artifacts must "+
			"agree on the site set", argStats.TotalSites, vcfStats.NumSites)
	}
	if argStats.TotalDropped != 0 {
		t.Errorf("dropped = %d, want 0: three sites cannot overflow a 10 bp bit window",
			argStats.TotalDropped)
	}

	// And the file itself is well-formed and genuinely segregating, with no leftover
	// founder rows.
	names, _, pos, alleles := readSites(t, argStats.Regions[0].Path)
	if len(names) != argStats.NumHaplotypes {
		t.Errorf("NAMES has %d entries, want %d", len(names), argStats.NumHaplotypes)
	}
	if len(pos) != vcfStats.NumSites {
		t.Fatalf("site rows = %d, want %d", len(pos), vcfStats.NumSites)
	}
	for i := range pos {
		if i > 0 && pos[i] <= pos[i-1] {
			t.Fatalf("positions must be sorted and distinct: pos[%d]=%d follows %d",
				i, pos[i], pos[i-1])
		}
		if len(alleles[i]) != argStats.NumHaplotypes {
			t.Fatalf("site at %d has %d allele chars, want %d",
				pos[i], len(alleles[i]), argStats.NumHaplotypes)
		}
		if !strings.Contains(alleles[i], "A") || !strings.Contains(alleles[i], "T") {
			t.Errorf("site at %d is invariant (%q) and should not have been written",
				pos[i], alleles[i])
		}
	}
}

// IncludeFixed must not resurrect the founder walk on a pool-only run. The founder
// branch emits one record per genome bit when fixed sites are kept, so leaving it
// enabled with no bitfield would bury the pool records under a monomorphic row for
// every bit in the genome (10^6 of them on the real createdvar-scale model).
func TestPoolOnlyIncludeFixedStillWritesNoFounderSites(t *testing.T) {
	model, pop := poolOnlyFixture()
	model.ResultsDir = t.TempDir()

	stats, err := ExportVCF(model, pop, []int{1, 2, 3, 4},
		VCFOptions{IncludeMutations: true, IncludeFixed: true})
	if err != nil {
		t.Fatalf("ExportVCF failed: %v", err)
	}
	if stats.NumFounderSites != 0 {
		t.Errorf("founder sites = %d under IncludeFixed, want 0: there is no bitfield",
			stats.NumFounderSites)
	}
	// 105 is fixed in the sample, so keeping fixed sites adds it to the three
	// segregating mutations; 103 is still outside every arm and 104 still unpooled.
	if stats.NumMutationSites != 4 {
		t.Errorf("de-novo sites = %d under IncludeFixed, want 4 (100,101,102,105)",
			stats.NumMutationSites)
	}
}

// An individual carrying NEITHER substrate is still not exportable, and a sample made
// only of such individuals is still an error rather than an empty file.
func TestExportRejectsSampleWithNeitherSubstrate(t *testing.T) {
	model, pop := poolOnlyFixture()
	model.ResultsDir = t.TempDir()

	if _, err := ExportARGweaver(model, pop, []int{4}, ARGweaverOptions{BpPerBit: 10}); err == nil {
		t.Fatal("ExportARGweaver accepted a sample with no chromosomes and no mutations")
	}
}

// The unpooled id must not be what admits an individual: an id carried but absent from
// MutationPool has no Position and can contribute no site, so treating it as evidence of
// exportable content would put an all-ancestral column in both files.
func TestUnpooledMutationDoesNotQualifyAnIndividual(t *testing.T) {
	_, pop := poolOnlyFixture()
	pop.IndMutations[4] = map[int][]int{0: {104}, 1: {}} // 104 is deliberately not in the pool

	samples, withBitfield := exportSamples(pop, []int{1, 2, 3, 4})
	if withBitfield != 0 {
		t.Errorf("withBitfield = %d, want 0 on a pool-only fixture", withBitfield)
	}
	for _, id := range samples {
		if id == 4 {
			t.Error("individual 4 was admitted on an unpooled mutation id, which has no position")
		}
	}
	if len(samples) != 3 {
		t.Errorf("samples = %v, want the three pool carriers 1,2,3", samples)
	}
}
