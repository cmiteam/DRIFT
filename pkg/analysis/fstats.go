package analysis

// f-statistics: f2, f3, f4, and the D-statistic (ABBA-BABA) — roadmap §6b,
// "Differentiation & admixture statistics". These are the Patterson/Reich tools
// used to argue admixture, tree topology, and "Africa as an outgroup"; they are
// how the out-of-Africa claims are actually *made*, so reproducing them is
// non-negotiable for testing/countering OoA.
//
// POPULATIONS. As with Fst (fst.go), "populations" here are the discrete DEMES of
// the island model (num_demes / deme_migration_rate). Every sampled individual
// carries a Deme label; f-stats need >= 2 (f2), >= 3 (f3) or >= 4 (f4/D) eligible
// demes. This file reuses the fst.go substrate: demeMembers (group sampled ids by
// deme) and countSite (per-site derived/het counts), plus buildContigs/bitAt.
//
// WHAT A "SITE" AND "FREQUENCY" MEAN. A genome bit is one biallelic SNP: bit 0 =
// ancestral, bit 1 = derived. For deme X at a site, p_X is the SAMPLE frequency of
// the derived allele = (derived-allele count) / n_X, where n_X = 2 * (sampled
// diploids in X) is the number of sampled haplotypes. All statistics are computed
// as a ratio/mean of sums across the polymorphic sites (never a mean of per-site
// ratios), the standard low-bias approach.
//
// ESTIMATORS (sample-size-unbiased; Patterson et al. 2012, Peter 2016). Per site:
//
//   f2(A,B)   = (p_A - p_B)^2  - p_A(1-p_A)/(n_A-1) - p_B(1-p_B)/(n_B-1)
//   f3(T;B,C) = (p_T - p_B)(p_T - p_C) - p_T(1-p_T)/(n_T-1)
//   f4(A,B;C,D) = (p_A - p_B)(p_C - p_D)                       [no bias term needed]
//
//   The bias term subtracts each population's sampling variance; it is applied only
//   to frequencies that appear SQUARED (both appearances in f2; the target T in f3).
//   f4 needs none because the four demes are distinct so every cross term is a
//   product of independent unbiased estimates.
//
//   f2/f3/f4 are reported as the MEAN over the polymorphic sites used (Σ / S).
//
//   D-statistic (ABBA-BABA), slots (A,B,C,D) = (P1,P2,P3,P4), P4 the outgroup, in
//   the symmetric form that allows outgroup polymorphism:
//       ABBA = (1-p1)p2 p3(1-p4) + p1(1-p2)(1-p3)p4
//       BABA = p1(1-p2)p3(1-p4) + (1-p1)p2(1-p3)p4
//       D = Σ(ABBA - BABA) / Σ(ABBA + BABA)
//   Algebraically ABBA-BABA = (p2-p1)(p3-p4), so D and f4(A,B;C,D) share the same
//   per-site kernel with OPPOSITE sign (f4 kernel = (p1-p2)(p3-p4) = -(ABBA-BABA)).
//   D > 0 means an excess of ABBA — P2 and P3 share more derived alleles than
//   expected under the tree (((P1,P2),P3),P4): gene flow between P2 and P3.
//
// INTERPRETATION. f3(T;B,C) < 0 is a hallmark of ADMIXTURE in T (T is intermediate
// between B and C, which no unadmixed population can be). f4/D != 0 rejects the
// unrooted tree topology (A,B),(C,D) and is the admixture / tree-topology test.
//
// WHICH DEMES FILL THE SLOTS (hybrid, chosen 2026-07-06). f2 is symmetric and is
// always computed for every unordered deme pair. f3/f4/D are driven by the
// fstats_demes string param:
//   - empty (default): AUTO — every deme as an f3 target against every other pair,
//     and every 4-deme quartet in its 3 independent f4/D pairings.
//   - "T;B,C"  : a single f3(T; B, C).
//   - "A,B;C,D": a single f4(A,B;C,D) and its D-statistic.
// Wired into drift.go behind track_fstats (default off); requires track_DNA and
// num_demes >= 2. Sample size via fstats_sample_size (0 = all).

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
)

// F2Stat is f2 for one unordered deme pair (mean over sites used).
type F2Stat struct {
	A, B     int
	Value    float64
	NumSites int
}

// F3Stat is f3(Target; B, C) (mean over sites used). Value < 0 signals admixture
// in the target deme.
type F3Stat struct {
	Target, B, C int
	Value        float64
	NumSites     int
}

// F4Stat carries both f4(A,B;C,D) and the corresponding D-statistic for one ordered
// quadruple; they share the per-site kernel (see file header).
type F4Stat struct {
	A, B, C, D int     // deme slots (D is the outgroup slot P4)
	F4         float64 // f4(A,B;C,D)
	Dstat      float64 // D-statistic (ABBA-BABA)
	NumSites   int
}

// FStatsResult is the output of an f-statistics pass.
type FStatsResult struct {
	Demes     []int       // eligible demes (>= 2 sampled diploids), sorted
	DemeSizes map[int]int // sampled diploids per eligible deme
	F2        []F2Stat
	F3        []F3Stat
	F4        []F4Stat
	NumSites  int    // polymorphic sites used (shared denominator S)
	Mode      string // "auto" or "explicit" (for the log line)
}

// ComputeFStats runs the f-statistics pass over the given sampled individual IDs
// (see SampleIDs), grouping them into demes by Deme label and reusing the fst.go
// site-counting substrate. Returns an empty result when there are fewer than two
// eligible demes (or no polymorphic sites). The f3/f4/D combinations computed are
// driven by the fstats_demes param (see file header).
func ComputeFStats(model *core.Model, pop *core.Pop, ids []int) *FStatsResult {
	byDeme := demeMembers(pop, ids)

	// Eligible demes have >= 2 sampled diploids (the bias corrections divide by
	// n-1 in haplotypes; 2 diploids => n=4). Sorted for deterministic output.
	demes := make([]int, 0, len(byDeme))
	for d, members := range byDeme {
		if len(members) >= 2 {
			demes = append(demes, d)
		} else if len(members) > 0 {
			fmt.Printf("f-stats: deme %d has only %d sampled diploid(s); excluded (need >= 2)\n", d, len(members))
		}
	}
	sort.Ints(demes)

	res := &FStatsResult{Demes: demes, DemeSizes: map[int]int{}}
	for _, d := range demes {
		res.DemeSizes[d] = len(byDeme[d])
	}
	if len(demes) < 2 {
		fmt.Printf("f-stats: %d eligible deme(s); need >= 2 (set num_demes and run long enough for demes to persist)\n", len(demes))
		return res
	}

	// Deme -> index into the per-site frequency slice.
	idx := make(map[int]int, len(demes))
	demeN := make([]int, len(demes)) // sampled haplotypes per deme (constant)
	for i, d := range demes {
		idx[d] = i
		demeN[i] = 2 * len(byDeme[d])
	}

	// Plan which f3 / f4 combinations to compute (f2 is always all pairs).
	spec := strings.TrimSpace(model.StringParam("fstats_demes", ""))
	f3plan, f4plan, mode, warn := planFStats(spec, demes)
	if warn != "" {
		fmt.Printf("f-stats: %s\n", warn)
	}
	res.Mode = mode

	// Accumulators (ratio/mean of sums across sites).
	type f2acc struct{ sum float64 }
	f2acc_ := make([]f2acc, 0)
	f2pairs := make([][2]int, 0) // (demeA, demeB)
	for i := 0; i < len(demes); i++ {
		for j := i + 1; j < len(demes); j++ {
			f2pairs = append(f2pairs, [2]int{demes[i], demes[j]})
			f2acc_ = append(f2acc_, f2acc{})
		}
	}
	type f3acc struct{ sum float64 }
	f3acc_ := make([]f3acc, len(f3plan))
	type f4acc struct {
		f4sum        float64
		abbaS, babaS float64
	}
	f4acc_ := make([]f4acc, len(f4plan))

	// Walk every genome position via the chromosome arms (same scan as VCF / Fst),
	// tallying per-deme derived counts once per site.
	contigs := buildContigs(model)
	p := make([]float64, len(demes)) // reused per site: sample derived freq per deme
	S := 0
	for _, ct := range contigs {
		for _, arm := range ct.arms {
			start, length := arm[0], arm[1]
			for bit := start; bit < start+length; bit++ {
				totalDerived, totalAlleles := 0, 0
				for i, d := range demes {
					sc := countSite(pop, byDeme[d], bit)
					p[i] = float64(sc.derived) / float64(demeN[i])
					totalDerived += sc.derived
					totalAlleles += demeN[i]
				}
				// Skip sites monomorphic across the eligible demes.
				if totalDerived == 0 || totalDerived == totalAlleles {
					continue
				}
				S++

				for k, pair := range f2pairs {
					ia, ib := idx[pair[0]], idx[pair[1]]
					f2acc_[k].sum += f2Site(p[ia], demeN[ia], p[ib], demeN[ib])
				}
				for k, t := range f3plan {
					it, ib, ic := idx[t[0]], idx[t[1]], idx[t[2]]
					f3acc_[k].sum += f3Site(p[it], demeN[it], p[ib], p[ic])
				}
				for k, q := range f4plan {
					ia, ib, ic, id := idx[q[0]], idx[q[1]], idx[q[2]], idx[q[3]]
					f4acc_[k].f4sum += f4Site(p[ia], p[ib], p[ic], p[id])
					abba, baba := abbaBaba(p[ia], p[ib], p[ic], p[id])
					f4acc_[k].abbaS += abba
					f4acc_[k].babaS += baba
				}
			}
		}
	}

	res.NumSites = S
	if S == 0 {
		fmt.Println("f-stats: no polymorphic sites across the eligible demes")
		return res
	}

	fs := float64(S)
	for k, pair := range f2pairs {
		res.F2 = append(res.F2, F2Stat{A: pair[0], B: pair[1], Value: f2acc_[k].sum / fs, NumSites: S})
	}
	for k, t := range f3plan {
		res.F3 = append(res.F3, F3Stat{Target: t[0], B: t[1], C: t[2], Value: f3acc_[k].sum / fs, NumSites: S})
	}
	for k, q := range f4plan {
		res.F4 = append(res.F4, F4Stat{
			A: q[0], B: q[1], C: q[2], D: q[3],
			F4:       f4acc_[k].f4sum / fs,
			Dstat:    ratio(f4acc_[k].abbaS-f4acc_[k].babaS, f4acc_[k].abbaS+f4acc_[k].babaS),
			NumSites: S,
		})
	}
	return res
}

// --- Per-site kernels ---------------------------------------------------------

// f2Site is the unbiased per-site f2 for two demes with sample derived frequencies
// pA,pB and sampled haplotype counts nA,nB.
func f2Site(pA float64, nA int, pB float64, nB int) float64 {
	return (pA-pB)*(pA-pB) - pA*(1-pA)/float64(nA-1) - pB*(1-pB)/float64(nB-1)
}

// f3Site is the unbiased per-site f3(T;B,C); only the target T carries a bias term
// (it is the frequency that appears in both factors).
func f3Site(pT float64, nT int, pB, pC float64) float64 {
	return (pT-pB)*(pT-pC) - pT*(1-pT)/float64(nT-1)
}

// f4Site is the per-site f4(A,B;C,D); unbiased with no correction for four distinct
// demes.
func f4Site(pA, pB, pC, pD float64) float64 {
	return (pA - pB) * (pC - pD)
}

// abbaBaba returns the per-site ABBA and BABA weights (symmetric form allowing an
// polymorphic outgroup) for slots (p1,p2,p3,p4) = (A,B,C,D).
func abbaBaba(p1, p2, p3, p4 float64) (abba, baba float64) {
	abba = (1-p1)*p2*p3*(1-p4) + p1*(1-p2)*(1-p3)*p4
	baba = p1*(1-p2)*p3*(1-p4) + (1-p1)*p2*(1-p3)*p4
	return abba, baba
}

// --- Planning -----------------------------------------------------------------

// planFStats decides which f3 and f4/D combinations to compute. With an empty spec
// it auto-enumerates all of them; otherwise it parses a single explicit spec of the
// form "T;B,C" (f3) or "A,B;C,D" (f4/D). Returns f3 triples and f4 quadruples of
// DEME IDs, a mode string, and a warning message (empty if none).
func planFStats(spec string, demes []int) (f3 [][3]int, f4 [][4]int, mode, warn string) {
	eligible := make(map[int]bool, len(demes))
	for _, d := range demes {
		eligible[d] = true
	}

	if spec == "" {
		mode = "auto"
		// Every deme as an f3 target against every unordered pair of the others.
		for _, t := range demes {
			others := make([]int, 0, len(demes)-1)
			for _, d := range demes {
				if d != t {
					others = append(others, d)
				}
			}
			for i := 0; i < len(others); i++ {
				for j := i + 1; j < len(others); j++ {
					f3 = append(f3, [3]int{t, others[i], others[j]})
				}
			}
		}
		// Every 4-deme quartet in its 3 independent f4/D pairings.
		for a := 0; a < len(demes); a++ {
			for b := a + 1; b < len(demes); b++ {
				for c := b + 1; c < len(demes); c++ {
					for d := c + 1; d < len(demes); d++ {
						w, x, y, z := demes[a], demes[b], demes[c], demes[d]
						f4 = append(f4,
							[4]int{w, x, y, z},
							[4]int{w, y, x, z},
							[4]int{w, z, x, y})
					}
				}
			}
		}
		return f3, f4, mode, ""
	}

	mode = "explicit"
	parts := strings.Split(spec, ";")
	if len(parts) != 2 {
		return nil, nil, mode, fmt.Sprintf("fstats_demes %q malformed; expected \"T;B,C\" (f3) or \"A,B;C,D\" (f4/D) — skipping f3/f4", spec)
	}
	left, err1 := parseDemeList(parts[0])
	right, err2 := parseDemeList(parts[1])
	if err1 != nil || err2 != nil {
		return nil, nil, mode, fmt.Sprintf("fstats_demes %q has a non-integer deme — skipping f3/f4", spec)
	}
	all := append(append([]int{}, left...), right...)
	for _, d := range all {
		if !eligible[d] {
			return nil, nil, mode, fmt.Sprintf("fstats_demes %q names deme %d which is not eligible (need >= 2 sampled diploids) — skipping f3/f4", spec, d)
		}
	}
	switch {
	case len(left) == 1 && len(right) == 2:
		f3 = append(f3, [3]int{left[0], right[0], right[1]})
	case len(left) == 2 && len(right) == 2:
		f4 = append(f4, [4]int{left[0], left[1], right[0], right[1]})
	default:
		return nil, nil, mode, fmt.Sprintf("fstats_demes %q: expected 1 deme before ';' and 2 after (f3), or 2 and 2 (f4/D) — skipping f3/f4", spec)
	}
	return f3, f4, mode, ""
}

// parseDemeList parses a comma-separated list of deme integers, ignoring blanks.
func parseDemeList(s string) ([]int, error) {
	out := make([]int, 0, 4)
	for _, tok := range strings.Split(s, ",") {
		tok = strings.TrimSpace(tok)
		if tok == "" {
			continue
		}
		v, err := strconv.Atoi(tok)
		if err != nil {
			return nil, err
		}
		out = append(out, v)
	}
	return out, nil
}

// --- Output -------------------------------------------------------------------

// SaveFStats writes a single CSV holding every computed statistic (f2/f3/f4/D) in a
// uniform schema, plus a one-line summary. A no-op (nil) when the result has no
// statistics (fewer than two eligible demes or no polymorphic sites).
func SaveFStats(model *core.Model, result *FStatsResult) error {
	if result == nil || (len(result.F2) == 0 && len(result.F3) == 0 && len(result.F4) == 0) {
		return nil
	}

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	path := fmt.Sprintf("%s/%s_fstats_run%d_year%d.csv",
		model.ResultsDir, model.ModelName, run, year)

	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create f-stats file: %v", err)
	}
	defer file.Close()

	w := csv.NewWriter(file)
	defer w.Flush()

	// Uniform schema: unused deme slots are blank. Stat is one of f2/f3/f4/D.
	if err := w.Write([]string{"Stat", "DemeA", "DemeB", "DemeC", "DemeD", "Value", "NumSites"}); err != nil {
		return err
	}
	itoa := func(i int) string { return strconv.Itoa(i) }
	f6 := func(v float64) string { return fmt.Sprintf("%.6f", v) }

	for _, s := range result.F2 {
		if err := w.Write([]string{"f2", itoa(s.A), itoa(s.B), "", "", f6(s.Value), itoa(s.NumSites)}); err != nil {
			return err
		}
	}
	for _, s := range result.F3 {
		// DemeA is the target T in f3(T; B, C).
		if err := w.Write([]string{"f3", itoa(s.Target), itoa(s.B), itoa(s.C), "", f6(s.Value), itoa(s.NumSites)}); err != nil {
			return err
		}
	}
	for _, s := range result.F4 {
		if err := w.Write([]string{"f4", itoa(s.A), itoa(s.B), itoa(s.C), itoa(s.D), f6(s.F4), itoa(s.NumSites)}); err != nil {
			return err
		}
		if err := w.Write([]string{"D", itoa(s.A), itoa(s.B), itoa(s.C), itoa(s.D), f6(s.Dstat), itoa(s.NumSites)}); err != nil {
			return err
		}
	}

	fmt.Printf("f-stats (%s): %d demes, %d polymorphic sites — %d f2, %d f3, %d f4/D → %s\n",
		result.Mode, len(result.Demes), result.NumSites,
		len(result.F2), len(result.F3), len(result.F4), path)
	return nil
}
