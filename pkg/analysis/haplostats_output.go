package analysis

// §6g output: CSVs in selscan/rehh column conventions, plus optional selscan
// .hap/.map export so the same numbers can be reproduced with external tools.
//
//   <model>_ihs_*      per-core iHS scan  (selscan --ihs shape)
//   <model>_ehh_*      EHH-vs-cM decay curves for the top-|iHS| cores (rehh calc_ehh shape)
//   <model>_xpehh_*    cross-population XP-EHH per site (when two demes present)
//   <model>_lddecay_*  mean r^2 by cM-distance bin (the LD-decay curve)
//   <model>_run*.hap / .map   sampled pool haplotypes in selscan format (opt-in)

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"sort"
)

func isNaN(f float64) bool   { return math.IsNaN(f) }
func isInf(f float64) bool   { return math.IsInf(f, 0) }
func absF(f float64) float64 { return math.Abs(f) }

// SaveHaploStats writes the iHS / EHH-curve / XP-EHH / LD-decay CSVs and prints a
// one-line summary (matching the other analyses' end-of-run reporting).
func SaveHaploStats(model *core.Model, res *HaploStatsResult) error {
	if res == nil {
		return nil
	}
	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	base := fmt.Sprintf("%s/%s", model.ResultsDir, model.ModelName)

	if err := writeIHS(fmt.Sprintf("%s_ihs_run%d_year%d.csv", base, run, year), res); err != nil {
		return err
	}
	if err := writeEHHCurves(fmt.Sprintf("%s_ehh_run%d_year%d.csv", base, run, year), res); err != nil {
		return err
	}
	if len(res.XPEHH) > 0 {
		if err := writeXPEHH(fmt.Sprintf("%s_xpehh_run%d_year%d.csv", base, run, year), res); err != nil {
			return err
		}
	}
	if err := writeLDDecay(fmt.Sprintf("%s_lddecay_run%d_year%d.csv", base, run, year), res); err != nil {
		return err
	}

	// Headline: number of |StdIHS|>2 cores (the conventional selection flag).
	nOut := 0
	for _, c := range res.Cores {
		if !isNaN(c.StdIHS) && absF(c.StdIHS) > 2 {
			nOut++
		}
	}
	xp := ""
	if res.DemeA >= 0 {
		xp = fmt.Sprintf(", XP-EHH demes %d/%d over %d sites", res.DemeA, res.DemeB, len(res.XPEHH))
	}
	fmt.Printf("HaploStats: %d cores (iHS), %d with |iHS|>2%s\n", res.NumCores, nOut, xp)
	return nil
}

func writeIHS(path string, res *HaploStatsResult) error {
	cores := make([]CoreEHH, len(res.Cores))
	copy(cores, res.Cores)
	sort.Slice(cores, func(i, j int) bool {
		if cores[i].Chrom != cores[j].Chrom {
			return cores[i].Chrom < cores[j].Chrom
		}
		if cores[i].Position != cores[j].Position {
			return cores[i].Position < cores[j].Position
		}
		return cores[i].MutID < cores[j].MutID
	})
	return withCSV(path, func(w *csv.Writer) error {
		if err := w.Write([]string{"MutID", "Chrom", "Position", "DerivedCount", "Nhap", "DerivedFreq", "IHH_Derived", "IHH_Ancestral", "UnstdIHS", "StdIHS"}); err != nil {
			return err
		}
		for _, c := range cores {
			if err := w.Write([]string{
				itoa(c.MutID), itoa(c.Chrom), itoa(c.Position), itoa(c.DerivedCount), itoa(c.Nhap),
				ftoa(c.DerivedFreq), ftoa(c.IHHDerived), ftoa(c.IHHAncestral), ftoa(c.UnstdIHS), ftoa(c.StdIHS),
			}); err != nil {
				return err
			}
		}
		return nil
	})
}

func writeEHHCurves(path string, res *HaploStatsResult) error {
	return withCSV(path, func(w *csv.Writer) error {
		if err := w.Write([]string{"MutID", "Chrom", "Position", "CM_from_core", "EHH_Derived", "EHH_Ancestral"}); err != nil {
			return err
		}
		for _, p := range res.Curves {
			if err := w.Write([]string{
				itoa(p.MutID), itoa(p.Chrom), itoa(p.Position), ftoa(p.CM), ftoa(p.EHHDerived), ftoa(p.EHHAncestral),
			}); err != nil {
				return err
			}
		}
		return nil
	})
}

func writeXPEHH(path string, res *HaploStatsResult) error {
	sites := make([]XPEHHSite, len(res.XPEHH))
	copy(sites, res.XPEHH)
	sort.Slice(sites, func(i, j int) bool {
		if sites[i].Chrom != sites[j].Chrom {
			return sites[i].Chrom < sites[j].Chrom
		}
		return sites[i].Position < sites[j].Position
	})
	return withCSV(path, func(w *csv.Writer) error {
		if err := w.Write([]string{"MutID", "Chrom", "Position", "DemeA", "DemeB", "iES_A", "iES_B", "UnstdXPEHH", "StdXPEHH"}); err != nil {
			return err
		}
		for _, s := range sites {
			if err := w.Write([]string{
				itoa(s.MutID), itoa(s.Chrom), itoa(s.Position), itoa(res.DemeA), itoa(res.DemeB),
				ftoa(s.IESA), ftoa(s.IESB), ftoa(s.UnstdXPEHH), ftoa(s.StdXPEHH),
			}); err != nil {
				return err
			}
		}
		return nil
	})
}

func writeLDDecay(path string, res *HaploStatsResult) error {
	return withCSV(path, func(w *csv.Writer) error {
		if err := w.Write([]string{"LoCM", "HiCM", "MeanR2", "NumPairs"}); err != nil {
			return err
		}
		for _, b := range res.LDDecay {
			if err := w.Write([]string{ftoa(b.LoCM), ftoa(b.HiCM), ftoa(b.MeanR2), itoa(b.NumPairs)}); err != nil {
				return err
			}
		}
		return nil
	})
}

// ExportHapMap writes the sampled pool haplotypes as selscan .hap/.map files
// (opt-in via haplostats_export_hap). The .map has one row per site
// (chrom, locusID, geneticPos_cM, physicalPos); the .hap has one row per
// haplotype (2*nsamples) of space-separated 0/1 alleles in .map column order.
// Deterministic (no RNG). Rebuilds the workspace so it is independent of the
// stats pass.
func ExportHapMap(model *core.Model, pop *core.Pop, ids []int) error {
	arms := armRangesOf(model)
	if len(arms) == 0 {
		return nil
	}
	gm := haploGeneticMap(model)
	cw := buildChromHaps(pop, ids, arms, gm)

	// Global, deterministic site order: chromosome then position then mutID.
	type gsite struct {
		chrom, pos, mutID int
		cm                float64
	}
	var sites []gsite
	for _, ch := range sortedChromKeys(cw) {
		for _, s := range cw[ch].sites {
			sites = append(sites, gsite{ch, s.pos, s.mutID, s.cumCM})
		}
	}
	sort.Slice(sites, func(i, j int) bool {
		if sites[i].chrom != sites[j].chrom {
			return sites[i].chrom < sites[j].chrom
		}
		if sites[i].pos != sites[j].pos {
			return sites[i].pos < sites[j].pos
		}
		return sites[i].mutID < sites[j].mutID
	})

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	base := fmt.Sprintf("%s/%s_run%d_year%d", model.ResultsDir, model.ModelName, run, year)

	// .map
	mapFile, err := os.Create(base + ".map")
	if err != nil {
		return fmt.Errorf("failed to create hap-map file: %v", err)
	}
	for _, s := range sites {
		fmt.Fprintf(mapFile, "%d\tmut%d\t%.6f\t%d\n", s.chrom, s.mutID, s.cm, s.pos)
	}
	mapFile.Close()

	// .hap — one row per haplotype (strand), aligned across chromosomes. A site's
	// mutID is chromosome-local in the workspace, so look it up in the owning
	// chromosome's membership.
	hapFile, err := os.Create(base + ".hap")
	if err != nil {
		return fmt.Errorf("failed to create hap file: %v", err)
	}
	defer hapFile.Close()
	nhap := 2 * len(ids)
	for h := 0; h < nhap; h++ {
		for i, s := range sites {
			if i > 0 {
				fmt.Fprint(hapFile, " ")
			}
			bit := 0
			if cw[s.chrom].member[h][s.mutID] {
				bit = 1
			}
			fmt.Fprint(hapFile, bit)
		}
		fmt.Fprintln(hapFile)
	}
	fmt.Printf("HaploStats: exported %d haplotypes x %d sites -> %s.hap/.map\n", nhap, len(sites), base)
	return nil
}

// ---- small helpers ----

func withCSV(path string, f func(*csv.Writer) error) error {
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create %s: %v", path, err)
	}
	defer file.Close()
	w := csv.NewWriter(file)
	defer w.Flush()
	return f(w)
}

func itoa(i int) string { return fmt.Sprintf("%d", i) }

func ftoa(f float64) string {
	if isNaN(f) {
		return "NA"
	}
	if isInf(f) {
		if f > 0 {
			return "Inf"
		}
		return "-Inf"
	}
	return fmt.Sprintf("%.6g", f)
}
