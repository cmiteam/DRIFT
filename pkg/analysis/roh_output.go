package analysis

// ROH output (roadmap §4). Writes two CSVs mirroring the other analyses'
// end-of-run reporting:
//
//   <model>_roh_per_individual_run%d_year%d.csv — one row per sampled individual:
//     N_ROH, SROH (bits + cM), F_ROH (bits + cM), het count, and the short/medium/
//     long class breakdown. This per-individual table is itself a scientific product
//     (real ROH studies plot N_ROH vs SROH per individual to separate consanguinity
//     from background relatedness).
//
//   <model>_roh_distribution_run%d_year%d.csv — the pooled run-length distribution
//     (the headline product): log2-binned run counts across all sampled individuals,
//     plus the per-class totals. This is the "ROH length distribution" the roadmap
//     names, comparable to the equivalent plot on real 1000G/HGDP panels.
//
// A one-line summary is printed to stdout. No-op (nil error) when no individuals
// carried a bitfield.

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"os"
	"sort"
)

var rohClassNames = [numROHClass]string{"short", "medium", "long"}

// SaveROH writes the per-individual and pooled-distribution CSVs and prints a
// summary. A no-op when the result is nil or holds no sampled individuals.
func SaveROH(model *core.Model, result *ROHResult) error {
	if result == nil || result.NumSampled == 0 {
		fmt.Println("ROH: no individuals with a tracked bitfield (requires track_DNA and a seeded genome); skipped")
		return nil
	}

	if err := saveROHPerIndividual(model, result); err != nil {
		return err
	}
	if err := saveROHDistribution(model, result); err != nil {
		return err
	}

	mapNote := "uniform cM map (recomb_cM_per_bit fallback)"
	if result.HasMap {
		mapNote = "genetic map"
	}
	fmt.Printf("ROH: %d individuals, mean F_ROH=%.4f (bits) / %.4f (cM), mean N_ROH=%.1f "+
		"[short %.1f / medium %.1f / long %.1f], mean het=%.4f — %s\n",
		result.NumSampled, result.MeanFROHBits, result.MeanFROHCM, result.MeanNumRuns,
		result.MeanClassCount[ROHShort], result.MeanClassCount[ROHMedium], result.MeanClassCount[ROHLong],
		result.MeanHetFrac, mapNote)
	return nil
}

func saveROHPerIndividual(model *core.Model, result *ROHResult) error {
	path := rohPath(model, "per_individual")
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create ROH per-individual file: %v", err)
	}
	defer file.Close()

	w := csv.NewWriter(file)
	defer w.Flush()

	if err := w.Write([]string{
		"ID", "HetSites", "N_ROH", "SROH_bits", "SROH_cM", "F_ROH_bits", "F_ROH_cM",
		"N_short", "N_medium", "N_long", "bits_short", "bits_medium", "bits_long",
	}); err != nil {
		return err
	}
	for _, ind := range result.Individuals {
		if err := w.Write([]string{
			fmt.Sprintf("%d", ind.ID),
			fmt.Sprintf("%d", ind.HetSites),
			fmt.Sprintf("%d", ind.NumRuns),
			fmt.Sprintf("%d", ind.TotalBits),
			fmt.Sprintf("%.4f", ind.TotalCM),
			fmt.Sprintf("%.6f", ind.FROHBits),
			fmt.Sprintf("%.6f", ind.FROHCM),
			fmt.Sprintf("%d", ind.ClassCount[ROHShort]),
			fmt.Sprintf("%d", ind.ClassCount[ROHMedium]),
			fmt.Sprintf("%d", ind.ClassCount[ROHLong]),
			fmt.Sprintf("%d", ind.ClassBits[ROHShort]),
			fmt.Sprintf("%d", ind.ClassBits[ROHMedium]),
			fmt.Sprintf("%d", ind.ClassBits[ROHLong]),
		}); err != nil {
			return err
		}
	}
	return nil
}

func saveROHDistribution(model *core.Model, result *ROHResult) error {
	path := rohPath(model, "distribution")
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create ROH distribution file: %v", err)
	}
	defer file.Close()

	w := csv.NewWriter(file)
	defer w.Flush()

	// Header block: config + totals as comment-style leading rows would break the
	// CSV, so the length distribution is the table and the config rides the summary
	// print. Columns: the log2 length bin, its bit range, and the pooled count.
	if err := w.Write([]string{"bin_log2", "bits_low", "bits_high", "count"}); err != nil {
		return err
	}
	bins := make([]int, 0, len(result.LengthHist))
	for b := range result.LengthHist {
		bins = append(bins, b)
	}
	sort.Ints(bins)
	for _, b := range bins {
		low := 1 << uint(b)
		high := 1 << uint(b+1) // exclusive
		if err := w.Write([]string{
			fmt.Sprintf("%d", b),
			fmt.Sprintf("%d", low),
			fmt.Sprintf("%d", high),
			fmt.Sprintf("%d", result.LengthHist[b]),
		}); err != nil {
			return err
		}
	}
	// Per-class pooled totals (bin_log2 = -1 sentinel rows).
	for k := 0; k < numROHClass; k++ {
		if err := w.Write([]string{
			"-1", "class", rohClassNames[k], fmt.Sprintf("%d", result.ClassTotals[k]),
		}); err != nil {
			return err
		}
	}
	return nil
}

// (bin b covers run lengths [2^b, 2^(b+1)); logBin from ibd.go supplies b.)

// rohPath builds a results-directory path for one of the ROH outputs, matching the
// other analyses' <model>_<kind>_run%d_year%d.csv naming.
func rohPath(model *core.Model, kind string) string {
	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	return fmt.Sprintf("%s/%s_roh_%s_run%d_year%d.csv", model.ResultsDir, model.ModelName, kind, run, year)
}
