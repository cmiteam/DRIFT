package analysis

// Ne-through-time output (roadmap §6d). Writes the per-window Ne series that
// CaptureNe accumulated in pop.NeHistory, plus a one-row summary that contrasts the
// trajectory with a single long-term number — the natural companion to the §6h
// coalescent Ne (theta_W/(2*mu), emergent Ne/N ≈ 0.17). Mirrors the SFS-time-series
// emitter (SaveSFSTimeSeries): rows are captured windows, in year order.

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"sort"
)

// neCell formats an Ne value for CSV: a finite value at two decimals, or empty for a
// non-finite (unresolved) estimate.
func neCell(v float64) string {
	if math.IsInf(v, 0) || math.IsNaN(v) {
		return ""
	}
	return fmt.Sprintf("%.2f", v)
}

// SaveNeTimeSeries writes the Ne-through-time series and a summary row for the run.
// A no-op (nil) when no windows were captured (fewer than two, so no temporal
// interval exists). The series file is <model>_ne_timeseries_run%d.csv with one row
// per captured window; the summary file is <model>_ne_summary_run%d.csv.
func SaveNeTimeSeries(model *core.Model, pop *core.Pop) error {
	if len(pop.NeHistory) < 2 {
		fmt.Printf("Ne-through-time: %d window(s) captured; need >= 2 for a temporal interval (raise the run length or lower ne_interval)\n", len(pop.NeHistory))
		return nil
	}

	// Windows in year order; the first window has no earlier sample (its estimates
	// are zero-valued) and is skipped.
	years := make([]int, 0, len(pop.NeHistory))
	for y := range pop.NeHistory {
		years = append(years, y)
	}
	sort.Ints(years)

	// Union of deme keys across windows, for stable per-deme columns.
	demeSet := map[int]bool{}
	for _, y := range years {
		for d := range pop.NeHistory[y].DemeNe {
			demeSet[d] = true
		}
	}
	demes := make([]int, 0, len(demeSet))
	for d := range demeSet {
		demes = append(demes, d)
	}
	sort.Ints(demes)

	run := model.FreeParameters["run"]
	seriesPath := fmt.Sprintf("%s/%s_ne_timeseries_run%d.csv", model.ResultsDir, model.ModelName, run)
	file, err := os.Create(seriesPath)
	if err != nil {
		return fmt.Errorf("failed to create Ne time-series file: %v", err)
	}
	defer file.Close()
	w := csv.NewWriter(file)
	defer w.Flush()

	header := []string{"Year", "CensusN", "IntervalGen", "TemporalNe", "Ne_over_N", "TemporalNumSites", "LD_Ne", "LD_NumPairs"}
	for _, d := range demes {
		header = append(header, fmt.Sprintf("Deme%d_Ne", d))
	}
	if err := w.Write(header); err != nil {
		return err
	}

	// Accumulators for the summary. The right way to average Ne over time is the
	// harmonic mean (drift compounds like 1/Ne), over resolved windows only.
	var invSum float64
	var resolved int
	var censusSum float64
	var censusCount int

	for _, y := range years {
		s := pop.NeHistory[y]
		if s.IntervalGen <= 0 {
			continue // the first (reference) window has no interval
		}
		neOverN := ""
		if s.CensusN > 0 && !math.IsInf(s.TemporalNe, 0) && !math.IsNaN(s.TemporalNe) {
			neOverN = fmt.Sprintf("%.4f", s.TemporalNe/float64(s.CensusN))
		}
		row := []string{
			fmt.Sprintf("%d", s.Year),
			fmt.Sprintf("%d", s.CensusN),
			fmt.Sprintf("%.3f", s.IntervalGen),
			neCell(s.TemporalNe),
			neOverN,
			fmt.Sprintf("%d", s.TemporalNumSites),
			neCell(s.LDNe),
			fmt.Sprintf("%d", s.LDNumPairs),
		}
		for _, d := range demes {
			if v, ok := s.DemeNe[d]; ok {
				row = append(row, neCell(v))
			} else {
				row = append(row, "")
			}
		}
		if err := w.Write(row); err != nil {
			return err
		}

		if !math.IsInf(s.TemporalNe, 0) && !math.IsNaN(s.TemporalNe) && s.TemporalNe > 0 {
			invSum += 1.0 / s.TemporalNe
			resolved++
		}
		censusSum += float64(s.CensusN)
		censusCount++
	}

	// Summary row: harmonic-mean temporal Ne over resolved windows vs mean census N.
	harmonicNe := math.Inf(1)
	if resolved > 0 && invSum > 0 {
		harmonicNe = float64(resolved) / invSum
	}
	meanCensus := 0.0
	if censusCount > 0 {
		meanCensus = censusSum / float64(censusCount)
	}
	summaryPath := fmt.Sprintf("%s/%s_ne_summary_run%d.csv", model.ResultsDir, model.ModelName, run)
	if err := writeNeSummary(summaryPath, harmonicNe, meanCensus, resolved, censusCount); err != nil {
		return err
	}

	fmt.Printf("Ne-through-time: %d resolved window(s), harmonic-mean temporal Ne = %s, mean census N = %.0f (Ne/N = %s) -> %s\n",
		resolved, neCell(harmonicNe), meanCensus, neOverNStr(harmonicNe, meanCensus), seriesPath)
	return nil
}

// neOverNStr formats the harmonic-mean Ne / mean census ratio, or "n/a".
func neOverNStr(ne, censusN float64) string {
	if censusN <= 0 || math.IsInf(ne, 0) || math.IsNaN(ne) {
		return "n/a"
	}
	return fmt.Sprintf("%.3f", ne/censusN)
}

// writeNeSummary writes the one-row Ne summary contrasting the harmonic-mean temporal
// Ne with the mean census size — the counterpart to the §6h long-term coalescent Ne.
func writeNeSummary(path string, harmonicNe, meanCensus float64, resolved, windows int) error {
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create Ne summary file: %v", err)
	}
	defer file.Close()
	w := csv.NewWriter(file)
	defer w.Flush()

	if err := w.Write([]string{"HarmonicMeanTemporalNe", "MeanCensusN", "Ne_over_N", "ResolvedWindows", "TotalWindows"}); err != nil {
		return err
	}
	return w.Write([]string{
		neCell(harmonicNe),
		fmt.Sprintf("%.1f", meanCensus),
		neOverNStr(harmonicNe, meanCensus),
		fmt.Sprintf("%d", resolved),
		fmt.Sprintf("%d", windows),
	})
}
