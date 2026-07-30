package analysis

// Geographic Fst — spatial population differentiation (roadmap §3, the readout that
// makes barriers / habitat / spatial mating "show structure end-to-end").
//
// The §6b Fst (fst.go) partitions the sampled individuals by their aspatial Deme label
// (individual.Deme) — the discrete island-model subpopulation each individual belongs to.
// But DRIFT's spatial machinery (movement barriers, habitat suitability, distance mating)
// operates purely on the map coordinates (individual.Lat / individual.Lon) and never
// touches the Deme label: a barrier that splits a single-deme population in two produces
// real allele-frequency divergence ACROSS THE MAP while every individual keeps Deme = 0.
// Deme-Fst is therefore blind to it (one deme ⇒ empty result). Geographic Fst closes that
// gap by partitioning the SAME estimator core (computeFstFromPartition) on GEOGRAPHY: the
// map is diced into a grid of rectangular regions and each region is treated as a
// population. This is what lets a barrier, a habitat gradient, or isolation-by-distance
// register as a measurable Fst without any deme labels at all.
//
// Configured by two params (Analysis group): `track_geographic_fst` enables the pass, and
// `geographic_fst_grid` sets the region grid as "<rows>x<cols>" (e.g. "2x2", "4x3"); a
// bare "<n>" means n×n, and an unset/malformed grid defaults to 2x2. Requires track_map (a
// loaded map, so Lat/Lon are meaningful and the grid has real dimensions) and track_DNA
// (chromosomes to compute frequencies over); it reuses the Fst sampling knob
// fst_sample_size. The pass is read-only and end-of-run, drawing no RNG, so leaving
// track_geographic_fst off is a strict no-op (byte-identical).
//
// The region axis is deliberately a coarse rectangular grid — the simplest geography that
// makes the point. A finer partition (per-cell, or by connected habitable component, or a
// continuous isolation-by-distance / Fst-by-distance curve) can slot in behind
// geographicMembers later without touching the estimator core.

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// GeoGrid is a parsed geographic_fst_grid: the region grid dimensions.
type GeoGrid struct {
	Rows int
	Cols int
}

// parseGeoGrid parses a "<rows>x<cols>" grid spec (case-insensitive 'x'); a bare "<n>"
// means n×n. A missing/malformed/non-positive spec falls back to 2x2 (the smallest grid
// that yields differentiation) so an enabled pass always has a usable grid.
func parseGeoGrid(spec string) GeoGrid {
	def := GeoGrid{Rows: 2, Cols: 2}
	s := strings.TrimSpace(strings.ToLower(spec))
	if s == "" {
		return def
	}
	parts := strings.SplitN(s, "x", 2)
	rows, err1 := strconv.Atoi(strings.TrimSpace(parts[0]))
	if err1 != nil || rows < 1 {
		return def
	}
	cols := rows // bare "<n>" => n×n
	if len(parts) == 2 {
		c, err2 := strconv.Atoi(strings.TrimSpace(parts[1]))
		if err2 != nil || c < 1 {
			return def
		}
		cols = c
	}
	return GeoGrid{Rows: rows, Cols: cols}
}

// regionOf maps a (lat, lon) cell to a region id in [0, Rows*Cols) for the given map
// dimensions. Coordinates are clamped into range so an out-of-bounds individual still
// lands in an edge region rather than being dropped. Region id = latBin*Cols + lonBin.
func (g GeoGrid) regionOf(lat, lon, mapH, mapW int) int {
	latBin := lat * g.Rows / max(mapH, 1)
	lonBin := lon * g.Cols / max(mapW, 1)
	if latBin < 0 {
		latBin = 0
	} else if latBin >= g.Rows {
		latBin = g.Rows - 1
	}
	if lonBin < 0 {
		lonBin = 0
	} else if lonBin >= g.Cols {
		lonBin = g.Cols - 1
	}
	return latBin*g.Cols + lonBin
}

// geographicMembers groups the sampled ids that carry two strand copies by the map region
// their (Lat, Lon) falls into. Returns nil when no map is loaded (regions are undefined
// without map dimensions). Mirrors demeMembers, but keyed on geography rather than the
// Deme label.
func geographicMembers(model *core.Model, pop *core.Pop, ids []int, grid GeoGrid) map[int][]int {
	if len(model.Map) == 0 || len(model.Map[0]) == 0 {
		return nil
	}
	mapH := len(model.Map)
	mapW := len(model.Map[0])
	byRegion := make(map[int][]int)
	for _, id := range ids {
		c, ok := pop.Chromosomes[id]
		if !ok || len(c) < 2 || c[0] == nil || c[1] == nil {
			continue
		}
		lat := pop.IndData[id][individual.Lat]
		lon := pop.IndData[id][individual.Lon]
		r := grid.regionOf(lat, lon, mapH, mapW)
		byRegion[r] = append(byRegion[r], id)
	}
	return byRegion
}

// ComputeGeographicFst runs the differentiation pass over map regions rather than deme
// labels, reusing the Fst estimator core. Returns an empty result when no map is loaded or
// fewer than two regions hold >= 2 sampled diploids.
func ComputeGeographicFst(model *core.Model, pop *core.Pop, ids []int) *FstResult {
	grid := parseGeoGrid(model.StringParam("geographic_fst_grid", ""))
	byRegion := geographicMembers(model, pop, ids, grid)
	if byRegion == nil {
		fmt.Printf("Geographic Fst: no map loaded (track_map=1 required); skipped\n")
		return &FstResult{DemeSizes: map[int]int{}}
	}
	return computeFstFromPartition(model, pop, byRegion, "region")
}

// SaveGeographicFst writes the pairwise geographic-Fst table (both estimators) plus a
// global-θ row to a CSV, labelling the population axis as regions. A no-op (nil error)
// when the result has no pairs (fewer than two eligible regions). Mirrors SaveFst.
func SaveGeographicFst(model *core.Model, result *FstResult) error {
	if result == nil || len(result.Pairs) == 0 {
		return nil
	}

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	path := fmt.Sprintf("%s/%s_geographic_fst_run%d_year%d.csv",
		model.ResultsDir, model.ModelName, run, year)

	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create geographic Fst file: %v", err)
	}
	defer file.Close()

	w := csv.NewWriter(file)
	defer w.Flush()

	if err := w.Write([]string{"RegionA", "RegionB", "Fst_Hudson", "Fst_WC", "NumSites"}); err != nil {
		return err
	}
	for _, p := range result.Pairs {
		row := []string{
			fmt.Sprintf("%d", p.DemeA),
			fmt.Sprintf("%d", p.DemeB),
			fmt.Sprintf("%.6f", p.Hudson),
			fmt.Sprintf("%.6f", p.WC),
			fmt.Sprintf("%d", p.NumSites),
		}
		if err := w.Write(row); err != nil {
			return err
		}
	}
	// Global Weir & Cockerham θ across all eligible regions (RegionB = -1 sentinel).
	if err := w.Write([]string{
		"ALL", "-1", "",
		fmt.Sprintf("%.6f", result.GlobalWC),
		fmt.Sprintf("%d", result.NumSites),
	}); err != nil {
		return err
	}

	fmt.Printf("Geographic Fst: %d regions, %d pairs, %d polymorphic sites (global W&C θ = %.4f) → %s\n",
		len(result.Demes), len(result.Pairs), result.NumSites, result.GlobalWC, path)
	return nil
}
