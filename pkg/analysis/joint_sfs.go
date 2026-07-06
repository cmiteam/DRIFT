package analysis

// Joint / multi-population site frequency spectrum (2D / 3D) — roadmap §6b, the
// last differentiation-statistics item. The single-population SFS (sfs.go) counts,
// for each derived-allele count k, how many sites segregate at that count. The
// JOINT SFS generalizes this to a matrix indexed by the derived count in EACH of
// two (2D) or three (3D) demes simultaneously. It is the object the out-of-Africa
// demographic model is actually fit to — dadi, fastsimcoal2 and momi all optimize
// split times, sizes and migration against the joint SFS of African / European /
// Asian samples — so emitting it lets DRIFT be scored with the field's own
// inference tools.
//
// POPULATIONS. As elsewhere in §6b, "populations" are the island-model demes; this
// reuses the fst.go substrate (demeMembers, countSite) and the buildContigs/bitAt
// genome scan. A cell (i,j[,k]) counts the sites at which the first axis deme has i
// derived alleles in its sample, the second j, the third k, out of that deme's
// sampled haplotypes (n = 2 x sampled diploids). As in sfs.go, only sites
// SEGREGATING within the combined sample of the axis demes are counted (the
// all-ancestral and all-derived corners carry no information and would dwarf the
// spectrum).
//
// WHICH DEMES FORM THE AXES. The joint_sfs_demes param names 2 or 3 demes
// (e.g. "0,1" or "0,1,2"); when empty the first three eligible demes are used (or
// both, if only two exist). Wired into drift.go behind track_joint_sfs (default
// off); requires track_DNA and num_demes >= 2. Sample size via
// joint_sfs_sample_size (0 = all).
//
// OUTPUT. A long-format CSV — one row per non-empty cell, columns being the derived
// count on each axis plus the site count. Long format loads directly into
// numpy/dadi (reshape by the recorded axis sizes) and generalizes cleanly from 2D
// to 3D.

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
)

// JointSFSResult is the output of a joint-SFS pass.
type JointSFSResult struct {
	Demes    []int          // the demes forming the axes, in order
	Sizes    []int          // sampled haplotypes per axis deme (axis length = Sizes[i]+1)
	Dim      int            // 2 or 3
	Cells    map[[3]int]int // (i,j,k) -> site count; k is 0 (unused) for 2D
	NumSites int
}

// ComputeJointSFS builds the 2D/3D joint SFS over the given sampled individual IDs
// (see SampleIDs). Returns an empty result (Dim 0) when fewer than two axis demes
// are available. The axis demes are chosen from the joint_sfs_demes param, or the
// first eligible demes when it is empty.
func ComputeJointSFS(model *core.Model, pop *core.Pop, ids []int) *JointSFSResult {
	byDeme := demeMembers(pop, ids)

	// Eligible demes need at least one sampled diploid (n >= 2) to form an axis.
	elig := make([]int, 0, len(byDeme))
	for d, members := range byDeme {
		if len(members) >= 1 {
			elig = append(elig, d)
		}
	}
	sort.Ints(elig)

	res := &JointSFSResult{Cells: map[[3]int]int{}}

	axisDemes, warn := planJointSFS(strings.TrimSpace(model.StringParam("joint_sfs_demes", "")), elig)
	if warn != "" {
		fmt.Printf("joint SFS: %s\n", warn)
	}
	if len(axisDemes) < 2 {
		fmt.Printf("joint SFS: need >= 2 eligible demes (have %d); set num_demes and run long enough for demes to persist\n", len(elig))
		return res
	}

	res.Demes = axisDemes
	res.Dim = len(axisDemes)
	res.Sizes = make([]int, len(axisDemes))
	totalAlleles := 0
	for i, d := range axisDemes {
		res.Sizes[i] = 2 * len(byDeme[d])
		totalAlleles += res.Sizes[i]
	}

	contigs := buildContigs(model)
	for _, ct := range contigs {
		for _, arm := range ct.arms {
			start, length := arm[0], arm[1]
			for bit := start; bit < start+length; bit++ {
				var key [3]int
				combined := 0
				for i, d := range axisDemes {
					dc := countSite(pop, byDeme[d], bit).derived
					key[i] = dc
					combined += dc
				}
				// Skip sites monomorphic across the combined axis-deme sample.
				if combined == 0 || combined == totalAlleles {
					continue
				}
				res.Cells[key]++
				res.NumSites++
			}
		}
	}
	return res
}

// planJointSFS picks the 2 or 3 axis demes. An explicit joint_sfs_demes list wins
// (validated against the eligible set); otherwise the first three eligible demes
// (or both, if only two) are used.
func planJointSFS(spec string, elig []int) (axes []int, warn string) {
	eligible := make(map[int]bool, len(elig))
	for _, d := range elig {
		eligible[d] = true
	}
	if spec != "" {
		list, err := parseDemeList(spec)
		if err != nil {
			return nil, fmt.Sprintf("joint_sfs_demes %q has a non-integer deme — falling back to defaults", spec)
		}
		if len(list) < 2 || len(list) > 3 {
			return nil, fmt.Sprintf("joint_sfs_demes %q must name 2 or 3 demes — falling back to defaults", spec)
		}
		for _, d := range list {
			if !eligible[d] {
				return nil, fmt.Sprintf("joint_sfs_demes %q names deme %d which has no sampled diploids — falling back to defaults", spec, d)
			}
		}
		return list, ""
	}
	if len(elig) >= 3 {
		return elig[:3], ""
	}
	if len(elig) == 2 {
		return elig[:2], ""
	}
	return nil, ""
}

// SaveJointSFS writes the joint spectrum as a long-format CSV (one row per non-empty
// cell) plus a one-line summary. A no-op (nil) when the result has no axes or no
// segregating sites.
func SaveJointSFS(model *core.Model, result *JointSFSResult) error {
	if result == nil || result.Dim < 2 || result.NumSites == 0 {
		return nil
	}

	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	path := fmt.Sprintf("%s/%s_joint_sfs_run%d_year%d.csv",
		model.ResultsDir, model.ModelName, run, year)

	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create joint SFS file: %v", err)
	}
	defer file.Close()

	w := csv.NewWriter(file)
	defer w.Flush()

	// Header: one derived-count column per axis deme, named by deme id, then Count.
	header := make([]string, 0, result.Dim+1)
	for _, d := range result.Demes {
		header = append(header, fmt.Sprintf("Derived_deme%d", d))
	}
	header = append(header, "Count")
	if err := w.Write(header); err != nil {
		return err
	}

	// Deterministic cell order (lexicographic by the axis indices).
	keys := make([][3]int, 0, len(result.Cells))
	for k := range result.Cells {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(a, b int) bool {
		for i := 0; i < 3; i++ {
			if keys[a][i] != keys[b][i] {
				return keys[a][i] < keys[b][i]
			}
		}
		return false
	})

	for _, k := range keys {
		row := make([]string, 0, result.Dim+1)
		for i := 0; i < result.Dim; i++ {
			row = append(row, strconv.Itoa(k[i]))
		}
		row = append(row, strconv.Itoa(result.Cells[k]))
		if err := w.Write(row); err != nil {
			return err
		}
	}

	sizes := make([]string, len(result.Sizes))
	for i, s := range result.Sizes {
		sizes[i] = strconv.Itoa(s)
	}
	fmt.Printf("joint SFS (%dD): demes %v, axis sizes [%s], %d segregating sites, %d non-empty cells → %s\n",
		result.Dim, result.Demes, strings.Join(sizes, "x"), result.NumSites, len(result.Cells), path)
	return nil
}
