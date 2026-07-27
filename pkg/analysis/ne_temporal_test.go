package analysis

import (
	"bytes"
	"math"
	"os"
	"path/filepath"
	"testing"

	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
)

// neInd builds one individual (with a deme label) and a 128-bit chromosome (two
// 64-bit words per strand), setting the given (strand, bit) positions to 1.
func neInd(deme int, set ...[2]int) ([]int, [][]uint64) {
	data := individual.MakeIndData()
	data[individual.Deme] = deme
	chr := [][]uint64{{0, 0}, {0, 0}}
	for _, sb := range set {
		strand, bit := sb[0], sb[1]
		chr[strand][bit/64] |= uint64(1) << (uint(bit) % 64)
	}
	return data, chr
}

// makeNeModel is makeIBDTestModel plus the parameters CaptureNe reads.
func makeNeModel(t *testing.T) *core.Model {
	m := makeIBDTestModel()
	m.Parameters = map[string]float64{
		"num_demes":       1,
		"ne_sample_size":  0,
		"generation_time": 1,
		"ne_cM_per_bit":   0,
		"ld_max_pairs":    10000,
	}
	m.ResultsDir = t.TempDir()
	m.ModelName = "NeTest"
	return m
}

// TestTemporalNe_RatioOfSums checks the Waples plan-II estimator against an
// explicit two-locus hand computation, confirming Fc is a ratio of sums (not a mean
// of per-locus ratios) and the 1/n0 + 1/nt sampling correction.
//
// prev: site 0 freq 0.5, site 64 freq 0.4, haploid n0 = 20.
// curr (built below): site 0 derived 5/20 = 0.25, site 64 derived 12/20 = 0.6, nt = 20.
func TestTemporalNe_RatioOfSums(t *testing.T) {
	// 10 diploids giving derived counts 5 at bit 0 and 12 at bit 64.
	d1, c1 := neInd(0, [2]int{0, 0}, [2]int{1, 0}, [2]int{0, 64}, [2]int{1, 64}) // bit0 hom, bit64 hom
	d2, c2 := neInd(0, [2]int{0, 0}, [2]int{1, 0}, [2]int{0, 64}, [2]int{1, 64})
	d3, c3 := neInd(0, [2]int{0, 0}, [2]int{0, 64}, [2]int{1, 64}) // bit0 het, bit64 hom
	d4, c4 := neInd(0, [2]int{0, 64}, [2]int{1, 64})               // bit64 hom
	d5, c5 := neInd(0, [2]int{0, 64}, [2]int{1, 64})
	d6, c6 := neInd(0, [2]int{0, 64}, [2]int{1, 64})
	d7, c7 := neInd(0)
	d8, c8 := neInd(0)
	d9, c9 := neInd(0)
	d10, c10 := neInd(0)
	pop, ids := buildFstPop([][2]interface{}{
		{d1, c1}, {d2, c2}, {d3, c3}, {d4, c4}, {d5, c5},
		{d6, c6}, {d7, c7}, {d8, c8}, {d9, c9}, {d10, c10},
	})

	prev := &core.FreqSnapshot{
		Year:     0,
		Groups:   map[int]map[int]float64{-1: {0: 0.5, 64: 0.4}},
		HaploidN: map[int]int{-1: 20},
	}

	got, used := temporalNe(pop, prev, -1, ids, 1.0)
	if used != 2 {
		t.Fatalf("used sites = %d, want 2", used)
	}

	// Explicit hand computation with the observed current frequencies.
	x0, y0 := 0.5, 0.25
	x1, y1 := 0.4, 0.6
	num := (x0-y0)*(x0-y0) + (x1-y1)*(x1-y1)
	den := ((x0+y0)/2 - x0*y0) + ((x1+y1)/2 - x1*y1)
	fc := num / den
	corrected := fc - 1.0/20 - 1.0/20
	want := 1.0 / (2 * corrected)
	assertClose(t, "TemporalNe", got, want)
}

// TestTemporalNe_Unresolved: with tiny samples the sampling correction exceeds Fc,
// so the drift signal is below the noise floor and Ne is reported as +Inf.
func TestTemporalNe_Unresolved(t *testing.T) {
	// nt = 4 (2 diploids), one het at bit 0 => derived 1/4 = 0.25.
	d1, c1 := neInd(0, [2]int{0, 0})
	d2, c2 := neInd(0)
	pop, ids := buildFstPop([][2]interface{}{{d1, c1}, {d2, c2}})

	prev := &core.FreqSnapshot{
		Year:     0,
		Groups:   map[int]map[int]float64{-1: {0: 0.5}},
		HaploidN: map[int]int{-1: 4},
	}
	got, _ := temporalNe(pop, prev, -1, ids, 1.0)
	if !math.IsInf(got, 1) {
		t.Fatalf("TemporalNe = %v, want +Inf (unresolved)", got)
	}
}

// TestTemporalNe_Distinguishable: a bigger allele-frequency change over the same
// interval must yield a SMALLER inferred Ne — the core signal a bottleneck produces.
func TestTemporalNe_Distinguishable(t *testing.T) {
	prev := &core.FreqSnapshot{
		Year:     0,
		Groups:   map[int]map[int]float64{-1: {0: 0.5}},
		HaploidN: map[int]int{-1: 20},
	}

	// Medium change: bit-0 derived 5/20 = 0.25.
	med := make([][2]interface{}, 0, 10)
	dm1, cm1 := neInd(0, [2]int{0, 0}, [2]int{1, 0})
	dm2, cm2 := neInd(0, [2]int{0, 0}, [2]int{1, 0})
	dm3, cm3 := neInd(0, [2]int{0, 0})
	med = append(med, [2]interface{}{dm1, cm1}, [2]interface{}{dm2, cm2}, [2]interface{}{dm3, cm3})
	for i := 0; i < 7; i++ {
		d, c := neInd(0)
		med = append(med, [2]interface{}{d, c})
	}
	popMed, idsMed := buildFstPop(med)
	neMed, _ := temporalNe(popMed, prev, -1, idsMed, 1.0)

	// Big change: bit-0 derived 1/20 = 0.05.
	big := make([][2]interface{}, 0, 10)
	db1, cb1 := neInd(0, [2]int{0, 0})
	big = append(big, [2]interface{}{db1, cb1})
	for i := 0; i < 9; i++ {
		d, c := neInd(0)
		big = append(big, [2]interface{}{d, c})
	}
	popBig, idsBig := buildFstPop(big)
	neBig, _ := temporalNe(popBig, prev, -1, idsBig, 1.0)

	if math.IsInf(neMed, 0) || math.IsInf(neBig, 0) {
		t.Fatalf("both estimates should resolve: neMed=%v neBig=%v", neMed, neBig)
	}
	if !(neBig < neMed) {
		t.Errorf("expected larger drift => smaller Ne: neBig=%.4f, neMed=%.4f", neBig, neMed)
	}
}

// TestCaptureNe_EndToEnd drives CaptureNe over two windows on a synthetic pop:
// the first seeds the reference snapshot; after mutating frequencies the second
// produces a resolved temporal estimate. SaveNeTimeSeries then writes the series
// and summary CSVs.
func TestCaptureNe_EndToEnd(t *testing.T) {
	model := makeNeModel(t)

	// 12 diploids, all heterozygous at bit 0 (freq 0.5) plus varied bit-64 state so
	// there are multiple segregating sites.
	specs := make([][2]interface{}, 0, 12)
	for i := 0; i < 12; i++ {
		var d []int
		var c [][]uint64
		if i%2 == 0 {
			d, c = neInd(0, [2]int{0, 0}, [2]int{0, 64}) // het at bit0, het at bit64
		} else {
			d, c = neInd(0, [2]int{0, 0}) // het at bit0 only
		}
		specs = append(specs, [2]interface{}{d, c})
	}
	pop, _ := buildFstPop(specs)

	// Window 1 (reference).
	model.FreeParameters["year"] = 0
	CaptureNe(model, pop)
	if pop.NePrevSnap == nil {
		t.Fatal("CaptureNe did not seed the previous snapshot")
	}

	// Drift bit-0 toward fixation between windows: clear the derived allele on a few
	// individuals' second strand so the frequency changes.
	changed := 0
	for id := range pop.Chromosomes {
		if changed >= 6 {
			break
		}
		pop.Chromosomes[id][0][0] = 0 // clear bit 0 on strand 0
		changed++
	}

	// Window 2.
	model.FreeParameters["year"] = 10
	CaptureNe(model, pop)

	if len(pop.NeHistory) != 2 {
		t.Fatalf("NeHistory has %d windows, want 2", len(pop.NeHistory))
	}
	w2 := pop.NeHistory[10]
	if w2 == nil || w2.IntervalGen != 10 {
		t.Fatalf("window 2 interval = %v, want 10", w2)
	}
	if w2.CensusN != 12 {
		t.Errorf("CensusN = %d, want 12", w2.CensusN)
	}

	if err := SaveNeTimeSeries(model, pop); err != nil {
		t.Fatalf("SaveNeTimeSeries: %v", err)
	}
	for _, suffix := range []string{"_ne_timeseries_run1.csv", "_ne_summary_run1.csv"} {
		p := filepath.Join(model.ResultsDir, model.ModelName+suffix)
		if _, err := os.Stat(p); err != nil {
			t.Errorf("expected output %s: %v", p, err)
		}
	}
}

// TestCaptureNe_PerDeme: with the island model active, per-deme temporal Ne is
// populated for each deme (composition with §6b).
func TestCaptureNe_PerDeme(t *testing.T) {
	model := makeNeModel(t)
	model.Parameters["num_demes"] = 2

	specs := make([][2]interface{}, 0, 8)
	for i := 0; i < 8; i++ {
		deme := i % 2
		d, c := neInd(deme, [2]int{0, 0}) // het at bit 0 in every individual
		specs = append(specs, [2]interface{}{d, c})
	}
	pop, _ := buildFstPop(specs)

	model.FreeParameters["year"] = 0
	CaptureNe(model, pop)
	model.FreeParameters["year"] = 5
	CaptureNe(model, pop)

	w := pop.NeHistory[5]
	if w == nil {
		t.Fatal("missing window at year 5")
	}
	if len(w.DemeNe) != 2 {
		t.Fatalf("DemeNe has %d demes, want 2", len(w.DemeNe))
	}
	if _, ok := w.DemeNe[0]; !ok {
		t.Error("missing per-deme Ne for deme 0")
	}
	if _, ok := w.DemeNe[1]; !ok {
		t.Error("missing per-deme Ne for deme 1")
	}
}

// TestCaptureNe_NoRNG guards byte-identity: with ne_sample_size = 0 the capture is a
// deterministic full-population scan that must not touch the shared RNG stream, so a
// track_Ne run stays byte-identical to a track_Ne-off run.
func TestCaptureNe_NoRNG(t *testing.T) {
	model := makeNeModel(t)
	d1, c1 := neInd(0, [2]int{0, 0})
	d2, c2 := neInd(0, [2]int{0, 0})
	d3, c3 := neInd(0)
	pop, _ := buildFstPop([][2]interface{}{{d1, c1}, {d2, c2}, {d3, c3}})

	utils.SeedRNG(20260727)
	before := utils.RNGState()
	model.FreeParameters["year"] = 0
	CaptureNe(model, pop)
	model.FreeParameters["year"] = 10
	CaptureNe(model, pop)
	after := utils.RNGState()

	if !bytes.Equal(before, after) {
		t.Error("CaptureNe consumed RNG with ne_sample_size = 0; must be RNG-free")
	}
}
