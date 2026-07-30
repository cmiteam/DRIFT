package analysis

import (
	"testing"

	"drift/pkg/core"
	"drift/pkg/individual"
)

// makeGeoTestModel returns an IBD-style 128-bit model (arms 0:[0,64], 1:[64,64]) with a
// 4x4 all-Land map and the given geographic_fst_grid spec.
func makeGeoTestModel(grid string) *core.Model {
	m := makeIBDTestModel()
	m.ModelName = "GeoFstTest"
	m.StringParams = map[string]string{"geographic_fst_grid": grid}
	m.Map = [][]int{
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
	}
	return m
}

// geoInd builds an all-Deme-0 individual placed at (lat, lon) with a 128-bit chromosome.
func geoInd(lat, lon int, c0w0, c0w1, c1w0, c1w1 uint64) ([]int, [][]uint64) {
	data := individual.MakeIndData()
	data[individual.Deme] = 0 // deliberately one deme: geography is the ONLY axis
	data[individual.Lat] = lat
	data[individual.Lon] = lon
	chr := [][]uint64{{c0w0, c0w1}, {c1w0, c1w1}}
	return data, chr
}

func TestParseGeoGrid(t *testing.T) {
	cases := []struct {
		spec       string
		rows, cols int
	}{
		{"2x2", 2, 2},
		{"4x3", 4, 3},
		{"5", 5, 5},       // bare n => n×n
		{"", 2, 2},        // default
		{"garbage", 2, 2}, // malformed => default
		{"3x0", 2, 2},     // non-positive cols => default
		{"0x3", 2, 2},     // non-positive rows => default
		{" 3 X 2 ", 3, 2}, // whitespace + uppercase X tolerated
	}
	for _, c := range cases {
		g := parseGeoGrid(c.spec)
		if g.Rows != c.rows || g.Cols != c.cols {
			t.Errorf("parseGeoGrid(%q) = %dx%d, want %dx%d", c.spec, g.Rows, g.Cols, c.rows, c.cols)
		}
	}
}

func TestRegionOf(t *testing.T) {
	g := GeoGrid{Rows: 1, Cols: 2} // left/right split on a width-4 map
	// lon 0,1 => region 0; lon 2,3 => region 1; lat is irrelevant (1 row).
	cases := []struct {
		lat, lon, want int
	}{
		{0, 0, 0}, {3, 1, 0}, {0, 2, 1}, {2, 3, 1},
		{0, 99, 1}, // clamped past the right edge
		{0, -5, 0}, // clamped past the left edge
	}
	for _, c := range cases {
		if got := g.regionOf(c.lat, c.lon, 4, 4); got != c.want {
			t.Errorf("regionOf(lat=%d,lon=%d) = %d, want %d", c.lat, c.lon, got, c.want)
		}
	}
	// A 2x2 grid ids run row-major: latBin*Cols + lonBin.
	g2 := GeoGrid{Rows: 2, Cols: 2}
	if r := g2.regionOf(3, 3, 4, 4); r != 3 { // bottom-right
		t.Errorf("2x2 bottom-right region = %d, want 3", r)
	}
	if r := g2.regionOf(0, 3, 4, 4); r != 1 { // top-right
		t.Errorf("2x2 top-right region = %d, want 1", r)
	}
}

func TestGeographicMembers_BinsByGeography(t *testing.T) {
	model := makeGeoTestModel("1x2")
	grid := GeoGrid{Rows: 1, Cols: 2}

	dL, cL := geoInd(0, 0, bit0, 0, bit0, 0) // region 0
	dR, cR := geoInd(0, 3, 0, 0, 0, 0)       // region 1
	pop, ids := buildFstPop([][2]interface{}{{dL, cL}, {dR, cR}})

	// An individual with no chromosomes must be dropped.
	noChr := individual.MakeIndData()
	noChr[individual.Lon] = 0
	pop.IndData[99] = noChr
	ids = append(ids, 99)

	byRegion := geographicMembers(model, pop, ids, grid)
	if len(byRegion[0]) != 1 || len(byRegion[1]) != 1 {
		t.Fatalf("region buckets = %v, want one member each in regions 0 and 1", byRegion)
	}
	if byRegion[0][0] != 1 || byRegion[1][0] != 2 {
		t.Errorf("wrong members: %v", byRegion)
	}
}

// TestComputeGeographicFst_SplitByGeographyNotDeme is the core payoff test: a population
// that is a SINGLE deme (so deme-Fst sees no structure) but is split across two map
// regions with fully-differentiated genotypes yields geographic Fst ≈ 1. This is exactly
// the structure a movement barrier or habitat gradient creates — invisible to deme-Fst,
// measured by geographic Fst.
func TestComputeGeographicFst_SplitByGeographyNotDeme(t *testing.T) {
	model := makeGeoTestModel("1x2")

	// Left region (lon 0-1): homozygous derived. Right region (lon 2-3): homozygous
	// ancestral. All individuals are Deme 0.
	la, ca := geoInd(0, 0, bit0, 0, bit0, 0)
	lb, cb := geoInd(1, 1, bit0, 0, bit0, 0)
	ra, cra := geoInd(2, 2, 0, 0, 0, 0)
	rb, crb := geoInd(3, 3, 0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{{la, ca}, {lb, cb}, {ra, cra}, {rb, crb}})

	// Deme-Fst sees a single deme => empty (no structure).
	if demeRes := ComputeFst(model, pop, ids); len(demeRes.Pairs) != 0 {
		t.Fatalf("deme-Fst should be empty (one deme), got %d pairs", len(demeRes.Pairs))
	}

	// Geographic Fst recovers the complete differentiation across the two regions.
	geoRes := ComputeGeographicFst(model, pop, ids)
	if geoRes.NumSites != 1 {
		t.Fatalf("geographic Fst NumSites = %d, want 1", geoRes.NumSites)
	}
	if len(geoRes.Pairs) != 1 {
		t.Fatalf("expected 1 region pair, got %d", len(geoRes.Pairs))
	}
	p := geoRes.Pairs[0]
	assertClose(t, "Hudson", p.Hudson, 1.0)
	assertClose(t, "WC", p.WC, 1.0)
	assertClose(t, "GlobalWC", geoRes.GlobalWC, 1.0)
}

// TestComputeGeographicFst_Homogeneous: same genotypes in both regions => no
// differentiation (W&C θ = 0).
func TestComputeGeographicFst_Homogeneous(t *testing.T) {
	model := makeGeoTestModel("1x2")
	la, ca := geoInd(0, 0, bit0, 0, 0, 0) // heterozygous
	lb, cb := geoInd(1, 1, bit0, 0, 0, 0)
	ra, cra := geoInd(2, 2, bit0, 0, 0, 0)
	rb, crb := geoInd(3, 3, bit0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{{la, ca}, {lb, cb}, {ra, cra}, {rb, crb}})

	geoRes := ComputeGeographicFst(model, pop, ids)
	assertClose(t, "WC", geoRes.Pairs[0].WC, 0.0)
}

// TestComputeGeographicFst_NoMapGuard: with no map loaded, geographic regions are
// undefined => empty result, no panic.
func TestComputeGeographicFst_NoMapGuard(t *testing.T) {
	model := makeGeoTestModel("2x2")
	model.Map = nil
	da, cra := geoInd(0, 0, bit0, 0, bit0, 0)
	db, crb := geoInd(0, 3, 0, 0, 0, 0)
	pop, ids := buildFstPop([][2]interface{}{{da, cra}, {db, crb}})

	res := ComputeGeographicFst(model, pop, ids)
	if len(res.Pairs) != 0 {
		t.Fatalf("expected empty result with no map, got %d pairs", len(res.Pairs))
	}
}
