package simulation

import (
	"sort"
	"testing"

	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
)

// buildDemePop builds a Pop whose individuals carry the given deme labels (ids 1..n).
func buildDemePop(demes []int) *core.Pop {
	pop := &core.Pop{IndData: map[int][]int{}}
	for i, d := range demes {
		data := individual.MakeIndData()
		data[individual.Deme] = d
		pop.IndData[i+1] = data
	}
	return pop
}

// demeCounts tallies how many individuals sit in each deme.
func demeCounts(pop *core.Pop) map[int]int {
	counts := map[int]int{}
	for _, data := range pop.IndData {
		counts[data[individual.Deme]]++
	}
	return counts
}

// TestMigrateDemesMatrix_AllMoveDirectional: a rate-1.0 edge 0>1 moves every resident of
// deme 0 into deme 1; deme 1 has no outgoing edge so it stays. Distinguishable outcome.
func TestMigrateDemesMatrix_AllMoveDirectional(t *testing.T) {
	utils.SeedRNG(4242)
	pop := buildDemePop([]int{0, 0, 0, 1, 1})
	mm := core.MigrationMatrix{Active: true, Out: map[int][]core.MigrationEdge{
		0: {{Dest: 1, Rate: 1.0}},
	}}
	migrateDemesMatrix(pop, &mm)
	counts := demeCounts(pop)
	if counts[0] != 0 || counts[1] != 5 {
		t.Fatalf("after 0>1:1.0, counts = %v, want deme0=0 deme1=5", counts)
	}
}

// TestMigrateDemesMatrix_SteppingStone: 0>1:1.0 and 1>2:1.0 shift each deme one step
// downstream in a single application (0→1, 1→2, 2 stays). Because moves are snapshotted
// before applying, an individual that was in deme 0 does not cascade 0→1→2 in one year.
func TestMigrateDemesMatrix_SteppingStone(t *testing.T) {
	utils.SeedRNG(1)
	pop := buildDemePop([]int{0, 0, 1, 1, 1, 2}) // 2 in d0, 3 in d1, 1 in d2
	mm := core.MigrationMatrix{Active: true, Out: map[int][]core.MigrationEdge{
		0: {{Dest: 1, Rate: 1.0}},
		1: {{Dest: 2, Rate: 1.0}},
	}}
	migrateDemesMatrix(pop, &mm)
	counts := demeCounts(pop)
	// d0 emptied → d1; original d1 → d2; original d2 stays in d2.
	if counts[0] != 0 || counts[1] != 2 || counts[2] != 4 {
		t.Fatalf("stepping-stone counts = %v, want d0=0 d1=2 d2=4", counts)
	}
}

// TestMigrateDemesMatrix_Asymmetric: with 0>1:1.0 but no 1>0 edge, all flow is one-way —
// everyone ends in deme 1, none return. Demonstrates the asymmetry the scalar island
// model (uniform, symmetric) cannot express.
func TestMigrateDemesMatrix_Asymmetric(t *testing.T) {
	utils.SeedRNG(7)
	pop := buildDemePop([]int{0, 0, 0, 1, 1})
	mm := core.MigrationMatrix{Active: true, Out: map[int][]core.MigrationEdge{
		0: {{Dest: 1, Rate: 1.0}},
		// deme 1: no outgoing edge (sink)
	}}
	migrateDemesMatrix(pop, &mm)
	if c := demeCounts(pop); c[0] != 0 || c[1] != 5 {
		t.Fatalf("asymmetric sink counts = %v, want all in deme 1", c)
	}
}

// TestMigrateDemesMatrix_PartialRateDeterministic: a fractional rate moves some but not
// all, and the same seed reproduces the exact partition (determinism under a fixed seed).
func TestMigrateDemesMatrix_PartialRateDeterministic(t *testing.T) {
	mm := core.MigrationMatrix{Active: true, Out: map[int][]core.MigrationEdge{
		0: {{Dest: 1, Rate: 0.5}},
	}}
	run := func() []int {
		utils.SeedRNG(2024)
		pop := buildDemePop([]int{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}) // 10 in deme 0
		migrateDemesMatrix(pop, &mm)
		ids := make([]int, 0)
		for id, data := range pop.IndData {
			if data[individual.Deme] == 1 {
				ids = append(ids, id)
			}
		}
		sort.Ints(ids)
		return ids
	}
	first := run()
	if len(first) == 0 || len(first) == 10 {
		t.Fatalf("partial rate moved %d/10; expected a strict subset", len(first))
	}
	second := run()
	if len(first) != len(second) {
		t.Fatalf("non-deterministic: %d vs %d moved", len(first), len(second))
	}
	for i := range first {
		if first[i] != second[i] {
			t.Fatalf("non-deterministic movers: %v vs %v", first, second)
		}
	}
}

// TestMigrateDemesMatrix_OneDrawPerIndividual: exactly one RandFloat64 is consumed per
// individual regardless of how many outgoing edges its deme has (a deme with no edges
// still draws), so the RNG stream stays aligned. Verified by advancing a reference stream
// by len(pop) draws and confirming the post-migration RNG state matches.
func TestMigrateDemesMatrix_OneDrawPerIndividual(t *testing.T) {
	pop := buildDemePop([]int{0, 1, 2, 0, 1}) // demes 1,2 have no outgoing edges
	mm := core.MigrationMatrix{Active: true, Out: map[int][]core.MigrationEdge{
		0: {{Dest: 1, Rate: 0.3}},
	}}

	utils.SeedRNG(99)
	migrateDemesMatrix(pop, &mm)
	afterMigration := utils.RandFloat64()

	// Reference: 5 individuals ⇒ 5 draws, then the 6th draw must match.
	utils.SeedRNG(99)
	for i := 0; i < len(pop.IndData); i++ {
		utils.RandFloat64()
	}
	reference := utils.RandFloat64()

	if afterMigration != reference {
		t.Fatalf("migration consumed the wrong number of draws: next value %v != %v", afterMigration, reference)
	}
}
