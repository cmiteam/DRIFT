package methods

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"testing"
)

// straitMatingModel builds a 7x16 map with an OpenWater(3) strait down columns 9-11
// (land elsewhere) and the given movement_barriers spec, plus a distance-mating
// configuration. Used to drive the real MatingDistance module across the strait.
func straitMatingModel(barriers string) *core.Model {
	m := &core.Model{
		Parameters: map[string]float64{
			"max_mating_distance": 5,
			"track_map":           1,
		},
		StringParams:   map[string]string{"movement_barriers": barriers},
		FreeParameters: map[string]int{"year": 100},
		Map:            make([][]int, 7),
	}
	for lat := 0; lat < 7; lat++ {
		row := make([]int, 16)
		for lon := 0; lon < 16; lon++ {
			if lon >= 9 && lon <= 11 {
				row[lon] = 3 // water strait
			} else {
				row[lon] = 1 // land
			}
		}
		m.Map[lat] = row
	}
	return m
}

// TestMatingBarrierBlocksStrait drives the real MatingDistance module over a
// population of men on the left shore (col 8) with eligible women both on the same
// shore (col 7, Manhattan 1) and across the strait (col 12, Manhattan 4 — well within
// max_mating_distance=5). Without a barrier, men marry across the strait; with the
// Water strait a barrier, every cross-strait pairing is blocked and men can only
// marry the same-shore women. This exercises the mating-side crossing test through the
// real influence-grid module, the mating analogue of the Wander dispersal block.
func TestMatingBarrierBlocksStrait(t *testing.T) {
	const nEach = 30

	// crossStraitMarriages runs MatingDistance once and returns (cross-strait, same-shore)
	// marriage counts. Men live at col 8; a marriage is cross-strait iff the man's wife
	// lives at col > 11 (the right shore).
	run := func(barriers string) (cross, same int) {
		m := straitMatingModel(barriers)
		utils.SeedRNG(7)
		pop := &core.Pop{IndData: map[int][]int{}, Tracking: map[string]int{}}

		var men, women []int
		id := 1
		// Men on the left shore at (3, 8).
		for i := 0; i < nEach; i++ {
			d := individual.MakeIndData()
			d[individual.Sex] = 0
			d[individual.Lat], d[individual.Lon] = 3, 8
			pop.IndData[id] = d
			men = append(men, id)
			id++
		}
		// Same-shore women at (3, 7).
		var sameShore, crossShore []int
		for i := 0; i < nEach; i++ {
			d := individual.MakeIndData()
			d[individual.Sex] = 1
			d[individual.Lat], d[individual.Lon] = 3, 7
			pop.IndData[id] = d
			women = append(women, id)
			sameShore = append(sameShore, id)
			id++
		}
		// Cross-strait women at (3, 12).
		for i := 0; i < nEach; i++ {
			d := individual.MakeIndData()
			d[individual.Sex] = 1
			d[individual.Lat], d[individual.Lon] = 3, 12
			pop.IndData[id] = d
			women = append(women, id)
			crossShore = append(crossShore, id)
			id++
		}

		MatingDistance(m, pop, men, women)

		crossSet := make(map[int]bool)
		for _, w := range crossShore {
			crossSet[w] = true
		}
		sameSet := make(map[int]bool)
		for _, w := range sameShore {
			sameSet[w] = true
		}
		for _, manID := range men {
			wife := pop.IndData[manID][individual.MarriageState]
			if wife == -1 {
				continue
			}
			if crossSet[wife] {
				cross++
			} else if sameSet[wife] {
				same++
			}
		}
		return cross, same
	}

	offCross, offSame := run("")
	onCross, onSame := run("3")

	if offCross == 0 {
		t.Fatalf("without a barrier, expected some cross-strait marriages, got 0 (same=%d)", offSame)
	}
	if onCross != 0 {
		t.Errorf("with a Water barrier, expected zero cross-strait marriages, got %d", onCross)
	}
	if onSame == 0 {
		t.Errorf("with a barrier, men should still marry same-shore women, got 0")
	}
	t.Logf("cross-strait marriages: OFF=%d ON=%d ; same-shore: OFF=%d ON=%d", offCross, onCross, offSame, onSame)
}
