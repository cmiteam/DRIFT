package utils

import (
	"testing"

	"drift/pkg/core"
)

// countModel returns a tiled-arms model with mutation_count_model set (and
// linkage_model=arm so the mask polarity matches meiosis — irrelevant to Count,
// but keeps the fixture realistic).
func countModel(mode string) *core.Model {
	m := tiledArmsModel("arm")
	if mode != "" {
		m.StringParams["mutation_count_model"] = mode
	}
	return m
}

// TestMutationCountLegacyAsymmetry documents the historical bug: a strand-0
// mutation inherited to a child does NOT bump its reference count under the
// default "legacy" model, while a strand-1 mutation does.
func TestMutationCountLegacyAsymmetry(t *testing.T) {
	// Parent (id 1) carries mutation 3 on strand 0 and mutation 5 on strand 1, both
	// at positions that the all-ones mask keeps on this gamete under arm polarity
	// (arm: keep0 where mask==1, keep1 where mask==0). Use a mask that keeps BOTH.
	// mask bit at pos 3 (=1) keeps strand0 under arm (keep0=1); pos 5 (=0) keeps
	// strand1 under arm (keep1=0). So set bit 3, clear bit 5.
	mask := []uint64{1 << 3} // bit 3 set, bit 5 clear
	pop := inheritPop([]int{3}, []int{5})
	InheritMutations(countModel("legacy"), pop, mask, 1, 2, 0)

	// Both were inherited to child 2 strand 0.
	got := pop.IndMutations[2][0]
	if len(got) != 2 {
		t.Fatalf("child inherited %v, want both 3 and 5", got)
	}
	// Legacy: strand-0 mutation (3) Count NOT incremented; strand-1 mutation (5) IS.
	if c := pop.MutationPool[3].Count; c != 1 {
		t.Errorf("legacy: strand-0 mutation Count=%d, want 1 (not incremented — the bug)", c)
	}
	if c := pop.MutationPool[5].Count; c != 2 {
		t.Errorf("legacy: strand-1 mutation Count=%d, want 2", c)
	}
}

// TestMutationCountRefcountSymmetry: under "refcount" BOTH inherited copies bump
// Count, so it tracks the true living-copy count.
func TestMutationCountRefcountSymmetry(t *testing.T) {
	mask := []uint64{1 << 3} // keep strand0 mut 3 (bit set) and strand1 mut 5 (bit clear)
	pop := inheritPop([]int{3}, []int{5})
	InheritMutations(countModel("refcount"), pop, mask, 1, 2, 0)

	if len(pop.IndMutations[2][0]) != 2 {
		t.Fatalf("child inherited %v, want both", pop.IndMutations[2][0])
	}
	if c := pop.MutationPool[3].Count; c != 2 {
		t.Errorf("refcount: strand-0 mutation Count=%d, want 2 (now symmetric)", c)
	}
	if c := pop.MutationPool[5].Count; c != 2 {
		t.Errorf("refcount: strand-1 mutation Count=%d, want 2", c)
	}
}

// TestMutationCountRefcountSurvivesGC is the payoff: a lineage inherited to a
// child and then losing its ORIGINAL carrier survives under refcount (Count still
// counts the child's copy) but is prematurely garbage-collected under legacy
// (Count already under-counted, so the single decrement drops it to 0). This is
// exactly why fixed load is invisible under legacy.
func TestMutationCountRefcountSurvivesGC(t *testing.T) {
	for _, mode := range []string{"legacy", "refcount"} {
		mask := []uint64{1 << 3} // keep strand-0 mutation 3 on the child
		pop := inheritPop([]int{3}, []int{})
		// Parent (1) and child (2) both now carry mutation 3.
		InheritMutations(countModel(mode), pop, mask, 1, 2, 0)

		// Simulate the parent dying: decrement its one carried copy (mirrors death.go).
		mut := pop.MutationPool[3]
		mut.Count--
		if mut.Count <= 0 {
			delete(pop.MutationPool, 3)
		} else {
			pop.MutationPool[3] = mut
		}

		_, stillResident := pop.MutationPool[3]
		childCarries := len(pop.IndMutations[2][0]) == 1
		if !childCarries {
			t.Fatalf("%s: child should carry mutation 3", mode)
		}
		switch mode {
		case "legacy":
			if stillResident {
				t.Errorf("legacy: mutation 3 unexpectedly survived (bug is it should be GC'd while child carries it)")
			}
		case "refcount":
			if !stillResident {
				t.Errorf("refcount: mutation 3 was GC'd while the child still carries it")
			}
		}
	}
}
