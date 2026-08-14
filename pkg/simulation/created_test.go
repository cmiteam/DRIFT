package simulation

import (
	"drift/pkg/analysis"
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"sort"
	"testing"
)

// Created founder alleles with germline heterogeneity (TMR4A.md W4), transmission half.
// The claim being tested is the objection in §0 made concrete: ONE created couple
// transmitting MORE THAN FOUR alleles at a locus, which strict diploidy cannot do — and
// doing it without disturbing recombination, the mask, or anything else.

const createdTestGenomeBits = 256

// createdModel is pedModel plus the created founder-allele model and one chromosome.
func createdModel(alleles int) *core.Model {
	m := pedModel(1)
	m.Parameters["track_DNA"] = 1
	m.Parameters["track_arg"] = 1
	m.Parameters["founder_alleles_per_locus"] = float64(alleles)
	m.StringParams["founder_allele_model"] = "created"
	m.StringParams["arg_loci"] = "0-255:8"
	m.ChromosomeArms = map[int]map[int][]int{1: {0: {0, 128}, 1: {128, 128}}}
	m.FreeParameters["genome_bits"] = createdTestGenomeBits
	return m
}

// createdPop seeds one couple with a created-allele pool. The allele sequences are drawn
// so that any two differ at roughly half their sites, which makes "which allele did this
// bit come from" answerable and keeps the meiosis check below from passing vacuously.
func createdPop(model *core.Model, alleles int) *core.Pop {
	pop := pedPop(model, 1) // founders 0 (male) and 1 (female), married
	words := (createdTestGenomeBits + 63) / 64
	pop.CreatedAlleleSeqs = make([][]uint64, alleles)
	for a := 0; a < alleles; a++ {
		seq := make([]uint64, words)
		for bit := 0; bit < createdTestGenomeBits; bit++ {
			if utils.RandFloat64() < 0.5 {
				seq[bit/64] |= 1 << (uint(bit) % 64)
			}
		}
		pop.CreatedAlleleSeqs[a] = seq
	}
	pop.AssignCreatedPools([]int{0, 1}, alleles)
	for _, id := range []int{0, 1} {
		pool := pop.CreatedPool(id)
		second := 0
		if len(pool) > 1 {
			second = 1
		}
		pop.Chromosomes[id] = [][]uint64{
			append([]uint64(nil), pop.CreatedAlleleSeqs[pool[0]]...),
			append([]uint64(nil), pop.CreatedAlleleSeqs[pool[second]]...),
		}
		pop.IndData[id][individual.AlleleCount] =
			utils.CountSetBits(pop.Chromosomes[id][0]) + utils.CountSetBits(pop.Chromosomes[id][1])
	}
	return pop
}

// runCreatedBirths advances births only — no deaths, so the founding couple survives to
// produce many gametes, which is what exposes germline heterogeneity.
func runCreatedBirths(model *core.Model, pop *core.Pop, years int) {
	for y := 1; y <= years; y++ {
		model.FreeParameters["year"] = y
		birthStandard(model, pop)
	}
}

// THE LOAD-BEARING TEST. A single created couple must transmit more than four distinct
// alleles at a locus. Under ordinary diploidy four is a hard ceiling — two parents, two
// strands each — and TMR4A's "4" assumes exactly that.
func TestCreatedCoupleTransmitsMoreThanFourAlleles(t *testing.T) {
	const alleles = 10
	utils.SeedRNG(20260813)
	model := createdModel(alleles)
	pop := createdPop(model, alleles)
	runCreatedBirths(model, pop, 40)

	loci := core.EnsureARGLoci(model)
	if len(pop.IndData) <= 2 {
		t.Fatal("no children were born; the test proved nothing")
	}

	transmitted := map[int]bool{}
	perParent := map[int]map[int]bool{0: {}, 1: {}}
	children := 0
	for id := range pop.IndData {
		if id < 2 {
			continue // the founders themselves
		}
		children++
		for g := 0; g < 2; g++ {
			for i := range loci {
				a, ok := pop.CreatedProvenance(id, g, i, len(loci))
				if !ok {
					t.Fatalf("child %d gamete %d locus %d: no created provenance", id, g, i)
				}
				transmitted[a] = true
				perParent[g][a] = true
			}
		}
	}
	t.Logf("%d children; the couple transmitted %d distinct created alleles "+
		"(father %d, mother %d)", children, len(transmitted), len(perParent[0]), len(perParent[1]))

	if len(transmitted) <= 4 {
		t.Errorf("the couple transmitted only %d distinct created alleles; germline "+
			"heterogeneity must exceed the four a strictly diploid pair can carry", len(transmitted))
	}
	for g, name := range map[int]string{0: "father", 1: "mother"} {
		if len(perParent[g]) <= 2 {
			t.Errorf("%s transmitted %d distinct alleles; a heterogeneous germline must "+
				"exceed the two somatic strands", name, len(perParent[g]))
		}
	}
	// Every transmitted allele must come from that parent's own pool: the partition is
	// disjoint, so a father transmitting one of the mother's alleles would be a bug.
	for g, parent := range []int{0, 1} {
		pool := map[int]bool{}
		for _, a := range pop.CreatedPool(parent) {
			pool[a] = true
		}
		for a := range perParent[g] {
			if !pool[a] {
				t.Errorf("parent %d transmitted allele %d, which is not in its pool", parent, a)
			}
		}
	}
}

// The genome must actually match the recorded provenance: at every focal locus the child's
// bit equals that bit of the created allele the record names. This is the W4 counterpart
// of TestARGRecordedStrandMatchesMeiosisResult — it checks the sequence the child really
// received, not just that the bookkeeping is self-consistent.
func TestCreatedGenomeMatchesRecordedAllele(t *testing.T) {
	const alleles = 8
	utils.SeedRNG(777)
	model := createdModel(alleles)
	pop := createdPop(model, alleles)
	runCreatedBirths(model, pop, 30)

	loci := core.EnsureARGLoci(model)
	checked := 0
	for id := range pop.IndData {
		if id < 2 {
			continue
		}
		chrom, ok := pop.Chromosomes[id]
		if !ok {
			continue
		}
		for g := 0; g < 2; g++ {
			for i, pos := range loci {
				a, ok := pop.CreatedProvenance(id, g, i, len(loci))
				if !ok {
					t.Fatalf("child %d gamete %d locus %d: no provenance", id, g, pos)
				}
				got := (chrom[g][pos/64] >> (uint(pos) % 64)) & 1
				want := (pop.CreatedAlleleSeqs[a][pos/64] >> (uint(pos) % 64)) & 1
				if got != want {
					t.Fatalf("child %d gamete %d locus %d: genome bit %d but created allele %d has %d",
						id, g, pos, got, a, want)
				}
				checked++
			}
		}
	}
	if checked == 0 {
		t.Fatal("no created gametes were checked; the test proved nothing")
	}
	t.Logf("verified %d (child, gamete, locus) triples against their created allele", checked)
}

// Recombination is untouched: a created gamete is still a two-allele mosaic, so adjacent
// focal loci usually agree and occasionally do not. A gamete carrying exactly one allele
// everywhere would mean the germ cell mechanism had replaced meiosis rather than fed it.
func TestCreatedGameteIsStillARecombinedMosaic(t *testing.T) {
	const alleles = 6
	utils.SeedRNG(31337)
	model := createdModel(alleles)
	model.StringParams["arg_loci"] = "0-255:1" // every bit, to see crossovers
	pop := createdPop(model, alleles)
	runCreatedBirths(model, pop, 30)

	loci := core.EnsureARGLoci(model)
	mosaics, single, switches, comparisons := 0, 0, 0, 0
	for id := range pop.IndData {
		if id < 2 {
			continue
		}
		for g := 0; g < 2; g++ {
			seen := map[int]bool{}
			prev, ok := pop.CreatedProvenance(id, g, 0, len(loci))
			if !ok {
				continue
			}
			seen[prev] = true
			for i := 1; i < len(loci); i++ {
				cur, _ := pop.CreatedProvenance(id, g, i, len(loci))
				seen[cur] = true
				if cur != prev {
					switches++
				}
				prev = cur
				comparisons++
			}
			if len(seen) > 2 {
				t.Errorf("child %d gamete %d drew from %d alleles; a germ cell carries two",
					id, g, len(seen))
			}
			if len(seen) == 2 {
				mosaics++
			} else {
				single++
			}
		}
	}
	if mosaics == 0 {
		t.Error("no gamete was a two-allele mosaic: recombination is not reaching the created path")
	}
	if switches == 0 {
		t.Error("provenance never switched allele along the chromosome")
	}
	rate := float64(switches) / float64(comparisons)
	if rate > 0.5 {
		t.Errorf("adjacent bits switch created allele %.3f of the time; linkage is broken", rate)
	}
	t.Logf("gametes: %d two-allele mosaics, %d single-allele; adjacent-bit switch rate %.4f",
		mosaics, single, rate)
}

// THE §2 GATE. With founder_allele_model unset, nothing about the run changes — no RNG is
// drawn for germ cells, no created records are written, and the population is identical.
func TestCreatedModelOffIsByteIdentical(t *testing.T) {
	run := func(created bool) ([]byte, map[int][]int, int) {
		utils.SeedRNG(4242)
		model := createdModel(6)
		if !created {
			delete(model.StringParams, "founder_allele_model")
		}
		pop := createdPop(model, 6)
		// createdPop consumes RNG identically in both arms (the allele sequences are
		// built either way), so any divergence below comes from the birth path.
		runCreatedBirths(model, pop, 25)
		return utils.RNGState(), pop.IndData, len(pop.CreatedARGDB)
	}
	offState, offInds, offRecs := run(false)
	onState, onInds, onRecs := run(true)

	if offRecs != 0 {
		t.Errorf("created provenance written with the model off: %d records", offRecs)
	}
	if onRecs == 0 {
		t.Error("no created provenance written with the model on")
	}
	if string(offState) == string(onState) {
		t.Error("the created model drew no extra RNG; germ-cell selection should consume two " +
			"draws per created gamete, so this test is not exercising the path")
	}
	if len(offInds) == 0 || len(onInds) == 0 {
		t.Fatal("no population")
	}
	// The off arm must match a run with no created state at all — the real byte-identity
	// claim, checked against a pop that never had CreatedAlleleSeqs.
	utils.SeedRNG(4242)
	model := createdModel(6)
	delete(model.StringParams, "founder_allele_model")
	pop := createdPop(model, 6)
	pop.CreatedAlleleSeqs = nil
	pop.FounderAllelePool = nil
	runCreatedBirths(model, pop, 25)
	if string(utils.RNGState()) != string(offState) {
		t.Error("the presence of created state perturbed a diploid-model run")
	}
}

// The walker must terminate at created ALLELES, not at founder strands: a lineage reaching
// created allele 7 has reached something with no ancestor and no coalescence time with
// allele 3. With A > 4 the sample can hold more than four created alleles at a locus, so
// TMR-4-A is censored while TMR-A-A is not — the exact reading TMR4A.md §0 says the metric
// forecloses.
func TestCreatedWalkerReportsMoreThanFourAlleles(t *testing.T) {
	const alleles = 10
	utils.SeedRNG(90210)
	model := createdModel(alleles)
	model.Parameters["tmrka_k"] = 4
	pop := createdPop(model, alleles)
	runCreatedBirths(model, pop, 40)

	ids := make([]int, 0, len(pop.IndData))
	for id := range pop.IndData {
		if id >= 2 { // sample the children, not the founders
			ids = append(ids, id)
		}
	}
	sort.Ints(ids)
	if len(ids) < 8 {
		t.Fatalf("only %d children; too few to sample", len(ids))
	}

	res := analysis.ComputeTMRKA(model, pop, ids)
	if res == nil || len(res.Loci) == 0 {
		t.Fatal("walk produced nothing")
	}
	aboveFour := 0
	for _, l := range res.Loci {
		if l.CreatedAllelesSurviving < 1 {
			t.Fatalf("locus %d: %d created alleles", l.Position, l.CreatedAllelesSurviving)
		}
		if l.CreatedAllelesSurviving > alleles {
			t.Errorf("locus %d: %d created alleles but only %d were created",
				l.Position, l.CreatedAllelesSurviving, alleles)
		}
		if l.UnlabelledSources != 0 {
			t.Errorf("locus %d: %d unlabelled sources under the created model",
				l.Position, l.UnlabelledSources)
		}
		if l.CreatedAllelesSurviving > 4 {
			aboveFour++
			// More than four created alleles means no TMR-4-A can exist: those lineages
			// never coalesce.
			if l.HasTMRKA {
				t.Errorf("locus %d reports a TMR-4-A of %d with %d created alleles surviving — "+
					"a censored walk was reported as a date",
					l.Position, l.TMRKAYear, l.CreatedAllelesSurviving)
			}
			if l.Outcome == analysis.OutcomeCoalesced {
				t.Errorf("locus %d: outcome %q with %d created alleles surviving",
					l.Position, l.Outcome, l.CreatedAllelesSurviving)
			}
		}
	}
	if aboveFour == 0 {
		t.Error("no locus retained more than four created alleles; one couple transmitting " +
			"more than four is the entire point of W4")
	}
	t.Logf("%d/%d loci retain more than four created alleles in a sample of %d children",
		aboveFour, len(res.Loci), len(ids))
}
