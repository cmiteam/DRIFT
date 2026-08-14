package simulation

import (
	"drift/pkg/core"
	"drift/pkg/individual"
	"drift/pkg/utils"
	"fmt"
	"math"
	"sort"
)

// birthStandard is the default birth algorithm, registered as "standard" in
// module_registry.go. Callers use the Birth dispatcher; add alternatives by
// registering another modules.BirthFunc under a different name.
func birthStandard(model *core.Model, pop *core.Pop) {

	// First, find eligible females and roll the dice
	var currentInds []int
	for ind := range pop.IndData {
		currentInds = append(currentInds, ind)
	}
	// Deterministic processing order (map iteration order is randomized), so a
	// fixed rng_seed reproduces the run.
	sort.Ints(currentInds)

	for _, ind := range currentInds {
		// skip males
		if pop.IndData[ind][individual.Sex] == 0 {
			continue
		}
		// skip unmarried women
		if pop.IndData[ind][individual.MarriageState] == -1 {
			continue
		}
		age := model.FreeParameters["year"] - pop.IndData[ind][individual.BirthYear]
		// skip adolescent girls
		if age < int(model.Parameters["maturity"]) {
			continue
		}
		// skip women in menopause
		if float64(age) > float64(pop.IndData[ind][individual.Lifespan])*model.Parameters["menopause"] {
			continue
		}
		// skip women with young children
		if int(pop.IndData[ind][individual.LastBirthYear])+int(model.Parameters["spacing"]) >= model.FreeParameters["year"] {
			continue
		}
		// Failed to get pregnant this year
		if utils.RandIntn(int(model.Parameters["birth_prob"])) != 0 {
			continue
		}

		// Next, put 'em in the oven
		mom := ind
		dad := pop.IndData[ind][individual.MarriageState]
		if _, dadExists := pop.IndData[dad]; !dadExists {
			fmt.Printf("WARNING: Woman %d claims to be married to non-existent man %d\n", ind, dad)
			// Clear the invalid marriage state
			pop.IndData[ind][individual.MarriageState] = -1
			continue
		}
		// Fitness selection is governed by selection_mode (roadmap §2): it can act on
		// fecundity (here), viability (death.go), both, or neither.

		// the average of the maternal and paternal fitness affects birth probability
		fitness := 1.0
		if fecunditySelection(model) {
			pfit := float64(pop.IndData[dad][individual.Fitness])
			mfit := float64(pop.IndData[mom][individual.Fitness])
			fitness = (pfit + mfit) / 2
			fitness = fitness / model.Parameters["mu_scale_factor"]
		}
		chance := utils.RandFloat64()
		if chance < fitness {
			model.FreeParameters["indID"] += 1
			child := model.FreeParameters["indID"]
			createChild(model, pop, dad, mom, child)
			pop.Tracking["births"]++
			// If DNA or mutations are being tracked, bitmasks will be created that
			// will be used to control meiosis and mutation inheritance. These will
			// be used for both meiosis and mutation inheritance, so we will set
			// them up once and use them at will.
			// track_arg joins the mask-building condition because strand provenance
			// (TMR4A.md W2) is READ OFF these masks. On a run that already tracks DNA
			// or mutations that costs nothing and changes no RNG draw; on a run that
			// tracks neither it causes masks to be drawn that otherwise would not be,
			// which changes the stream — opt-in, and documented in pkg/core/arg.go.
			argOn := core.ARGFocalOn(model)
			if model.Parameters["track_DNA"] > 0 || model.Parameters["track_mutations"] > 0 || argOn {
				var genomemask0, genomemask1 []uint64
				var centsmask0, centsmask1 []uint64
				genomemask0, centsmask0 = createMask(model, 0)
				genomemask1, centsmask1 = createMask(model, 1)

				// Capture which parental strand each focal locus came from, before the
				// masks are consumed and discarded. Requires track_pedigree: without the
				// parent links these records are unreadable, and nothing would ever
				// prune them.
				// CREATED-FOUNDER GERMLINE (TMR4A.md W4), hard-gated and implemented in
				// created_meiosis.go. When the created model is active and this parent is
				// a founder with a created-allele pool, the germ cell producing this
				// gamete carries two alleles drawn from that pool instead of the parent's
				// two somatic strands — which is how one couple transmits more than four
				// alleles at a locus. Everything downstream (mask, recombination,
				// polarity) is unchanged. Two RNG draws per created gamete, in a fixed
				// paternal-then-maternal order; nothing at all under the default model.
				var germ [2][2]int
				var germOK [2]bool
				germ[0], germOK[0] = createdGermCell(model, pop, dad)
				germ[1], germOK[1] = createdGermCell(model, pop, mom)

				if argOn && model.Parameters["track_pedigree"] == 1 {
					loci := core.EnsureARGLoci(model)
					pop.ARGCapture(child, loci, genomemask0, genomemask1)
					// A gamete drawn from a created pool also records WHICH created
					// allele it took, in a separate record the walker prefers over the
					// strand bit above.
					pop.CreatedCapture(child, loci,
						[2][]uint64{genomemask0, genomemask1}, germ, germOK)
				}

				// Add tracked DNA
				if model.Parameters["track_DNA"] > 0 {
					// only create a child's chromosomes if there is something to track in at least one parent
					if pop.IndData[dad][individual.AlleleCount] > 0 || pop.IndData[mom][individual.AlleleCount] > 0 || germOK[0] || germOK[1] {
						pop.Chromosomes[child] = [][]uint64{make([]uint64, (model.FreeParameters["genome_bits"]+63)/64), make([]uint64, (model.FreeParameters["genome_bits"]+63)/64)}
					}
					numSetBits := 0
					// only go through meiosis if there is a set bit in mom or dad — or,
					// under the created model, if this parent presented a germ cell, whose
					// alleles are the source instead of the parent's somatic strands.
					if germOK[0] {
						createdMeiosis(pop, genomemask0, germ[0], child, 0)
						numSetBits += utils.CountSetBits(pop.Chromosomes[child][0])
					} else if pop.IndData[dad][individual.AlleleCount] > 0 {
						meiosis(pop, genomemask0, dad, child, 0)
						numSetBits += utils.CountSetBits(pop.Chromosomes[child][0])
					}
					if germOK[1] {
						createdMeiosis(pop, genomemask1, germ[1], child, 1)
						numSetBits += utils.CountSetBits(pop.Chromosomes[child][1])
					} else if pop.IndData[mom][individual.AlleleCount] > 0 {
						meiosis(pop, genomemask1, mom, child, 1)
						numSetBits += utils.CountSetBits(pop.Chromosomes[child][1])
					}
					pop.IndData[child][individual.AlleleCount] = numSetBits
					// delete the child's chromosomes if they inherited zero set bits
					if pop.IndData[child][individual.AlleleCount] < 1 {
						delete(pop.Chromosomes, child)
					} else {
						pop.IndData[child][individual.NumBlocks] = utils.CountContiguousBlocks(model, pop, child, 0)
						pop.IndData[child][individual.NumBlocks] += utils.CountContiguousBlocks(model, pop, child, 1)
					}

					// inherit centromeres if mom or dad have a set bit in their centromeres
					if pop.IndData[dad][individual.NumCentromeres] > 0 || pop.IndData[mom][individual.NumCentromeres] > 0 {
						inheritCentromeres(model, pop, centsmask0, centsmask1, dad, mom, child)
					}

					// track avenues of descent from the seed individual(s)
					pop.IndData[child][individual.YGens] = -1
					if pop.IndData[dad][individual.YGens] > -1 && pop.IndData[child][individual.Sex] == 0 {
						pop.IndData[child][individual.YGens] = pop.IndData[dad][individual.YGens] + 1
					}
					pop.IndData[child][individual.MtGens] = -1
					if pop.IndData[mom][individual.MtGens] > -1 {
						pop.IndData[child][individual.MtGens] = pop.IndData[mom][individual.MtGens] + 1
					}

					pop.IndData[child][individual.MinGenealoGens] = -1
					minGenealo := pop.IndData[dad][individual.MinGenealoGens]
					if pop.IndData[mom][individual.MinGenealoGens] > minGenealo {
						minGenealo = pop.IndData[mom][individual.MinGenealoGens]
					}
					if minGenealo > -1 {
						pop.IndData[child][individual.MinGenealoGens] = minGenealo + 1
					}

					pop.IndData[child][individual.MaxGenealoGens] = -1
					maxGenealo := pop.IndData[dad][individual.MaxGenealoGens]
					if pop.IndData[mom][individual.MaxGenealoGens] > maxGenealo {
						maxGenealo = pop.IndData[mom][individual.MaxGenealoGens]
					}
					if maxGenealo > -1 {
						pop.IndData[child][individual.MaxGenealoGens] = maxGenealo + 1
					}
				}

				// Assign mutations, both inherited and de novo
				if model.Parameters["track_mutations"] > 0 {
					utils.InheritMutations(model, pop, genomemask0, dad, child, 0)
					utils.InheritMutations(model, pop, genomemask1, mom, child, 1)
					utils.GenerateNewMutations(model, pop, child)
					numMutations, mutationFitness := utils.CountFitnessAndMutations(pop, child, model)
					fitness := mutationFitness
					pop.IndData[child][individual.Fitness] = int(float64(fitness) * model.Parameters["mu_scale_factor"])
					pop.IndData[child][individual.NumMutations] = numMutations
				}
			}

			// Diploid pedigree retention (TMR4A.md W1). Unconditional on child sex —
			// that asymmetry is exactly what makes MaleDB/FemaleDB below a pair of
			// single-sex chains rather than a pedigree, and the TMR-K-A walker needs
			// both parents of every child. Parents are ensured first so a founder gets
			// a node the moment it reproduces (and so the child's increment lands on
			// something); IndData is still live for both, since a parent cannot be
			// dead. No RNG, so a track_pedigree run stays byte-identical.
			if model.Parameters["track_pedigree"] == 1 {
				pop.PedEnsure(dad, pop.IndData[dad][individual.BirthYear], pop.IndData[dad][individual.Sex])
				pop.PedEnsure(mom, pop.IndData[mom][individual.BirthYear], pop.IndData[mom][individual.Sex])
				pop.PedRecordBirth(child, dad, mom,
					pop.IndData[child][individual.BirthYear],
					pop.IndData[child][individual.Sex])
			}

			if model.Parameters["track_coalescence"] == 1 {
				if pop.IndData[child][individual.Sex] == 0 {
					pop.IndData[dad][individual.Sons]++
					pop.MaleDB[child] = core.Ancestor{
						ID:        dad,
						BirthYear: pop.IndData[dad][individual.BirthYear],
					}
				} else {
					pop.IndData[mom][individual.Daughters]++
					pop.FemaleDB[child] = core.Ancestor{
						ID:        mom,
						BirthYear: pop.IndData[mom][individual.BirthYear],
					}
				}
			}
		}
	}
}

func createChild(model *core.Model, pop *core.Pop, dad, mom, child int) {

	// potential lifespan is the average of the parents X the lifespan drop per generation, but it bottoms out at min_lifespan
	lifespan := int(float64(pop.IndData[dad][individual.Lifespan]+pop.IndData[mom][individual.Lifespan]) / 2.0 * model.Parameters["lifespan_drop"])
	if lifespan < int(model.Parameters["min_lifespan"]) {
		lifespan = int(model.Parameters["min_lifespan"])
	}
	childLat, childLon := utils.Wander(model, pop.IndData[dad][individual.Lat], pop.IndData[dad][individual.Lon])

	pop.IndData[child] = individual.MakeIndData()
	pop.IndData[child][individual.Dad] = dad
	pop.IndData[child][individual.Mom] = mom
	// Deme is inherited from one parent (pairs with individual-level migration in
	// Mating). Paternal by default; set deme_inheritance="maternal" to switch. In a
	// patrilocal island model the father's deme is the natural default, and it makes
	// deme membership consistent with the Y-chromosome line. With num_demes=1 every
	// deme is 0 so the choice has no effect (the single-population no-op is preserved).
	demeParent := dad
	if model.StringParam("deme_inheritance", "paternal") == "maternal" {
		demeParent = mom
	}
	pop.IndData[child][individual.Deme] = pop.IndData[demeParent][individual.Deme]
	pop.IndData[child][individual.Sex] = utils.RandIntn(2)
	pop.IndData[child][individual.BirthYear] = model.FreeParameters["year"]
	pop.IndData[child][individual.Lifespan] = lifespan
	pop.IndData[child][individual.Lat] = childLat
	pop.IndData[child][individual.Lon] = childLon
	pop.IndData[child][individual.NumCentromeres] = 0 //TODO guessing

	// Track lines of descent from the seed individual(s)
	pop.IndData[child][individual.YGens] = -1
	if pop.IndData[dad][individual.YGens] >= 0 && pop.IndData[child][individual.Sex] == 0 {
		pop.IndData[child][individual.YGens] = pop.IndData[dad][individual.YGens] + 1
	}

	pop.IndData[child][individual.MtGens] = -1
	if pop.IndData[mom][individual.MtGens] >= 0 {
		pop.IndData[child][individual.MtGens] = pop.IndData[mom][individual.MtGens] + 1
	}

	pop.IndData[child][individual.MinGenealoGens] = -1
	minGenealo := pop.IndData[dad][individual.MinGenealoGens]
	if pop.IndData[mom][individual.MinGenealoGens] < minGenealo && pop.IndData[mom][individual.MinGenealoGens] >= 0 {
		minGenealo = pop.IndData[mom][individual.MinGenealoGens]
	} else if minGenealo < 0 && pop.IndData[mom][individual.MinGenealoGens] >= 0 {
		minGenealo = pop.IndData[mom][individual.MinGenealoGens]
	}
	if minGenealo >= 0 {
		pop.IndData[child][individual.MinGenealoGens] = minGenealo + 1
	}

	pop.IndData[child][individual.MaxGenealoGens] = -1
	maxGenealo := pop.IndData[dad][individual.MaxGenealoGens]
	if pop.IndData[mom][individual.MaxGenealoGens] > maxGenealo {
		maxGenealo = pop.IndData[mom][individual.MaxGenealoGens]
	}
	if maxGenealo >= 0 {
		pop.IndData[child][individual.MaxGenealoGens] = maxGenealo + 1
	}

	pop.IndData[mom][individual.LastBirthYear] = model.FreeParameters["year"]
	if numBirths := pop.IndData[mom][individual.NumBirths]; numBirths == individual.EmptyField {
		pop.IndData[mom][individual.NumBirths] = 0
	}
	pop.IndData[mom][individual.NumBirths]++
}

func createMask(model *core.Model, sex int) ([]uint64, []uint64) {
	// SAFETY CHECK - prevent infinite loop
	genomeBits := model.FreeParameters["genome_bits"]
	if genomeBits <= 0 {
		fmt.Printf("ERROR: genome_bits is %d, cannot create mask\n", genomeBits)
		return []uint64{}, []uint64{0}
	}

	// Real-recombination-map meiosis (roadmap §6a), opt-in via recombination_model.
	// The default "legacy" path below is the original single-interior-segment model
	// and is left byte-for-byte untouched (same 3 RNG draws per chromosome in sorted
	// order), so strictly-neutral runs stay byte-identical. Only "map" mode places
	// crossovers per a genetic map; it consumes a DIFFERENT (documented, deterministic)
	// RNG stream, which is fine because it is opt-in. Falls through to legacy if no map
	// can be built (e.g. no arms). Centromere handling is unchanged in both modes.
	if model.StringParam("recombination_model", "legacy") == "map" {
		if gm := utils.EnsureGeneticMap(model); gm != nil {
			return createMaskFromMap(model, gm), []uint64{0}
		}
	}

	genomeArrSize := (model.FreeParameters["genome_bits"] + 63) / 64
	genomemask := make([]uint64, genomeArrSize)
	centromask := []uint64{0}

	// Iterate chromosomes in a fixed order: this loop draws RNG per chromosome,
	// and Go map iteration order is randomized, so unsorted iteration would make
	// meiosis non-reproducible even under a fixed rng_seed.
	chroms := make([]int, 0, len(model.ChromosomeArms))
	for chrom := range model.ChromosomeArms {
		chroms = append(chroms, chrom)
	}
	sort.Ints(chroms)

	//    fmt.Printf("\n=== DEBUG createMask ===\n")
	//    fmt.Printf("genome_bits: %d\n", model.FreeParameters["genome_bits"])
	//    fmt.Printf("genomeArrSize: %d\n", genomeArrSize)
	//    fmt.Printf("Number of chromosomes: %d\n", len(model.ChromosomeArms))

	for _, chrom := range chroms {
		//        fmt.Printf("\nChromosome %d:\n", chrom)

		// Check if arms exist
		arm0, ok0 := model.ChromosomeArms[chrom][0]
		arm1, ok1 := model.ChromosomeArms[chrom][1]

		if !ok0 || !ok1 {
			fmt.Printf("  SKIPPING - missing arm\n")
			continue
		}

		pstart := arm0[0]
		plen := arm0[1]
		qstart := arm1[0]
		qlen := arm1[1]

		if plen <= 0 || qlen <= 0 {
			fmt.Printf("  SKIPPING - invalid length\n")
			continue
		}

		ploc := utils.RandIntn(plen)
		qloc := utils.RandIntn(qlen)
		whichCopy := utils.RandIntn(2)

		if whichCopy == 1 {
			startBit := pstart + ploc
			endBit := qstart + qloc

			for i := startBit; i < endBit; i++ {
				genomemask[i/64] |= (1 << (i % 64))
			}
			// centromere stuff...
		}
	}

	return genomemask, centromask
}

// createMaskFromMap builds a meiosis mask by placing crossovers along a real genetic
// map (roadmap §6a), replacing the legacy single-interior-segment model. Per
// chromosome (in sorted key order, so the run stays reproducible under a fixed seed):
//
//   - crossover COUNT = 1 obligate crossover (crossover assurance — every bivalent
//     recombines, chosen with Rob) + Poisson(max(0, totalCM - recomb_obligate_cM)/100)
//     extra crossovers, so the count scales with the chromosome's genetic length while
//     small chromosomes still get their one guaranteed crossover;
//   - crossover POSITIONS are drawn uniformly in genetic (cM) space and inverted to
//     bit positions via the map's CDF, so they concentrate where cM/bit is high;
//   - the STARTING homolog is chosen 50/50 (independent assortment — fixes the legacy
//     model's defect where a chromosome's telomeres always came from copy 1).
//
// The mask alternates parental copy at each crossover: where the current copy is 0 the
// bit is set (meiosis then takes copy 0 there), matching the legacy convention
// (mask bit 1 => copy 0). RNG order per chromosome: Poisson(count) → count× position
// draws → start-copy — a fixed order documented here so the stream is deterministic.
func createMaskFromMap(model *core.Model, gm *core.GeneticMap) []uint64 {
	genomeArrSize := (model.FreeParameters["genome_bits"] + 63) / 64
	mask := make([]uint64, genomeArrSize)

	baseline := 50.0 // obligate-crossover cM reservation (extras are Poisson beyond it)
	if v, ok := model.Parameters["recomb_obligate_cM"]; ok {
		baseline = v
	}

	chroms := make([]int, 0, len(model.ChromosomeArms))
	for chrom := range model.ChromosomeArms {
		chroms = append(chroms, chrom)
	}
	sort.Ints(chroms)

	for _, chrom := range chroms {
		if !gm.Has(chrom) {
			continue
		}
		totalCM := gm.TotalCM(chrom)
		start, end := gm.Span(chrom)

		var xovers []int
		if totalCM > 0 {
			nExtra := utils.RandPoisson(math.Max(0, totalCM-baseline) / 100.0)
			nX := 1 + nExtra // obligate crossover + Poisson extras
			for k := 0; k < nX; k++ {
				xovers = append(xovers, gm.BitAtCM(chrom, utils.RandFloat64()*totalCM))
			}
			sort.Ints(xovers)
		}
		// A chromosome with totalCM == 0 draws no crossovers but still assorts.
		cur := utils.RandIntn(2) // starting homolog (independent assortment)

		prev := start
		for _, xb := range xovers {
			if cur == 0 {
				setMaskBits(mask, prev, xb)
			}
			cur ^= 1
			prev = xb
		}
		if cur == 0 {
			setMaskBits(mask, prev, end)
		}
	}
	return mask
}

// setMaskBits sets bits [lo, hi) in the genome word array.
func setMaskBits(mask []uint64, lo, hi int) {
	for i := lo; i < hi; i++ {
		mask[i/64] |= 1 << (uint(i) % 64)
	}
}

// Meiosis simulates genetic recombination during gamete formation.
// By applying a combination of AND, OR, and NOT between the genome mask and the
// two parental chromosome copies, the child can get a haploid, recombined
// version of a parent's genome. This happens once for the father and once for
// the mother, so chromosomes[child][0] = paternal inheritance and
// chromosomes[child][1] = maternal inheritance

func meiosis(pop *core.Pop, mask []uint64, parent int, child int, copy int) {
	pop.Chromosomes[child][copy] = meiosisFrom(mask, pop.Chromosomes[parent][0], pop.Chromosomes[parent][1])
}

// meiosisFrom is the meiosis kernel over two explicit source strands: where the mask bit
// is SET the child takes src0, where it is CLEAR it takes src1. Factored out so the
// created-founder germline path (TMR4A.md W4, created_meiosis.go) applies the identical
// polarity to a germ cell's two created alleles instead of a parent's two somatic strands.
// One definition, so the two paths cannot drift apart — and so the strand-provenance
// convention recorded in pkg/core/arg.go has exactly one thing to agree with.
func meiosisFrom(mask, src0, src1 []uint64) []uint64 {
	childCopy := make([]uint64, len(mask))
	for i := 0; i < len(mask); i++ {
		childCopy[i] = (mask[i] & src0[i]) | (^mask[i] & src1[i])
	}
	return childCopy
}

func inheritCentromeres(model *core.Model, pop *core.Pop, centsmask0 []uint64, centsmask1 []uint64, dad int, mom int, child int) {
	if pop.IndData[dad][individual.NumCentromeres] > 0 || pop.IndData[mom][individual.NumCentromeres] > 0 {
		_, dadExists := pop.Centromeres[dad]
		_, momExists := pop.Centromeres[mom]

		var blankCentromeres [2][]uint64
		blankCentromeres[0] = make([]uint64, len(centsmask0))
		blankCentromeres[1] = make([]uint64, len(centsmask1))
		pop.Centromeres[child] = blankCentromeres

		for chrom, _ := range model.ChromosomeArms {
			idx := chrom / 64
			bit := chrom % 64
			if dadExists {
				if centsmask0[idx]&(1<<bit) == 0 {
					pop.Centromeres[child][0][idx] |= (pop.Centromeres[dad][0][idx] & (1 << bit))
				} else {
					pop.Centromeres[child][0][idx] |= (pop.Centromeres[dad][1][idx] & (1 << bit))
				}
			}
			if momExists {
				if centsmask1[idx]&(1<<bit) == 0 {
					pop.Centromeres[child][1][idx] |= (pop.Centromeres[mom][0][idx] & (1 << bit))
				} else {
					pop.Centromeres[child][1][idx] |= (pop.Centromeres[mom][1][idx] & (1 << bit))
				}
			}
		}

		centromereCount := utils.CountSetBits(pop.Centromeres[child][0])
		centromereCount += utils.CountSetBits(pop.Centromeres[child][1])
		pop.IndData[child][individual.NumCentromeres] = centromereCount
		if centromereCount == 0 {
			delete(pop.Centromeres, child) // child has been a disappointment
		}
	}
}
