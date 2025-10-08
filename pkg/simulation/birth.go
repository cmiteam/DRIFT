package simulation

import (
	"drift/modules/individual"
	"drift/modules/mutation"
	"drift/pkg/core"
	"drift/pkg/utils"
	"fmt"
	"math/rand"
)

func Birth(model *core.Model, pop *core.Pop) {

	// First, find eligible females and roll the dice
	var currentInds []int
	for ind := range pop.IndData {
		currentInds = append(currentInds, ind)
	}

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
		if rand.Intn(int(model.Parameters["birth_prob"])) != 0 {
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
		// TO DO: fitness ALSO affects survivorship each year, work out a way to
		// use fitness for birth OR survivorship OR both

		// the average of the maternal and paternal fitness affects birth probability
		fitness := 1.0
		if model.Parameters["track_mutations"] == 1 {
			pfit := float64(pop.IndData[dad][individual.Fitness])
			mfit := float64(pop.IndData[mom][individual.Fitness])
			fitness = (pfit + mfit) / 2
			fitness = fitness / model.Parameters["mu_scale_factor"]
		}
		chance := rand.Float64()
		if chance < fitness {
			model.FreeParameters["indID"] += 1
			child := model.FreeParameters["indID"]
			createChild(model, pop, dad, mom, child)
			pop.Tracking["births"]++
			// If DNA or mutations are being tracked, bitmasks will be created that
			// will be used to control meiosis and mutation inheritance. These will
			// be used for both meiosis and mutation inheritance, so we will set
			// them up once and use them at will.
			if model.Parameters["track_DNA"] > 0 || model.Parameters["track_mutations"] > 0 {
				var genomemask0, genomemask1 []uint64
				var centsmask0, centsmask1 []uint64
				genomemask0, centsmask0 = createMask(model, 0)
				genomemask1, centsmask1 = createMask(model, 1)

				// Add tracked DNA
				if model.Parameters["track_DNA"] > 0 {
					// only create a child's chromosomes if there is something to track in at least one parent
					if pop.IndData[dad][individual.AlleleCount] > 0 || pop.IndData[mom][individual.AlleleCount] > 0 {
						pop.Chromosomes[child] = [][]uint64{make([]uint64, (model.FreeParameters["genome_bits"]+63)/64), make([]uint64, (model.FreeParameters["genome_bits"]+63)/64)}
					}
					numSetBits := 0
					// only go through meiosis if there is a set bit in mom or dad
					if pop.IndData[dad][individual.AlleleCount] > 0 {
						meiosis(pop, genomemask0, dad, child, 0)
						numSetBits += utils.CountSetBits(pop.Chromosomes[child][0])
					}
					if pop.IndData[mom][individual.AlleleCount] > 0 {
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
					mutation.InheritMutations(pop, genomemask0, dad, child, 0)
					mutation.InheritMutations(pop, genomemask1, mom, child, 1)
					mutation.GenerateNewMutations(model, pop, child)
					numMutations, mutationLoad := mutation.CountFitnessAndMutations(pop, child)
					fitness := 1 + mutationLoad
					pop.IndData[child][individual.Fitness] = int(float64(fitness) * model.Parameters["mu_scale_factor"])
					pop.IndData[child][individual.NumMutations] = numMutations
				}
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
	pop.IndData[child][individual.Sex] = rand.Intn(2)
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

	// masks are uint64 (8-byte unsigned integers with 64 bits of memory). It takes about 50 uint64 to code for one copy of a 3,100 bit genome
	// the centromere mask is a single uint64, therefore models with up to 64 chromosomes can be handled

	genomeArrSize := (model.FreeParameters["genome_bits"] + 63) / 64
	genomemask := make([]uint64, genomeArrSize)
	centromask := []uint64{0}

	for chrom, _ := range model.ChromosomeArms {
		// in biology, chromosomes generally have a shorter 'p' arm and a longer 'q' arm', the lengths were loaded previously
		// chromosomeArms[chrom][0] = p, chromosomeArms[chrom][1] = q
		// chromosomeArms[chrom][0][0] = start of p arm in bits, chromosomeArms[chrom][0][1] = length of p arm in bits
		pstart := model.ChromosomeArms[chrom][0][0]
		qstart := model.ChromosomeArms[chrom][1][0]
		plen := model.ChromosomeArms[chrom][0][1]
		qlen := model.ChromosomeArms[chrom][1][1]

		if plen <= 0 || qlen <= 0 {
			continue
		}

		// choose a random place on each chromosome arm and decide if the paternal
		// or maternal centromere will be inherited by the child
		ploc := rand.Intn(plen)
		qloc := rand.Intn(qlen)
		whichCopy := rand.Intn(2)

		if whichCopy == 1 {
			// example: 00001111x11110000, where 0 = paternal, 1 = maternal, and x = the centromere
			for i := pstart + ploc; i < qstart+qloc; i++ {
				genomemask[i/64] |= (1 << (i % 64))
			}
			// make sure we have enough space in the centromask
			for len(centromask) < chrom/64 {
				centromask = append(centromask, uint64(0))
			}
			centromask[chrom/64] |= (1 << (chrom % 64))
		} else {
			// example: 11110000x00001111
			for i := pstart; i < pstart+ploc; i++ {
				genomemask[i/64] |= (1 << (i % 64))
			}
			for i := qstart + qloc; i < qstart+qlen; i++ {
				genomemask[i/64] |= (1 << (i % 64))
			}
		}
	}

	if sex == 0 && model.Parameters["sex_chrom_index"] > 0 {
		// males don't inherit the father's X
		sexIdx := int(model.Parameters["sex_chrom_index"])
		xstart := model.ChromosomeArms[sexIdx][0][0]
		xend := model.ChromosomeArms[1][1][0] + model.ChromosomeArms[1][1][1]
		for i := xstart; i < xend; i++ {
			genomemask[i/64] &^= (1 << (i % 64))
		}
	}
	return genomemask, centromask
}

// Meiosis simulates genetic recombination during gamete formation.
// By applying a combination of AND, OR, and NOT between the genome mask and the
// two parental chromosome copies, the child can get a haploid, recombined
// version of a parent's genome. This happens once for the father and once for
// the mother, so chromosomes[child][0] = paternal inheritance and
// chromosomes[child][1] = maternal inheritance

func meiosis(pop *core.Pop, mask []uint64, parent int, child int, copy int) {
	parentCopy0 := pop.Chromosomes[parent][0]
	parentCopy1 := pop.Chromosomes[parent][1]
	childCopy := make([]uint64, len(mask))
	for i := 0; i < len(mask); i++ {
		childCopy[i] = (mask[i] & parentCopy0[i]) | (^mask[i] & parentCopy1[i])
	}
	pop.Chromosomes[child][copy] = childCopy
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

		centromereCount := countSetBits(pop.Centromeres[child][0])
		centromereCount += countSetBits(pop.Centromeres[child][1])
		pop.IndData[child][individual.NumCentromeres] = centromereCount
		if centromereCount == 0 {
			delete(pop.Centromeres, child) // child has been a disappointment
		}
	}
}
