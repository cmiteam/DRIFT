package birth

import (
	"drift/modules/mutation"
	"drift/types"
	"fmt"
	"math/rand"
	"strings"
)

func Birth(model *types.Model, pop *types.Pop, year int) {

	// First, find eligible females and roll the dice
	for ind := range pop.IndData {
		// skip males
		if pop.IndData[ind]["sex"] == 0 {
			continue
		}
		// skip unmarried women
		if pop.IndData[ind]["marriage_state"] == -1 {
			continue
		}
		age := year - pop.IndData[ind]["birth_year"]
		// skip adolescent girls
		if age < int(model.Parameters["maturity"]) {
			continue
		}
		// skip women in menopause
		if float64(age) > float64(pop.IndData[ind]["lifespan"])*model.Parameters["menopause"] {
			continue
		}
		// skip women with young children
		if int(pop.IndData[ind]["last_birth_year"])+int(model.Parameters["spacing"]) >= year {
			continue
		}
		// Failed to get pregnant this year
		if rand.Intn(int(model.Parameters["birth_prob"])) != 0 {
			continue
		}

		// Next, put 'em in the oven
		mom := ind
		dad := pop.IndData[ind]["marriage_state"]
		// TO DO: fitness ALSO affects survivorship each year, work out a way to
		// use fitness for birth OR survivorship OR both

		// the average of the maternal and paternal fitness affects birth probability
		fitness := 1.0
		if model.Parameters["track_mutations"] == 1 {
			pfit := float64(pop.IndData[dad]["fitness"])
			mfit := float64(pop.IndData[mom]["fitness"])
			fitness = (pfit + mfit) / 2
			fitness = fitness / model.Parameters["mu_Scale_factor"]
		}
		chance := rand.Float64()
		if chance < fitness {
			model.FreeParameters["indID"] += 1
			child := model.FreeParameters["indID"]
			createChild(model, pop, dad, mom, child, year)

			// If DNA or mutations are being tracked, bitmasks will be created that
			// will be used to control meiosis and mutation inheritance. These will
			// be used for both meiosis and mutation inheritance, so we will set
			// them up once and use them at will.

			if model.Parameters["track_DNA"] > 0 || model.Parameters["track_mutations"] > 0 {
				var genomemask1, genomemask2 []uint64
				var centsmask1, centsmask2 uint64
				genomemask1, centsmask1 = createMask(model, 0)
				genomemask2, centsmask2 = createMask(model, 1)

				// Add tracked DNA
				if model.Parameters["track_DNA"] > 0 {
					// only create a child's chromosomes if there is something to track at least one parent
					if pop.IndData[dad]["allele_count"] > 0 || pop.IndData[mom]["allele_count"] > 0 {
						pop.Chromosomes[child] = [][]uint64{make([]uint64, (model.FreeParameters["NumBits"]+63)/64), make([]uint64, (model.FreeParameters["NumBits"]+63)/64)}
					}
					numSetBits := 0
					// only go through meiosis if there is a set bit in mom or dad
					if pop.IndData[dad]["allele_count"] > 0 {
						meiosis(pop, genomemask1, dad, child, 0)
						numSetBits += countSetBits(pop.Chromosomes[child][0])
					}
					if pop.IndData[mom]["allele_count"] > 0 {
						meiosis(pop, genomemask2, mom, child, 1)
						numSetBits += countSetBits(pop.Chromosomes[child][1])
					}
					pop.IndData[child]["allele_count"] = numSetBits
					// delete the child's chromosomes if they inherited zero set bits
					if pop.IndData[child]["allele_count"] < 1 {
						delete(pop.Chromosomes, child)
					} else {
						pop.IndData[child]["num_blocks"] = countContiguousBlocks(model, pop, child, 0)
						pop.IndData[child]["num_blocks"] += countContiguousBlocks(model, pop, child, 1)
					}

					// inherit centromeres if mom or dad have a set bit in their centromeres
					if pop.IndData[dad]["centromere_count"] > 0 || pop.IndData[mom]["centromere_count"] > 0 {
						inheritCentromeres(model, pop, centsmask1, centsmask2, dad, mom, child)
					}

					// track avenues of descent from the seed individual(s)
					pop.IndData[child]["Y_gens"] = -1
					if pop.IndData[dad]["Y_gens"] > -1 && pop.IndData[child]["sex"] == 0 {
						pop.IndData[child]["Y_gens"] = pop.IndData[dad]["Y_gens"] + 1
					}
					pop.IndData[child]["mt_gens"] = -1
					if pop.IndData[mom]["mt_gens"] > -1 {
						pop.IndData[child]["mt_gens"] = pop.IndData[mom]["mt_gens"] + 1
					}

					pop.IndData[child]["min_genealo_gens"] = -1
					minGenealo := pop.IndData[dad]["min_genealo_gens"]
					if pop.IndData[mom]["min_genealo_gens"] > minGenealo {
						minGenealo = pop.IndData[mom]["min_genealo_gens"]
					}
					if minGenealo > -1 {
						pop.IndData[child]["min_genealo_gens"] = minGenealo + 1
					}

					pop.IndData[child]["max_genealo_gens"] = -1
					maxGenealo := pop.IndData[dad]["max_genealo_gens"]
					if pop.IndData[mom]["max_genealo_gens"] > maxGenealo {
						maxGenealo = pop.IndData[mom]["max_genealo_gens"]
					}
					if maxGenealo > -1 {
						pop.IndData[child]["max_genealo_gens"] = maxGenealo + 1
					}
				}

				// Assign mutations, both inherited and de novo
				if model.Parameters["track_mutations"] > 0 {
					mutation.InheritMutations(pop, genomemask1, dad, child, 0)
					mutation.InheritMutations(pop, genomemask2, mom, child, 1)
					mutation.GenerateNewMutations(model, pop, child)
					numMutations, mutationLoad := mutation.CountFitnessAndMutations(pop, child)
					fitness := 1 + mutationLoad
					pop.IndData[child]["fitness"] = int(float64(fitness) * model.Parameters["mu_scale_factor"])
					pop.IndData[child]["num_mutations"] = numMutations
				}
			}
			pop.Tracking["births"]++
		}
	}
}

func createChild(model *types.Model, pop *types.Pop, dad, mom, child int, year int) {

	// potential lifespan is the average of the parents X the lifespan drop per generation, but it bottoms out at min_lifespan
	lifespan := int((pop.IndData[dad]["lifespan"] + pop.IndData[mom]["lifespan"]) / 2 * int(model.Parameters["lifespan_drop"]))
	if lifespan < int(model.Parameters["min_lifespan"]) {
		lifespan = int(model.Parameters["min_lifespan"])
	}

	pop.IndData[child] = map[string]int{
		"dad":            dad,
		"mom":            mom,
		"sex":            rand.Intn(2),
		"birth_year":     year,
		"lifespan":       lifespan,
		"marriage_state": -1,
		"lat":            0,
		"lon":            0,
	}

	pop.IndData[mom]["last_birth_year"] = year
	if _, exists := pop.IndData[mom]["numbirths"]; !exists {
		pop.IndData[mom]["numbirths"] = 0
	}
	pop.IndData[mom]["numbirths"]++
}

func createMask(model *types.Model, sex int) ([]uint64, uint64) {

	// masks are uint64 (8-byte unsigned integers with 64 bits of memory). It takes about 50 uint64 to code for one copy of a 3,100 bit genome
	// the centromere mask is a single uint64, therefore models with up to 64 chromosomes can be handled

	numUint64s := (model.FreeParameters["numbits"] + 63) / 64
	genomemask := make([]uint64, numUint64s)
	var centromask uint64

	for chrom := 1; chrom < len(model.ChromosomeArms); chrom++ {
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
			centromask |= (1 << (chrom % 64))
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

	if sex == 0 {
		// males don't inherit the father's X chromosome
		xstart := model.ChromosomeArms[23][0][0]
		xend := model.ChromosomeArms[23][1][0] + model.ChromosomeArms[23][1][1]
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

func meiosis(pop *types.Pop, mask []uint64, parent int, child int, copy int) {
	parentCopy0 := pop.Chromosomes[parent][0]
	parentCopy1 := pop.Chromosomes[parent][1]
	childCopy := make([]uint64, len(mask))
	for i := 0; i < len(mask); i++ {
		childCopy[i] = (mask[i] & parentCopy0[i]) | (^mask[i] & parentCopy1[i])
	}
	pop.Chromosomes[child][copy] = childCopy
}

// countContiguousBlocks counts blocks of contiguous set bits
func countContiguousBlocks(model *types.Model, pop *types.Pop, ind int, copy int) int {
	blockCount := 0
	genomestring := uint64ArrayToBitString(pop.Chromosomes[ind][copy])
	for chrom := 1; chrom < len(model.ChromosomeArms); chrom++ {
		pstart := model.ChromosomeArms[chrom][0][0]
		plen := model.ChromosomeArms[chrom][0][1]
		if pstart+plen <= len(genomestring) {
			psegment := genomestring[pstart : pstart+plen]
			blockCount += countBlocksInRange(psegment)
		}

		qstart := model.ChromosomeArms[chrom][1][0]
		qlen := model.ChromosomeArms[chrom][1][1]
		if qstart+qlen <= len(genomestring) {
			qsegment := genomestring[qstart : qstart+qlen]
			blockCount += countBlocksInRange(qsegment)
		}
	}
	return blockCount
}

func countBlocksInRange(genomesegment string) int {
	blocks := strings.Split(genomesegment, "0")
	blockCount := 0
	for _, block := range blocks {
		if len(block) > 0 {
			blockCount++
		}
	}
	return blockCount
}

func uint64ArrayToBitString(genomesegment []uint64) string {
	var bitString strings.Builder
	for _, value := range genomesegment {
		bitString.WriteString(fmt.Sprintf("%064b", value))
	}
	return bitString.String()
}

func inheritCentromeres(model *types.Model, pop *types.Pop, centsmask1 uint64, centsmask2 uint64, dad int, mom int, child int) {

	if pop.IndData[dad]["num_centromeres"] > 0 || pop.IndData[mom]["num_centromeres"] > 0 {
		pop.Centromeres[child] = make([]uint64, 2)
		for i := 0; i < len(model.ChromosomeArms); i++ {
			if centsmask1&(1<<i) == 0 {
				pop.Centromeres[child][0] |= (pop.Centromeres[dad][0] & (1 << i))
			} else {
				pop.Centromeres[child][0] |= (pop.Centromeres[dad][1] & (1 << i))
			}
			if centsmask2&(1<<i) == 0 {
				pop.Centromeres[child][1] |= (pop.Centromeres[mom][0] & (1 << i))
			} else {
				pop.Centromeres[child][1] |= (pop.Centromeres[mom][1] & (1 << i))
			}
		}

		centromereCount := countSetBitsSingleVar(pop.Centromeres[child][0])
		centromereCount += countSetBitsSingleVar(pop.Centromeres[child][1])
		pop.IndData[child]["num_centromeres"] = centromereCount
		if centromereCount == 0 {
			delete(pop.Centromeres, child)
		}
	}
}

func countSetBits(bits []uint64) int {
	count := 0
	for _, b := range bits {
		count += countSetBitsSingleVar(b)
	}
	return count
}

func countSetBitsSingleVar(n uint64) int {
	count := 0
	for n > 0 {
		count += int(n & 1)
		n >>= 1
	}
	return count
}
