package analysis

// The W7 scaling rationale, computed rather than asserted — TMR4A.md W7 and §6 risk 1.
//
// WHAT THIS EXISTS TO SETTLE. §6 risk 1 reduces the old open-ended scaling problem to a
// documentation problem, and names exactly three things that must be defended in writing
// before any headline run: the bp-per-bit convention, mutations per individual per
// generation, and crossovers per meiosis. This module computes all three off the live
// model and writes them out with the derived quantities, so the rationale is a MEASURED
// artifact that ships with the data rather than a paragraph someone has to trust.
//
// WHY THE RATIO IS THE THING. ARGweaver takes -m (mutations per site per generation) and
// -r (recombinations per site per generation). DRIFT's genome is a bit array, not
// nucleotides, so "per site" only means something once a bp-per-bit convention B is
// declared. But per W7:
//
//	mutations are a per-INDIVIDUAL Poisson count spread uniformly over genome_bits;
//	crossovers per MEIOSIS come from the separate cM column.
//
// Neither anchor mentions B. So declaring B fixes the region length (genome_bits × B) and
// divides BOTH per-bit rates by B, leaving
//
//	θ/ρ = μ/r = (mutations per gamete) / (crossovers per meiosis)
//
// invariant. B buys resolution only; the ratio is set by two per-generation counts and is
// the thing that has to match human values. Human: μ ≈ 1.25e-8 and r ≈ 1.2e-8 per bp per
// generation ⇒ θ/ρ ≈ 1.04, or equivalently ≈ 38.8 mutations per gamete against ≈ 34
// crossovers per meiosis.
//
// TWO PLACES THE OBVIOUS COMPUTATION IS WRONG, both found here rather than argued:
//
//  1. **`mu` is per INDIVIDUAL, crossovers are per GAMETE.** GenerateNewMutations draws
//     Poisson(rate) once per newborn and puts each mutation on ONE of its two strands, so
//     the per-gamete mutation count is mu/2. Comparing mu directly against crossovers per
//     meiosis overstates θ/ρ by exactly 2.
//
//  2. **The realized crossover rate is NOT the cM map's rate.** createMaskFromMap gives
//     every mapped chromosome one OBLIGATE crossover plus Poisson((totalCM −
//     recomb_obligate_cM)/100) extras. On a genome whose arms carry only a few cM each —
//     which is every shipped DRIFT chromosome map — the obligate crossover dominates
//     completely and the realized rate can exceed the map's nominal rate by an order of
//     magnitude. Telling ARGweaver a ρ derived from the cM column would then be wrong by
//     that factor. Both numbers are reported and the report says which one was used.
//
// PER-REGION ρ, NOT A GENOME-WIDE AVERAGE. ARGweaver runs one region at a time, and the
// obligate crossover makes a short chromosome recombine far more per bp than a long one.
// So ρ is computed per chromosome from that chromosome's own realized crossover count. μ
// is genome-wide by construction (positions are drawn uniformly over genome_bits), which
// is also what ARGweaver's constant-μ assumption wants.

import (
	"drift/pkg/core"
	"drift/pkg/utils"
	"fmt"
	"io"
	"math"
)

// Human reference anchors, per bp per generation, used only for the comparison line in
// the report. Sources are the same standard values §1A quotes.
const (
	HumanMuPerBp  = 1.25e-8
	HumanRhoPerBp = 1.2e-8
	// Genome-wide per-generation counts implied by those rates over a 3.1 Gb haploid
	// genome and a ~3400 cM sex-averaged map.
	HumanMutationsPerGamete   = HumanMuPerBp * 3.1e9
	HumanCrossoversPerMeiosis = 34.0
)

// HumanThetaOverRho is the ratio the two DRIFT anchors have to reproduce.
const HumanThetaOverRho = HumanMuPerBp / HumanRhoPerBp

// RegionRates holds one chromosome's exported recombination anchor.
type RegionRates struct {
	Chrom   int
	Bits    int     // genome bits spanned by the chromosome
	BpLen   int     // Bits × BpPerBit — the ARGweaver REGION length
	TotalCM float64 // the cM map's genetic length for this chromosome

	// CrossoversPerMeiosis is the REALIZED expected count for this chromosome:
	// 1 obligate + max(0, TotalCM − recomb_obligate_cM)/100 extras.
	CrossoversPerMeiosis float64
	// CrossoversNominal is what the cM column alone implies (TotalCM/100). The gap
	// between the two is the obligate-crossover inflation.
	CrossoversNominal float64

	RhoPerBp float64 // CrossoversPerMeiosis / BpLen — what arg-sample -r must be told

	// cM per bit at the least- and most-recombining arm. ARGweaver assumes a constant
	// rate across the region, so a wide spread here is a modelling mismatch worth
	// knowing about before the run, not after.
	CMPerBitMin float64
	CMPerBitMax float64
}

// RateAnchors is the whole scaling rationale for one export.
type RateAnchors struct {
	// The two per-generation counts §6 risk 1 says must be defended.
	MutationsPerIndividual float64 // Σ class rates — expected de-novo per newborn
	MutationsPerGamete     float64 // half of it: each mutation lands on one strand
	CrossoversPerMeiosis   float64 // realized genome-wide, obligate included
	CrossoversNominal      float64 // what the cM map alone implies
	ObligateCrossovers     float64 // one per mapped chromosome; the inflation's source

	// The declared convention and what follows from it.
	GenomeBits    int
	BpPerBit      int
	RegionBpTotal int
	MuPerBp       float64 // genome-wide; uniform by construction

	// The scale-invariant quantity under test.
	ThetaOverRho float64
	HumanRatio   float64
	RatioVsHuman float64 // ThetaOverRho / HumanRatio; 1.0 = on target
	RecombModel  string
	LinkageModel string
	Regions      []RegionRates
	Warnings     []string
}

// ComputeRateAnchors derives the whole rationale from the model. Pure: no RNG, no
// mutation of the model beyond the cached genetic map EnsureGeneticMap already builds.
func ComputeRateAnchors(model *core.Model, bpPerBit int, contigs []vcfContig) *RateAnchors {
	a := &RateAnchors{
		GenomeBits:   model.FreeParameters["genome_bits"],
		BpPerBit:     bpPerBit,
		HumanRatio:   HumanThetaOverRho,
		RecombModel:  model.StringParam("recombination_model", "legacy"),
		LinkageModel: model.StringParam("linkage_model", "legacy"),
	}

	// Anchor 1 — mutations per individual per generation, summed over the class
	// spectrum. With mutation_classes unset this is the single point class whose rate is
	// the model's global mu, i.e. exactly the historical scalar.
	for _, class := range utils.MutationClasses(model) {
		a.MutationsPerIndividual += class.Rate
	}
	a.MutationsPerGamete = a.MutationsPerIndividual / 2

	// Anchor 2 — crossovers per meiosis, per chromosome and summed, from the same
	// arithmetic createMaskFromMap actually performs.
	obligateCM := 50.0
	if v, ok := model.Parameters["recomb_obligate_cM"]; ok {
		obligateCM = v
	}
	gm := utils.EnsureGeneticMap(model)

	for _, ct := range contigs {
		r := RegionRates{
			Chrom: ct.chrom,
			Bits:  ct.span,
			BpLen: ct.span * bpPerBit,
		}
		if gm != nil && gm.Has(ct.chrom) {
			r.TotalCM = gm.TotalCM(ct.chrom)
			r.CMPerBitMin, r.CMPerBitMax = armRateSpread(gm, ct)
		}
		if r.TotalCM > 0 {
			r.CrossoversNominal = r.TotalCM / 100
			r.CrossoversPerMeiosis = 1 + math.Max(0, r.TotalCM-obligateCM)/100
			a.ObligateCrossovers++
		}
		if r.BpLen > 0 {
			r.RhoPerBp = r.CrossoversPerMeiosis / float64(r.BpLen)
		}
		a.CrossoversPerMeiosis += r.CrossoversPerMeiosis
		a.CrossoversNominal += r.CrossoversNominal
		a.RegionBpTotal += r.BpLen
		a.Regions = append(a.Regions, r)
	}

	if a.RegionBpTotal > 0 {
		a.MuPerBp = a.MutationsPerGamete / float64(a.RegionBpTotal)
	}
	if a.CrossoversPerMeiosis > 0 {
		a.ThetaOverRho = a.MutationsPerGamete / a.CrossoversPerMeiosis
	}
	if a.HumanRatio > 0 {
		a.RatioVsHuman = a.ThetaOverRho / a.HumanRatio
	}

	a.Warnings = a.check()
	return a
}

// armRateSpread returns the least and greatest cM-per-bit across the chromosome's arms,
// read straight off the map's own cumulative function so the arm-cM fallback rule is not
// duplicated here and cannot drift from BuildGeneticMap's.
func armRateSpread(gm *core.GeneticMap, ct vcfContig) (lo, hi float64) {
	lo = math.Inf(1)
	for _, arm := range ct.arms {
		start, length := arm[0], arm[1]
		if length <= 0 {
			continue
		}
		cm := gm.CumCM(ct.chrom, start+length) - gm.CumCM(ct.chrom, start)
		rate := cm / float64(length)
		if rate < lo {
			lo = rate
		}
		if rate > hi {
			hi = rate
		}
	}
	if math.IsInf(lo, 1) {
		lo = 0
	}
	return lo, hi
}

// check collects everything that would make the exported data an artifact rather than a
// model. These are stated as warnings and not as refusals: a deliberately off-anchor run
// is a legitimate experiment (W9 sweeps exactly this kind of thing), but it must never be
// off-anchor by accident.
func (a *RateAnchors) check() []string {
	var w []string

	if a.RecombModel != "map" {
		w = append(w, fmt.Sprintf("recombination_model=%q, not \"map\": the legacy kernel puts one "+
			"contiguous segment per chromosome, so a locus's genealogy depends on where in the arm "+
			"it sits (TMR4A.md §8a measured >2x). ARGweaver assumes recombination is "+
			"position-homogeneous; this export would measure the engine's mask shape.", a.RecombModel))
	}
	if a.LinkageModel != "arm" {
		w = append(w, fmt.Sprintf("linkage_model=%q, not \"arm\": the de-novo pool and the founder "+
			"bitfield anti-segregate, so the exported haplotypes are an artifact of mixed polarity "+
			"(TMR4A.md §8c).", a.LinkageModel))
	}
	if a.MutationsPerIndividual <= 0 {
		w = append(w, "mutations per individual per generation is 0: there is no mutational signal "+
			"to infer a genealogy from. Set mu (and track_mutations).")
	}
	if a.CrossoversPerMeiosis <= 0 {
		w = append(w, "crossovers per meiosis is 0: every site is completely linked, so there is one "+
			"tree for the whole region and no ARG to infer.")
	}
	if a.RatioVsHuman > 0 && (a.RatioVsHuman < 0.5 || a.RatioVsHuman > 2) {
		w = append(w, fmt.Sprintf("theta/rho = %.4g is %.3gx the human value %.4g. The mainstream "+
			"comparison arm depends on this ratio being right — at %.3gx, each non-recombining block "+
			"carries %.3gx the human number of segregating sites and the per-block genealogy is "+
			"correspondingly better or worse determined than the published analysis's (§1A).",
			a.ThetaOverRho, a.RatioVsHuman, a.HumanRatio, a.RatioVsHuman, a.RatioVsHuman))
	}
	if a.CrossoversNominal > 0 && a.CrossoversPerMeiosis/a.CrossoversNominal > 1.5 {
		w = append(w, fmt.Sprintf("the obligate crossover dominates: %.4g realized crossovers per "+
			"meiosis against %.4g the cM map implies (%.3gx), because %v mapped chromosomes each get "+
			"one guaranteed crossover regardless of their genetic length. The exported rho uses the "+
			"REALIZED count, which is correct — but the cM column is no longer describing the "+
			"recombination the engine performs. Lengthen the map or lower recomb_obligate_cM to "+
			"close the gap.", a.CrossoversPerMeiosis, a.CrossoversNominal,
			a.CrossoversPerMeiosis/a.CrossoversNominal, int(a.ObligateCrossovers)))
	}
	if a.BpPerBit > 100 {
		w = append(w, fmt.Sprintf("bp per bit is %d, coarser than the ~100 bp W7 targets. A local "+
			"tree spans hundreds of bp to ~1 kb (§1A), so at this resolution breakpoints inside a "+
			"bit are unresolvable. Raise resolution by scaling the chromosome CSV's Start/Length — "+
			"that raises genome_bits without touching either rate anchor, so theta/rho is unchanged.",
			a.BpPerBit))
	}
	for _, r := range a.Regions {
		if r.CMPerBitMin > 0 && r.CMPerBitMax/r.CMPerBitMin > 2 {
			w = append(w, fmt.Sprintf("chromosome %d's arms differ in recombination rate by %.3gx "+
				"(%.4g to %.4g cM/bit). ARGweaver is told ONE rho for the region, so this "+
				"heterogeneity is invisible to it and will be absorbed into the inferred genealogy.",
				r.Chrom, r.CMPerBitMax/r.CMPerBitMin, r.CMPerBitMin, r.CMPerBitMax))
		}
	}
	return w
}

// WriteRationale renders the anchors as the standalone written justification W7 asks for
// ("a written scaling rationale that stands on its own"). Prose, not a CSV, because its
// audience is a reader deciding whether to believe the run.
func (a *RateAnchors) WriteRationale(w io.Writer, model *core.Model) {
	fmt.Fprintln(w, "ARGweaver export — scaling rationale (TMR4A.md W7, §6 risk 1)")
	fmt.Fprintf(w, "model=%s run=%d year=%d\n\n",
		model.ModelName, model.FreeParameters["run"], model.FreeParameters["year"])

	fmt.Fprintln(w, "THE BP-PER-BIT CONVENTION")
	fmt.Fprintf(w, "  DRIFT's genome is a bit array, not nucleotides. This export declares\n")
	fmt.Fprintf(w, "    1 genome bit = %d bp\n", a.BpPerBit)
	fmt.Fprintf(w, "  so the %d-bit genome is exported as %d bp of sequence.\n", a.GenomeBits, a.RegionBpTotal)
	fmt.Fprintln(w, "  This is a DECLARATION, not a measurement. It sets resolution and the absolute")
	fmt.Fprintln(w, "  per-bp rates below; it does NOT affect theta/rho, because it divides both.")
	fmt.Fprintln(w)

	fmt.Fprintln(w, "THE TWO PER-GENERATION RATE ANCHORS")
	fmt.Fprintf(w, "  mutations per individual per generation   %12.6g   (sum over the class spectrum)\n",
		a.MutationsPerIndividual)
	fmt.Fprintf(w, "  mutations per GAMETE                      %12.6g   (each lands on one strand)\n",
		a.MutationsPerGamete)
	fmt.Fprintf(w, "  crossovers per meiosis, realized          %12.6g   (%.0f obligate + Poisson extras)\n",
		a.CrossoversPerMeiosis, a.ObligateCrossovers)
	fmt.Fprintf(w, "  crossovers per meiosis, cM map alone      %12.6g\n", a.CrossoversNominal)
	fmt.Fprintln(w, "  The REALIZED count is what the engine performs and what the exported rho uses.")
	fmt.Fprintln(w)

	fmt.Fprintln(w, "THE SCALE-INVARIANT QUANTITY")
	fmt.Fprintf(w, "  theta/rho  = (mutations per gamete) / (crossovers per meiosis) = %.6g\n", a.ThetaOverRho)
	fmt.Fprintf(w, "  human      = %.4g per bp / %.4g per bp                          = %.6g\n",
		HumanMuPerBp, HumanRhoPerBp, a.HumanRatio)
	fmt.Fprintf(w, "  ratio      = %.4g x human\n", a.RatioVsHuman)
	fmt.Fprintf(w, "  human reference counts: %.4g mutations per gamete, %.4g crossovers per meiosis\n",
		HumanMutationsPerGamete, HumanCrossoversPerMeiosis)
	fmt.Fprintln(w)

	fmt.Fprintln(w, "WHAT arg-sample MUST BE TOLD")
	fmt.Fprintf(w, "  --mutrate    %.6g   (per site per generation, genome-wide uniform)\n", a.MuPerBp)
	fmt.Fprintln(w, "  --recombrate  per region, since the obligate crossover makes it chromosome-specific:")
	fmt.Fprintf(w, "    %-6s %10s %10s %12s %14s %14s\n", "chrom", "bits", "bp", "cM", "xovers/meiosis", "rho/bp")
	for _, r := range a.Regions {
		fmt.Fprintf(w, "    %-6d %10d %10d %12.4g %14.6g %14.6g\n",
			r.Chrom, r.Bits, r.BpLen, r.TotalCM, r.CrossoversPerMeiosis, r.RhoPerBp)
	}
	fmt.Fprintln(w)

	fmt.Fprintln(w, "ENGINE SETTINGS THIS EXPORT DEPENDS ON")
	fmt.Fprintf(w, "  recombination_model = %q  (must be \"map\"; see TMR4A.md §8a)\n", a.RecombModel)
	fmt.Fprintf(w, "  linkage_model       = %q  (must be \"arm\"; see TMR4A.md §8c)\n", a.LinkageModel)
	fmt.Fprintln(w)

	if len(a.Warnings) == 0 {
		fmt.Fprintln(w, "NO WARNINGS — the anchors above are internally consistent and near human values.")
		return
	}
	fmt.Fprintf(w, "WARNINGS (%d) — read these before treating any inference on this data as a result\n",
		len(a.Warnings))
	for i, msg := range a.Warnings {
		fmt.Fprintf(w, "  %d. %s\n", i+1, msg)
	}
}
