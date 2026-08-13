package analysis

// TMR-K-A output (TMR4A.md W5 / §4). Writes two CSVs and prints a summary.
//
//	<model>_tmrka_loci_run%d_year%d.csv — one row per focal locus: the §4 schema.
//	  This is the table that joins against ARGweaver's inferred TMR-K-A (W7) on
//	  `locus` to make the headline truth-vs-inference figure.
//
//	<model>_tmrka_trajectory_run%d_year%d.csv — the full lineage-count-versus-year
//	  trajectory, one row per coalescence per locus. TMR-K-A for ANY K is a lookup
//	  against this, which is what makes K a parameter rather than a design centre.
//
// WHAT THE OUTCOME COLUMN IS FOR. TMR4A.md §6 risk 5: "reporting a censored walk as a
// date would reproduce exactly the error the module exists to expose." So the date
// columns are written EMPTY, not zero and not a sentinel number, whenever the walk did
// not actually reach that many lineages — an empty cell cannot be averaged by accident,
// and the outcome column says why in words. The stdout summary breaks the loci down by
// outcome class for the same reason.

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
)

// SaveTMRKA writes the per-locus and trajectory CSVs and prints a summary. A no-op
// (nil error) when the walk produced nothing.
func SaveTMRKA(model *core.Model, res *TMRKAResult) error {
	if res == nil || len(res.Loci) == 0 {
		fmt.Println("TMR-K-A: no focal loci walked (requires track_pedigree + track_arg + arg_loci); skipped")
		return nil
	}
	if err := saveTMRKALoci(model, res); err != nil {
		return err
	}
	if err := saveTMRKATrajectory(model, res); err != nil {
		return err
	}

	// Summarise by outcome class. Averaging only over the loci that actually coalesced
	// is the point: a mean over censored loci treated as zeros would be meaningless.
	counts := map[string]int{}
	sumYears, nYears := 0, 0
	minFinal, maxFinal, sumFinal := 1<<62, 0, 0
	for i := range res.Loci {
		l := &res.Loci[i]
		counts[l.Outcome]++
		if l.HasTMRKA {
			sumYears += res.SampleYear - l.TMRKAYear
			nYears++
		}
		if l.FinalLineages < minFinal {
			minFinal = l.FinalLineages
		}
		if l.FinalLineages > maxFinal {
			maxFinal = l.FinalLineages
		}
		sumFinal += l.FinalLineages
	}

	fmt.Printf("TMR-%d-A: %d focal loci, n=%d individuals (%d lineages) — %d coalesced, %d censored at founding, %d with no coalescence\n",
		res.K, len(res.Loci), res.NumSampled, 2*res.NumSampled,
		counts[OutcomeCoalesced], counts[OutcomeCensored], counts[OutcomeNoCoal])
	if nYears > 0 {
		meanAge := float64(sumYears) / float64(nYears)
		fmt.Printf("  mean TMR-%d-A depth over coalesced loci: %.1f years (%.1f generations at %.0f yr/gen)\n",
			res.K, meanAge, meanAge/res.GenerationTime, res.GenerationTime)
	}
	// The headline falsifiable number (TMR4A.md §0): distinct founder haplotypes still
	// present in the sample. Compare against the founding couple's transmissible allele
	// count, NOT against an inferred date.
	fmt.Printf("  founder haplotypes surviving per locus: min %d, mean %.2f, max %d\n",
		minFinal, float64(sumFinal)/float64(len(res.Loci)), maxFinal)
	if res.MissingARG > 0 {
		fmt.Printf("  %d loci skipped: no strand provenance recorded\n", res.MissingARG)
	}
	return nil
}

// blankUnless writes an empty cell rather than a number when the quantity does not
// exist. See the file comment: a censored walk has no date, and an empty cell is the
// only representation of that which survives being loaded into a plotting script.
func blankUnless(ok bool, v int) string {
	if !ok {
		return ""
	}
	return strconv.Itoa(v)
}

func saveTMRKALoci(model *core.Model, res *TMRKAResult) error {
	path := tmrkaPath(model, "loci")
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create TMR-K-A loci file: %v", err)
	}
	defer file.Close()

	w := csv.NewWriter(file)
	defer w.Flush()

	if err := w.Write([]string{
		"locus", "k", "outcome", "tmrka_year", "tmrka_gens_ago", "tmrca_year",
		"sample_year", "initial_lineages", "lineages_at_founding", "coalescences",
		"founder_strands_surviving", "created_alleles_surviving",
		"mutational_diffs", "created_diffs",
	}); err != nil {
		return err
	}

	for i := range res.Loci {
		l := &res.Loci[i]
		gensAgo := ""
		if l.HasTMRKA && res.GenerationTime > 0 {
			gensAgo = fmt.Sprintf("%.2f", float64(res.SampleYear-l.TMRKAYear)/res.GenerationTime)
		}
		if err := w.Write([]string{
			strconv.Itoa(l.Position),
			strconv.Itoa(l.K),
			l.Outcome,
			blankUnless(l.HasTMRKA, l.TMRKAYear),
			gensAgo,
			blankUnless(l.HasTMRCA, l.TMRCAYear),
			strconv.Itoa(l.SampleYear),
			strconv.Itoa(l.InitialLineages),
			strconv.Itoa(l.FinalLineages),
			strconv.Itoa(l.Coalescences),
			strconv.Itoa(l.FounderStrandsSurviving),
			strconv.Itoa(l.CreatedAllelesSurviving),
			strconv.Itoa(l.MutationalDiffs),
			strconv.Itoa(l.CreatedDiffs),
		}); err != nil {
			return err
		}
	}
	return nil
}

func saveTMRKATrajectory(model *core.Model, res *TMRKAResult) error {
	path := tmrkaPath(model, "trajectory")
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create TMR-K-A trajectory file: %v", err)
	}
	defer file.Close()

	w := csv.NewWriter(file)
	defer w.Flush()

	if err := w.Write([]string{"locus", "year", "gens_ago", "lineages"}); err != nil {
		return err
	}
	for i := range res.Loci {
		l := &res.Loci[i]
		// Row 0 is the sample itself, so the trajectory is a complete step function
		// from 2n lineages down to whatever survived, not just the change points.
		rows := append([]TMRKAPoint{{Year: l.SampleYear, Lineages: l.InitialLineages}}, l.Trajectory...)
		for _, p := range rows {
			gens := ""
			if res.GenerationTime > 0 {
				gens = fmt.Sprintf("%.2f", float64(res.SampleYear-p.Year)/res.GenerationTime)
			}
			if err := w.Write([]string{
				strconv.Itoa(l.Position),
				strconv.Itoa(p.Year),
				gens,
				strconv.Itoa(p.Lineages),
			}); err != nil {
				return err
			}
		}
	}
	return nil
}

func tmrkaPath(model *core.Model, kind string) string {
	run := model.FreeParameters["run"]
	year := model.FreeParameters["year"]
	return fmt.Sprintf("%s/%s_tmrka_%s_run%d_year%d.csv", model.ResultsDir, model.ModelName, kind, run, year)
}
