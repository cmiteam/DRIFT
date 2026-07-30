package main

import (
	"drift/pkg/analysis"
	"drift/pkg/checkpoint"
	"drift/pkg/config"
	"drift/pkg/core"
	"drift/pkg/events"
	"drift/pkg/modules"
	"drift/pkg/simulation"
	"drift/pkg/utils"
	"drift/pkg/validation"
	"drift/pkg/visualization"
	"drift/pkg/webserver"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"time"
)

func writeProgress(resultsDir string, currentYear, totalYears, popSize, run, numRuns int, status string) {
	payload := map[string]interface{}{
		"current_generation": currentYear,
		"total_generations":  totalYears,
		"population_size":    popSize,
		"current_run":        run,
		"total_runs":         numRuns,
		"status":             status,
		"updated_at":         time.Now().Format(time.RFC3339),
	}
	data, err := json.Marshal(payload)
	if err != nil {
		return
	}
	_ = os.WriteFile(filepath.Join(resultsDir, "progress.json"), data, 0644)
}

// Main function does the following:
// 1. Parses command-line arguments
// 2. Initializes the model
// 3. Runs the model for the specified number of iterations
// 4. Saves the results
// 5. Prints the execution time

func main() {
	// Start the timer for execution time
	starttime := time.Now()

	// Parse the command-line arguments
	commands := config.ParseCommandLine()

	// Neutral-expectation validation harness (roadmap §6h). Standalone mode: drives
	// the real engine through a strictly neutral scenario over replicate seeds and
	// checks the emergent statistics (SFS, theta_W/theta_pi, Tajima's D, HWE) against
	// DRIFT's characterized neutral baseline. Exits non-zero if any check fails so it
	// can gate CI. Does not need a model.
	if commands.Validate {
		cfg := validation.DefaultNeutralConfig()
		fmt.Printf("Running neutral-expectation validation (§6h): R=%d replicates, N=%d, mu=%.3g, burn-in=%dy ...\n",
			cfg.Replicates, cfg.N, cfg.Mu, cfg.BurnInYears)
		res := validation.RunNeutral(cfg)
		tol := validation.CharacterizedTolerances()
		validation.PrintSummary(res, tol)

		if err := os.MkdirAll(commands.Results, 0755); err != nil {
			fmt.Fprintf(os.Stderr, "Error creating results dir: %v\n", err)
			os.Exit(1)
		}
		path := filepath.Join(commands.Results, "neutral_validation.csv")
		if err := validation.WriteAveragedReport(path, res, tol); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing validation report: %v\n", err)
			os.Exit(1)
		}
		fmt.Printf("Wrote validation report to %s\n", path)

		allPass := true
		for _, c := range res.Checks(tol) {
			if !c.Pass {
				allPass = false
			}
		}
		if !allPass {
			os.Exit(1)
		}
		return
	}

	// If web mode is enabled, start the web server
	if commands.WebMode {
		fmt.Println("Starting Drift web server...")
		err := webserver.Start()
		if err != nil {
			fmt.Fprintf(os.Stderr, "Web server error: %v\n", err)
			os.Exit(1)
		}
		return
	}

	// Initialize the model
	var model *core.Model
	var err error

	// Priority: User model > Base model > Legacy config
	if commands.Username != "" && commands.ModelName != "" {
		// Load user model using new ModelManager system
		model, err = config.LoadUserModel(commands.Username, commands.ModelName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error loading user model: %v\n", err)
			fmt.Fprintf(os.Stderr, "Try creating the model first with: drift -username=%s -model=%s -base-model=standard\n", commands.Username, commands.ModelName)
			os.Exit(1)
		}
		fmt.Printf("Loaded user model: %s/%s (based on %s)\n", commands.Username, commands.ModelName, model.BaseModelID)
	} else if commands.BaseModel != "" {
		// Create temporary model from base model (legacy support)
		model, err = config.InitializeModelWithBaseModel(commands.ConfigRoot, commands.BaseModel)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error initializing base model: %v\n", err)
			os.Exit(1)
		}
		fmt.Printf("Loaded base model: %s\n", commands.BaseModel)
	} else {
		// Legacy: load from config root (old system)
		model, err = config.InitializeModel(commands.ConfigRoot)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error initializing model: %v\n", err)
			os.Exit(1)
		}
		fmt.Println("Loaded model from config directory")
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error initializing model: %v\n", err)
		os.Exit(1)
	}

	// Fail fast if any selected phase module (birth/death/mating/seed/setup) isn't compiled in.
	if err := modules.ValidateStyles(model); err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		os.Exit(1)
	}

	// Setup user-specific directories if username is specified
	err = config.SetupUserDirectories(model)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error setting up user directories: %v\n", err)
		os.Exit(1)
	}

	// If output-dir was specified via command line, use it (overrides user directory)
	if commands.OutputDir != "" {
		model.ResultsDir = commands.OutputDir
	}

	// List the model's named saved states and exit, if requested.
	if commands.ListSaves {
		saves, err := checkpoint.List(model)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error listing saves: %v\n", err)
			os.Exit(1)
		}
		if len(saves) == 0 {
			fmt.Printf("No saved states in %s\n", checkpoint.SavesDir(model))
		} else {
			fmt.Printf("Saved states in %s:\n", checkpoint.SavesDir(model))
			for _, s := range saves {
				fmt.Printf("  %-20s year=%d run=%d n=%d scenario=%s seed=%d saved=%s\n",
					s.Name, s.Year, s.Run, s.NumLiving, s.Scenario, s.RNGSeed, s.SavedAt)
			}
		}
		return
	}

	// Setup directories
	err = config.SetupDirectories(model.ResultsDir)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		os.Exit(1)
	}

	// Create CSV file with headers (must be done after ResultsDir is finalized)
	err = simulation.SaveHeaders(model.ModelName, model.ResultsDir)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating results file: %v\n", err)
		os.Exit(1)
	}

	// VCF import (roadmap §6a): parse a real (or DRIFT-exported) phased VCF into
	// DRIFT structures and run the model's enabled analyses STATICALLY — no run loop
	// — so the same pipeline (Fst / f-stats / joint-SFS / §6g EHH-iHS-XP-EHH / §6h
	// validation / §6e dating) scores real 1000G/HGDP panels and simulated data
	// identically. The panel is loaded into BOTH the founder bitfield and the de-novo
	// pool consistently (see analysis.ImportVCF), so bitfield-reading and pool-reading
	// analyses both see it. Which analyses run is controlled by the model's existing
	// track_* params, exactly as for a live run. Runs before animation/profiling
	// setup, which a static import does not need.
	if commands.ImportVCF != "" {
		pol := analysis.DefaultImportPolicy()
		if commands.ImportPolarize != "" {
			pol.Polarize = commands.ImportPolarize
		}
		if commands.ImportPanel != "" {
			popMap, popNames, err := analysis.LoadPanelFile(commands.ImportPanel, 1)
			if err != nil {
				fmt.Fprintf(os.Stderr, "Error loading panel file: %v\n", err)
				os.Exit(1)
			}
			pol.PopMap = popMap
			fmt.Printf("Loaded panel %s: %d samples across %d populations %v\n",
				commands.ImportPanel, len(popMap), len(popNames), popNames)
		}
		pop := &core.Pop{
			IndData:      make(map[int][]int),
			Chromosomes:  make(map[int][][]uint64),
			Centromeres:  make(map[int][2][]uint64),
			IndMutations: make(map[int]map[int][]int),
			MutationPool: make(map[int]core.Mutation),
			MutationHist: make(map[int]int),
			Tracking:     make(map[string]int),
			AlleleFreqs:  make(map[int][]int16),
			MaleDB:       make(map[int]core.Ancestor),
			FemaleDB:     make(map[int]core.Ancestor),
		}
		stats, err := analysis.ImportVCF(commands.ImportVCF, model, pop, pol)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error importing VCF: %v\n", err)
			os.Exit(1)
		}
		fmt.Printf("\nImported %s\n  %d samples × %d biallelic-SNP sites (genome_bits=%d) across %d contig(s), %d deme(s)\n",
			stats.Path, stats.NumSamples, stats.NumSites, stats.GenomeBits, stats.NumContigs, stats.NumDemes)
		fmt.Printf("  polarization: %d AA-derived (%d flipped to REF) + %d REF-fallback | skipped %d multiallelic, %d indel, %d high-missing, %d monomorphic\n",
			stats.PolarizedAA, stats.PolarizedFlip, stats.PolarizedFall,
			stats.SkipMulti, stats.SkipIndel, stats.SkipMissing, stats.SkipMonomorph)
		if stats.UnphasedCalls > 0 {
			fmt.Printf("  NOTE: %d unphased genotype call(s) read in column order (DRIFT's strand model assumes phasing)\n", stats.UnphasedCalls)
		}
		model.FreeParameters["run"] = 0
		model.FreeParameters["year"] = 0
		runImportedAnalyses(model, pop)
		return
	}

	// Initialize animations if enabled
	animContainer, err := visualization.InitializeIfEnabled(model, commands.MapRoot)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error initializing animations: %v\n", err)
		os.Exit(1)
	}

	// Setup CPU profiling
	cleanup, err := config.StartCPUProfiling(commands.CpuProfile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		os.Exit(1)
	}
	defer cleanup()

	numRuns := int(model.Parameters["num_runs"])
	totalYears := int(model.Parameters["end_year"])

	// Establish a reproducible RNG seed. If rng_seed is unset (<= 0), derive one
	// from the clock and log it so the run can still be reproduced later by
	// setting rng_seed to that value. Each run in a multi-run gets a distinct,
	// deterministic seed (baseSeed + run).
	baseSeed := int64(model.Parameters["rng_seed"])
	if baseSeed <= 0 {
		baseSeed = time.Now().UnixNano()
		fmt.Printf("No rng_seed set; using generated seed %d (set rng_seed to this value to reproduce)\n", baseSeed)
	}

	// Optionally resume from a saved state instead of starting fresh. A named
	// state (-load-state) is resolved to a path in the model's saves/ directory;
	// -checkpoint-in takes a raw path.
	resumePath := commands.CheckpointIn
	if commands.LoadState != "" {
		resumePath = checkpoint.SavePath(model, commands.LoadState)
	}
	var resumeCP *checkpoint.Checkpoint
	if resumePath != "" {
		resumeCP, err = checkpoint.Load(resumePath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error loading saved state: %v\n", err)
			os.Exit(1)
		}
	}

	// maybeCheckpoint writes a checkpoint of the current state if enabled.
	maybeCheckpoint := func(pop *core.Pop) {
		if commands.CheckpointOut == "" {
			return
		}
		if err := checkpoint.Save(commands.CheckpointOut, model, pop); err != nil {
			log.Printf("Error writing checkpoint: %v", err)
		}
	}

	// runOne executes a single model run over pop, from startYear to totalYears.
	runOne := func(run, startYear int, pop *core.Pop) {
		model.FreeParameters["run"] = run
		writeProgress(model.ResultsDir, startYear, totalYears, len(pop.IndData), run, numRuns, "running")

		// Loop over the number years in each model run
		for year := startYear; year <= totalYears; year++ {
			model.FreeParameters["year"] = year
			// Advance any demographic scenario (§6c) for this year — fires splits,
			// migrates per the epoch matrix, and sets per-deme caps — before the
			// year's Birth/Mating/Death. No-op when no scenario is attached.
			if model.DemographyScheduler != nil {
				model.DemographyScheduler.Apply(model, pop)
			}
			// Fire any scriptable environmental events (§2) for this year — resets
			// the per-year event modifiers, then applies famines / migration pulses.
			// No-op when no events are scheduled (nil scheduler).
			if model.EventScheduler != nil {
				model.EventScheduler.Apply(model, pop)
			}
			if year >= int(model.Parameters["seed_year"]) && model.FreeParameters["seed"] == -1 {
				events.Seed(model, pop)
				model.FreeParameters["seed"] = 1
				fmt.Println("   Seeded population in year", year, "with seed style", int(model.Parameters["seed_style"]))
			}

			simulation.Birth(model, pop)
			simulation.Mating(model, pop)
			simulation.Death(model, pop)
			model.FreeParameters["last_pop_size"] = len(pop.IndData) // save pop size for future growth rate calculations

			// Escape clauses
			// Save and quit if population extinct
			if len(pop.IndData) <= 1 {
				simulation.Save(model, pop, animContainer)
				break
			}
			// Save and quit if genealo = pop size
			genealo := utils.CountGenealo(&pop.IndData)

			// Save and quit if genealo = 0
			if genealo == 0 {
				simulation.Save(model, pop, animContainer)
				break
			}
			if model.Parameters["track_map"] == 1 && year%int(model.Parameters["animation_save_interval"]) == 0 {
				err := visualization.AddAnimationFrames(model, pop, animContainer)
				if err != nil {
					log.Printf("Error adding animation frames: %v", err)
				}
			}
			// Normal save at save interval
			if year%int(model.Parameters["save_interval"]) == 0 {
				simulation.Save(model, pop, animContainer)
				writeProgress(model.ResultsDir, year, totalYears, len(pop.IndData), run, numRuns, "running")
			}
			// Capture an Ne-through-time window (§6d) at the ne_interval cadence.
			// Temporal (Waples) Ne between successive windows, contrasted end-of-run
			// with the §6h long-term coalescent Ne. Default off; a deterministic
			// full-population scan (no RNG) so a track_Ne run stays byte-identical to
			// a track_Ne-off run. ne_interval 0 falls back to save_interval; a
			// smaller value resolves short bottleneck windows.
			if model.Parameters["track_Ne"] == 1 {
				neInterval := int(model.Parameters["ne_interval"])
				if neInterval <= 0 {
					neInterval = int(model.Parameters["save_interval"])
				}
				if neInterval > 0 && year%neInterval == 0 {
					analysis.CaptureNe(model, pop)
				}
			}
			// Capture a genetic-load window (§6f) at the load_interval cadence.
			// Read-only full-population scan of the de-novo mutation pool (no RNG),
			// so a track_load run stays byte-identical to a track_load-off run.
			// load_interval 0 falls back to save_interval; a smaller value resolves
			// how load accumulates through a bottleneck / demographic epoch.
			if model.Parameters["track_load"] == 1 {
				loadInterval := int(model.Parameters["load_interval"])
				if loadInterval <= 0 {
					loadInterval = int(model.Parameters["save_interval"])
				}
				if loadInterval > 0 && year%loadInterval == 0 {
					analysis.CaptureLoad(model, pop)
				}
			}
			// Checkpoint at the configured interval.
			if commands.CheckpointInterval > 0 && year%commands.CheckpointInterval == 0 {
				maybeCheckpoint(pop)
			}

			// Honor a graceful stop-and-save request (e.g. from the GUI "Stop &
			// Save"): save this completed year's state under a name and finish the
			// run cleanly (end-of-run analyses/outputs still run below).
			if name, ok := checkpoint.StopRequested(model.ResultsDir); ok {
				checkpoint.ClearStop(model.ResultsDir)
				if name == "" {
					name = "stopped"
				}
				if path, err := checkpoint.SaveNamed(model, pop, name); err != nil {
					log.Printf("Error saving state on stop: %v", err)
				} else {
					fmt.Printf("Stop requested — saved run state %q (year %d, n=%d) to %s\n",
						name, year, len(pop.IndData), path)
				}
				break
			}
		}

		writeProgress(model.ResultsDir, totalYears, totalYears, len(pop.IndData), run, numRuns, "running")

		// Things to do at the end of a model run
		if model.Parameters["track_DNA"] == 1 {
			filename := fmt.Sprintf("%s/%s genome map.png", model.ResultsDir, model.ModelName)
			pixelSize := 4
			simulation.SaveGenomeMap(pop.Chromosomes, model.ChromosomeArms, filename, pixelSize, int(model.Parameters["NumBits"]))
		}
		if model.Parameters["track_map"] == 1 {
			//			print("Saving map...\n")
			err := visualization.SaveAllGIFs(animContainer, model.ResultsDir)
			if err != nil {
				log.Printf("Error saving animation: %v", err)
			}
		}
		if model.Parameters["save_detailed_SFS"] == 1 {
			analysis.SaveSFSTimeSeries(model, pop)
		}
		if model.Parameters["track_coalescence"] == 1 {
			yadamResult := analysis.FindYAdam(model, pop)
			mteveResult := analysis.FindMtEve(model, pop)
			fmt.Printf("\nY-Adam: ID %d, born year %d, %d generations back\n",
				yadamResult.YAdamID, yadamResult.YAdamBirthYear, yadamResult.GenerationsBack)
			fmt.Printf("Mt-Eve: ID %d, born year %d, %d generations back\n",
				mteveResult.MtEveID, mteveResult.MtEveBirthYear, mteveResult.GenerationsBack)
			fmt.Printf("Final MaleDB size: %d, FemaleDB size: %d\n", len(pop.MaleDB), len(pop.FemaleDB))
			fmt.Printf("Living individuals: %d\n", len(pop.IndData))
		}
		if model.Parameters["track_IBD"] == 1 {
			// Level 1 IBD prototype: co-inherited founder-derived segments.
			// Defaults off; requires track_DNA so chromosomes exist to compare.
			sampleSize := int(model.Parameters["ibd_sample_size"]) // 0 = use all individuals
			minBits := int(model.Parameters["ibd_min_bits"])       // "long" segment threshold
			ids := analysis.SampleIDs(pop, sampleSize)
			ibdResult := analysis.ExtractIBD(model, pop, ids, minBits)
			if err := analysis.SaveIBD(model, ibdResult); err != nil {
				log.Printf("Error saving IBD results: %v", err)
			}
		}
		if model.Parameters["export_VCF"] == 1 {
			// VCF export (§6a): standard interchange format so external tools
			// (PLINK, vcftools, ADMIXTURE) run on DRIFT output as on real data.
			// Requires track_DNA so chromosomes exist. vcf_sample_size 0 = all.
			sampleSize := int(model.Parameters["vcf_sample_size"])
			includeFixed := model.Parameters["vcf_include_fixed"] == 1
			ids := analysis.SampleIDs(pop, sampleSize)
			if _, err := analysis.ExportVCF(model, pop, ids, includeFixed); err != nil {
				log.Printf("Error exporting VCF: %v", err)
			}
		}
		if model.Parameters["track_Fst"] == 1 {
			// Fst between demes (§6b): Hudson's + Weir & Cockerham estimators.
			// Requires track_DNA (chromosomes) and num_demes >= 2 (structure to
			// measure). fst_sample_size 0 = all individuals.
			sampleSize := int(model.Parameters["fst_sample_size"])
			ids := analysis.SampleIDs(pop, sampleSize)
			fstResult := analysis.ComputeFst(model, pop, ids)
			if err := analysis.SaveFst(model, fstResult); err != nil {
				log.Printf("Error saving Fst results: %v", err)
			}
		}
		if model.Parameters["track_fstats"] == 1 {
			// f-statistics (§6b): f2, f3(admixture), f4/D (ABBA-BABA). Requires
			// track_DNA and num_demes >= 2. Which demes fill the f3/f4 slots is set
			// by fstats_demes (empty = auto-enumerate). fstats_sample_size 0 = all.
			sampleSize := int(model.Parameters["fstats_sample_size"])
			ids := analysis.SampleIDs(pop, sampleSize)
			fstatsResult := analysis.ComputeFStats(model, pop, ids)
			if err := analysis.SaveFStats(model, fstatsResult); err != nil {
				log.Printf("Error saving f-statistics: %v", err)
			}
		}
		if model.Parameters["track_joint_sfs"] == 1 {
			// Joint 2D/3D SFS (§6b): the object OoA demographic models are fit to.
			// Axis demes from joint_sfs_demes (empty = first eligible demes).
			// Requires track_DNA and num_demes >= 2. joint_sfs_sample_size 0 = all.
			sampleSize := int(model.Parameters["joint_sfs_sample_size"])
			ids := analysis.SampleIDs(pop, sampleSize)
			jsfs := analysis.ComputeJointSFS(model, pop, ids)
			if err := analysis.SaveJointSFS(model, jsfs); err != nil {
				log.Printf("Error saving joint SFS: %v", err)
			}
		}
		if model.Parameters["track_het_distance"] == 1 {
			// Heterozygosity vs distance-from-origin (§6c): the serial-founder
			// gradient (Ramachandran 2005). Per-deme He/Ho tagged with each deme's
			// serial-split rank + an He~rank slope, for the OoA-vs-Babel head-to-head.
			// Requires track_DNA and >= 2 eligible demes. het_distance_sample_size 0 = all.
			sampleSize := int(model.Parameters["het_distance_sample_size"])
			ids := analysis.SampleIDs(pop, sampleSize)
			hd := analysis.ComputeHetDistance(model, pop, ids)
			if err := analysis.SaveHetDistance(model, hd); err != nil {
				log.Printf("Error saving heterozygosity-distance results: %v", err)
			}
		}
		if model.Parameters["track_Ne"] == 1 {
			// Ne-through-time (§6d): write the per-window temporal (Waples) Ne series
			// (plus the caveated LD recent-Ne column) and a summary contrasting the
			// harmonic-mean temporal Ne with the census — the recent-history
			// counterpart to the §6h long-term coalescent Ne. Requires track_DNA so
			// the Chromosomes bitfield holds the standing variation the estimator
			// reads. A no-op when fewer than two windows were captured.
			if err := analysis.SaveNeTimeSeries(model, pop); err != nil {
				log.Printf("Error saving Ne time series: %v", err)
			}
		}
		if model.Parameters["track_dating"] == 1 {
			// Molecular-clock dating (§6e): date the autosomal coalescence (mutation-
			// pool theta_pi) under a grid of assumed mutation-rate / generation-time
			// values, and report the true genealogical Y-Adam/Mt-Eve dates alongside —
			// making the ~200 kya dates' assumption-dependence quantitative. Requires
			// track_mutations (the molecular signal). dating_sample_size 0 = all.
			sampleSize := int(model.Parameters["dating_sample_size"])
			ids := analysis.SampleLiving(pop, sampleSize)
			dating := analysis.ComputeDating(model, pop, ids)
			if err := analysis.SaveDating(model, dating); err != nil {
				log.Printf("Error saving dating results: %v", err)
			}
		}
		if model.Parameters["track_validation"] == 1 {
			// Neutral-expectation validation (§6h), single-realization diagnostic:
			// computes SFS/theta/Tajima's D/HWE from THIS run's mutation pool and
			// scores them against DRIFT's characterized neutral baseline. Requires
			// track_mutations. For the rigorous replicated pass/fail, run `-validate`
			// (pkg/validation) instead; a single run's Tajima's D is noisy (SD ~ 1).
			if err := analysis.SaveNeutralValidation(model, pop); err != nil {
				log.Printf("Error writing neutral validation report: %v", err)
			}
		}
		if model.Parameters["track_load"] == 1 {
			// Genetic load / DFE / meltdown characterization (§6f): write the per-
			// window load time series, the input-vs-realized (segregating vs fixed)
			// DFE, and a summary reporting the drift barrier and whether load is
			// rising (meltdown) or stationary (mutation-selection-drift balance).
			// Read-only over the de-novo mutation pool; requires track_mutations.
			if err := analysis.SaveLoadTimeSeries(model, pop); err != nil {
				log.Printf("Error saving genetic-load results: %v", err)
			}
		}
		if model.Parameters["track_haplostats"] == 1 {
			// Haplotype & selection statistics (§6g): EHH / iHS / XP-EHH / cM-binned
			// LD-decay over the de-novo mutation POOL (the LD-bearing, selection-
			// active substrate; the founder bitfield carries no sweep). Read-only;
			// requires track_mutations. haplostats_sample_size 0 = full-population
			// scan (no RNG); >0 subsamples via SampleLiving (opt-in, draws RNG).
			sampleSize := int(model.Parameters["haplostats_sample_size"])
			ids := analysis.SampleLiving(pop, sampleSize)
			cfg := analysis.DefaultHaploConfig()
			if v, ok := model.Parameters["haplostats_min_maf"]; ok && v > 0 {
				cfg.MinMAF = v
			}
			if v, ok := model.Parameters["haplostats_ehh_cutoff"]; ok && v > 0 {
				cfg.EHHCutoff = v
			}
			cfg.DemeA, cfg.DemeB = analysis.ParseDemePair(model.StringParam("haplostats_demes", ""))
			res := analysis.ComputeHaploStats(model, pop, ids, cfg)
			if err := analysis.SaveHaploStats(model, res); err != nil {
				log.Printf("Error saving haplotype-stats results: %v", err)
			}
			if model.Parameters["haplostats_export_hap"] == 1 {
				if err := analysis.ExportHapMap(model, pop, ids); err != nil {
					log.Printf("Error exporting hap/map: %v", err)
				}
			}
		}

		// Final checkpoint at the end of the run.
		maybeCheckpoint(pop)

		// Save the end-of-run state under a name, if requested, so it can be
		// loaded to start other runs (baseline / fork).
		if commands.SaveState != "" {
			name := commands.SaveState
			if numRuns > 1 {
				name = fmt.Sprintf("%s_run%d", name, run)
			}
			if path, err := checkpoint.SaveNamed(model, pop, name); err != nil {
				log.Printf("Error saving state %q: %v", name, err)
			} else {
				fmt.Printf("Saved run state %q (year %d, n=%d) to %s\n",
					name, model.FreeParameters["year"], len(pop.IndData), path)
			}
		}
	}

	if resumeCP != nil {
		// Resume a single continuation from the checkpoint (exact, or a fork if
		// fork-seed is set). Model config was loaded from files above; overlay the
		// saved counters + RNG position and continue from the next year.
		if err := checkpoint.Restore(model, resumeCP, int64(commands.ForkSeed)); err != nil {
			fmt.Fprintf(os.Stderr, "Error restoring checkpoint: %v\n", err)
			os.Exit(1)
		}
		mode := "exact resume"
		if commands.ForkSeed != 0 {
			mode = fmt.Sprintf("fork (fork_seed=%d)", commands.ForkSeed)
		}
		fmt.Printf("\nResuming run %d from checkpoint at year %d — %s\n", resumeCP.Run, resumeCP.Year, mode)
		runOne(resumeCP.Run, resumeCP.Year+1, resumeCP.Pop)
	} else {
		// Loop over the number of model runs
		for run := 1; run <= numRuns; run++ {
			runSeed := baseSeed + int64(run)
			utils.SeedRNG(runSeed)
			fmt.Printf("\nRun %d (rng_seed=%d)\n", run, runSeed)
			pop := config.InitializePop(model)
			runOne(run, 0, pop)
		}
	}

	// End of model runs
	elapsed := time.Since(starttime)
	fmt.Printf("Execution time: %s\n", elapsed)
	fmt.Print("\a")
}

// runImportedAnalyses runs the snapshot (single-generation) analyses on an imported
// VCF panel (roadmap §6a), gated on the same track_* params as a live run's
// end-of-run block. It deliberately omits the time-series analyses (§6d Ne, §6f
// load) — those need multiple captured windows and are meaningless for a static
// panel — and the genealogy analyses (coalescence/IBD, which need a tracked family
// tree the panel has none of). Bitfield-reading (Fst/f-stats/joint-SFS/het-distance)
// and pool-reading (§6h/§6e/§6g) analyses both apply because ImportVCF populated both.
func runImportedAnalyses(model *core.Model, pop *core.Pop) {
	if model.Parameters["export_VCF"] == 1 {
		// Round-trip re-export: read the imported panel back out as a normalized VCF.
		sampleSize := int(model.Parameters["vcf_sample_size"])
		includeFixed := model.Parameters["vcf_include_fixed"] == 1
		ids := analysis.SampleIDs(pop, sampleSize)
		if _, err := analysis.ExportVCF(model, pop, ids, includeFixed); err != nil {
			log.Printf("Error exporting VCF: %v", err)
		}
	}
	if model.Parameters["track_Fst"] == 1 {
		ids := analysis.SampleIDs(pop, int(model.Parameters["fst_sample_size"]))
		if err := analysis.SaveFst(model, analysis.ComputeFst(model, pop, ids)); err != nil {
			log.Printf("Error saving Fst results: %v", err)
		}
	}
	if model.Parameters["track_fstats"] == 1 {
		ids := analysis.SampleIDs(pop, int(model.Parameters["fstats_sample_size"]))
		if err := analysis.SaveFStats(model, analysis.ComputeFStats(model, pop, ids)); err != nil {
			log.Printf("Error saving f-statistics: %v", err)
		}
	}
	if model.Parameters["track_joint_sfs"] == 1 {
		ids := analysis.SampleIDs(pop, int(model.Parameters["joint_sfs_sample_size"]))
		if err := analysis.SaveJointSFS(model, analysis.ComputeJointSFS(model, pop, ids)); err != nil {
			log.Printf("Error saving joint SFS: %v", err)
		}
	}
	if model.Parameters["track_het_distance"] == 1 {
		ids := analysis.SampleIDs(pop, int(model.Parameters["het_distance_sample_size"]))
		if err := analysis.SaveHetDistance(model, analysis.ComputeHetDistance(model, pop, ids)); err != nil {
			log.Printf("Error saving heterozygosity-distance results: %v", err)
		}
	}
	if model.Parameters["track_dating"] == 1 {
		ids := analysis.SampleLiving(pop, int(model.Parameters["dating_sample_size"]))
		if err := analysis.SaveDating(model, analysis.ComputeDating(model, pop, ids)); err != nil {
			log.Printf("Error saving dating results: %v", err)
		}
	}
	if model.Parameters["track_validation"] == 1 {
		if err := analysis.SaveNeutralValidation(model, pop); err != nil {
			log.Printf("Error writing neutral validation report: %v", err)
		}
	}
	if model.Parameters["track_haplostats"] == 1 {
		ids := analysis.SampleLiving(pop, int(model.Parameters["haplostats_sample_size"]))
		cfg := analysis.DefaultHaploConfig()
		if v, ok := model.Parameters["haplostats_min_maf"]; ok && v > 0 {
			cfg.MinMAF = v
		}
		if v, ok := model.Parameters["haplostats_ehh_cutoff"]; ok && v > 0 {
			cfg.EHHCutoff = v
		}
		cfg.DemeA, cfg.DemeB = analysis.ParseDemePair(model.StringParam("haplostats_demes", ""))
		res := analysis.ComputeHaploStats(model, pop, ids, cfg)
		if err := analysis.SaveHaploStats(model, res); err != nil {
			log.Printf("Error saving haplotype-stats results: %v", err)
		}
		if model.Parameters["haplostats_export_hap"] == 1 {
			if err := analysis.ExportHapMap(model, pop, ids); err != nil {
				log.Printf("Error exporting hap/map: %v", err)
			}
		}
	}
}
