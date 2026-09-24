package config

import (
	"flag"
)

// Config holds all the command-line configuration options
type Config struct {
	ConfigRoot string
	MapRoot    string
	Results    string
	CpuProfile string
	WebMode    bool   // Enable web server mode
	OutputDir  string // Output directory for results (GUI mode)
	BaseModel  string // Base model to use (e.g., "standard", "flood")
	Username   string // Username for model loading
	ModelName  string // Model name to load

	CheckpointOut      string // path to write run-state checkpoints (empty = disabled)
	CheckpointInterval int    // years between checkpoints (0 = only at end of run)
	CheckpointIn       string // path to a checkpoint to resume from (empty = fresh run)
	ForkSeed           int    // if resuming and non-zero, re-seed the RNG (divergent fork)

	SaveState string // name to save the end-of-run state under (in the model's saves/ dir)
	LoadState string // name of a saved state to start from (resume/fork)
	ListSaves bool   // list named saves for the model and exit

	Validate bool // run the §6h neutral-expectation validation harness and exit

	// VCF import (roadmap §6a). ImportVCF is the path to a phased VCF (.vcf or
	// .vcf.gz) to parse into DRIFT structures and analyze statically (no run loop),
	// so the same pipeline runs on real 1000G/HGDP and simulated data. ImportPanel
	// is an optional sample→population panel file for deme assignment (Fst/f-stats).
	// ImportPolarize is the ancestral/derived policy ("aa" = AA INFO with REF
	// fallback, default; "ref" = always REF-ancestral). Empty ImportVCF = no import.
	ImportVCF      string
	ImportPanel    string
	ImportPolarize string

	// ARGweaver output parsing (TMR4A.md W7, the parser half). ARGweaverTrees is
	// the path to `arg-summarize -a <run>.bed.gz --tree` output; ARGweaverTruth is
	// the DRIFT <model>_tmrka_loci_<run>_<year>.csv the inference is judged against.
	// Standalone: no model is loaded and no run happens. Empty ARGweaverTrees = off.
	ARGweaverTrees     string
	ARGweaverTruth     string
	ARGweaverOut       string
	ARGweaverK         int
	ARGweaverBpPerBit  int
	ARGweaverGenTime   float64
	ARGweaverMinSample int
	ARGweaverMaxSample int

	// Raw ARGweaver .smc scanning — the genome-wide TMR-K-A distribution and its
	// spread across MCMC samples, for real data where there is no truth to join
	// against (Rasmussen et al. 2014). Comma-separated files and/or directories;
	// directories are walked for *.smc(.gz). Empty = off.
	ARGweaverSMC         string
	ARGweaverSMCOut      string
	ARGweaverSMCIterPath bool

	// Haplotype clustering + the sampling bootstrap (TMR4A.md W9c stage 1).
	// ClusterSites is an exported `.sites` file; the clustering runs on the window
	// 1..ClusterWindow bp so it measures the same object arg-sample was given via
	// `--region`. Standalone: no model, no run, no ARGweaver. Empty = off.
	ClusterSites     string
	ClusterOut       string
	ClusterWindow    int
	ClusterN         int
	ClusterDraws     int
	ClusterSeed      int64
	ClusterThreshold float64
	ClusterMetric    string
	ClusterMinHeight float64
}

// Default values for command-line parameters
const (
	defaultConfigRoot = "static"
	defaultMapRoot    = "maps"
	defaultResults    = "results"
	defaultCpuProfile = ""
)

// ParseCommandLine parses command-line arguments and returns a Config struct
func ParseCommandLine() *Config {
	// Define command-line parameters
	configRootArg := flag.String("config-root",
		defaultConfigRoot,
		"path to directory containing configuration files")
	mapRootArg := flag.String("map-root",
		defaultMapRoot,
		"path to directory containing map files")
	resultsArg := flag.String("results",
		defaultResults,
		"path to directory containing results files")
	cpuProfileArg := flag.String("cpu-profile",
		defaultCpuProfile,
		"path to CPU profile file (empty will disable CPU profiling)")
	webModeArg := flag.Bool("web",
		false,
		"run in web server mode")
	outputDirArg := flag.String("output-dir",
		"",
		"output directory for results (used by web GUI)")
	baseModelArg := flag.String("base-model",
		"",
		"base model to use (e.g., standard, flood, out_of_africa)")
	usernameArg := flag.String("username",
		"",
		"username for loading user models")
	modelNameArg := flag.String("model",
		"",
		"model name to load from user directory")
	checkpointOutArg := flag.String("checkpoint-out",
		"",
		"path to write run-state checkpoints (enables checkpointing)")
	checkpointIntervalArg := flag.Int("checkpoint-interval",
		0,
		"years between checkpoints (0 = only at end of each run)")
	checkpointInArg := flag.String("checkpoint-in",
		"",
		"path to a checkpoint to resume from")
	forkSeedArg := flag.Int("fork-seed",
		0,
		"when resuming, re-seed the RNG with this value for a divergent fork (0 = exact resume)")
	saveStateArg := flag.String("save-state",
		"",
		"name to save the end-of-run state under (in the model's saves/ directory)")
	loadStateArg := flag.String("load-state",
		"",
		"name of a saved state to start the run from (resume, or fork with -fork-seed)")
	listSavesArg := flag.Bool("list-saves",
		false,
		"list the model's named saved states and exit")
	validateArg := flag.Bool("validate",
		false,
		"run the neutral-expectation validation harness (roadmap §6h) and exit")
	importVCFArg := flag.String("import-vcf",
		"",
		"import a phased VCF (.vcf/.vcf.gz) into DRIFT structures, run the model's enabled analyses statically, and exit (roadmap §6a)")
	importPanelArg := flag.String("import-panel",
		"",
		"optional sample→population panel file (e.g. the 1000G .panel) for deme assignment on -import-vcf")
	importPolarizeArg := flag.String("import-polarize",
		"aa",
		"ancestral/derived polarization for -import-vcf: aa (AA INFO, REF fallback) or ref (always REF-ancestral)")

	argweaverTreesArg := flag.String("argweaver-trees",
		"",
		"parse `arg-summarize --tree` output into an inferred TMR-K-A and exit (TMR4A.md W7)")
	argweaverTruthArg := flag.String("argweaver-truth",
		"",
		"DRIFT tmrka_loci CSV to join the inferred TMR-K-A against (used with -argweaver-trees)")
	argweaverOutArg := flag.String("argweaver-out",
		"",
		"path to write the truth-vs-inference CSV (default: alongside -argweaver-trees)")
	argweaverKArg := flag.Int("argweaver-k",
		4,
		"K for the inferred TMR-K-A; 4 reproduces the published statistic")
	argweaverBpPerBitArg := flag.Int("argweaver-bp-per-bit",
		100,
		"bp-per-genome-bit convention the .sites file was exported under; must match argweaver_bp_per_bit")
	argweaverGenTimeArg := flag.Float64("argweaver-gen-time",
		25,
		"years per generation, for reporting inferred depths in years")
	argweaverMinSampleArg := flag.Int("argweaver-min-sample",
		0,
		"discard MCMC samples below this iteration (burn-in); establish the cut-off from the run's .log convergence check, not from this flag (0 keeps everything)")
	argweaverMaxSampleArg := flag.Int("argweaver-max-sample",
		0,
		"discard MCMC samples above this iteration; pair with -argweaver-min-sample for the early/late split that shows the answer has stopped moving (0 = no upper bound)")

	argweaverSMCArg := flag.String("argweaver-smc",
		"",
		"scan raw ARGweaver .smc(.gz) files for the genome-wide TMR-K-A distribution and its\n"+
			"\tspread across MCMC samples, then exit; comma-separated files and/or directories")
	argweaverSMCOutArg := flag.String("argweaver-smc-out",
		"",
		"path to write the per-MCMC-sample CSV (default: argweaver_smc_tmr<K>a.csv beside the first input)")
	argweaverSMCIterPathArg := flag.Bool("argweaver-smc-iter-from-path",
		false,
		"take the MCMC iteration from the whole path, not just the filename, for layouts\n"+
			"\tthat put the iteration in a directory (.../2400/chr1.smc.gz)")

	// Parse the command-line arguments
	clusterSitesArg := flag.String("cluster-sites",
		"", "Exported .sites file to cluster by pairwise divergence (TMR4A.md W9c stage 1)")
	clusterOutArg := flag.String("cluster-out",
		"", "Output prefix for the clustering CSVs (default: alongside the .sites file)")
	clusterWindowArg := flag.Int("cluster-window",
		2000000, "Cluster only sites at or below this bp position; 0 = the whole file")
	clusterNArg := flag.Int("cluster-n",
		20, "Individuals per bootstrap draw (20 = the study's vcf_sample_size)")
	clusterDrawsArg := flag.Int("cluster-draws",
		1000, "Number of bootstrap draws")
	clusterSeedArg := flag.Int64("cluster-seed",
		4242, "Seed for the resampling only; the realization is already fixed")
	clusterThresholdArg := flag.Float64("cluster-threshold",
		-1, "Single-linkage threshold; negative selects the parameter-free largest-gap cut")
	clusterMetricArg := flag.String("cluster-metric",
		"sites", "Distance unit: \"sites\" = differing-site count, comparable across sample "+
			"sizes; \"fraction\" = fraction of window sites, the legacy metric")
	clusterMinHeightArg := flag.Float64("cluster-min-height",
		-1, "Smallest merge height the largest-gap rule may cut above; negative = per-metric default")

	flag.Parse()

	// If output-dir is specified, use it for results
	resultsDir := *resultsArg
	if *outputDirArg != "" {
		resultsDir = *outputDirArg
	}

	// Return the configuration
	return &Config{
		ConfigRoot: *configRootArg,
		MapRoot:    *mapRootArg,
		Results:    resultsDir,
		CpuProfile: *cpuProfileArg,
		WebMode:    *webModeArg,
		OutputDir:  *outputDirArg,
		BaseModel:  *baseModelArg,
		Username:   *usernameArg,
		ModelName:  *modelNameArg,

		CheckpointOut:      *checkpointOutArg,
		CheckpointInterval: *checkpointIntervalArg,
		CheckpointIn:       *checkpointInArg,
		ForkSeed:           *forkSeedArg,

		SaveState: *saveStateArg,
		LoadState: *loadStateArg,
		ListSaves: *listSavesArg,

		Validate: *validateArg,

		ImportVCF:      *importVCFArg,
		ImportPanel:    *importPanelArg,
		ImportPolarize: *importPolarizeArg,

		ARGweaverTrees:     *argweaverTreesArg,
		ARGweaverTruth:     *argweaverTruthArg,
		ARGweaverOut:       *argweaverOutArg,
		ARGweaverK:         *argweaverKArg,
		ARGweaverBpPerBit:  *argweaverBpPerBitArg,
		ARGweaverGenTime:   *argweaverGenTimeArg,
		ARGweaverMinSample: *argweaverMinSampleArg,
		ARGweaverMaxSample: *argweaverMaxSampleArg,

		ARGweaverSMC:         *argweaverSMCArg,
		ARGweaverSMCOut:      *argweaverSMCOutArg,
		ARGweaverSMCIterPath: *argweaverSMCIterPathArg,

		ClusterSites:     *clusterSitesArg,
		ClusterOut:       *clusterOutArg,
		ClusterWindow:    *clusterWindowArg,
		ClusterN:         *clusterNArg,
		ClusterDraws:     *clusterDrawsArg,
		ClusterSeed:      *clusterSeedArg,
		ClusterThreshold: *clusterThresholdArg,
		ClusterMetric:    *clusterMetricArg,
		ClusterMinHeight: *clusterMinHeightArg,
	}
}
