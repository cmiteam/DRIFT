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

	// Parse the command-line arguments
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
	}
}
