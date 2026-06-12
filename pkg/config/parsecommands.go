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
	}
}
