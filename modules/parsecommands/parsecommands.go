package parsecommands

import (
	"flag"
)

// Config holds all the command-line configuration options
type Config struct {
	ConfigRoot string
	MapRoot    string
	Results    string
	CpuProfile string
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

	// Parse the command-line arguments
	flag.Parse()

	// Return the configuration
	return &Config{
		ConfigRoot: *configRootArg,
		MapRoot:    *mapRootArg,
		Results:    *resultsArg,
		CpuProfile: *cpuProfileArg,
	}
}
