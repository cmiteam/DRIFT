package initializedrift

import (
	"fmt"
	"os"
	"runtime/pprof"
)

// SetupDirectories creates necessary directories and handles setup tasks
func SetupDirectories(resultsPath string) error {
	// Create the results directory if it doesn't exist
	if _, err := os.Stat(resultsPath); os.IsNotExist(err) {
		err = os.Mkdir(resultsPath, 0755)
		if err != nil {
			return fmt.Errorf("error creating results directory: %v", err)
		}
	}
	return nil
}

// StartCPUProfiling starts CPU profiling if a profile path is provided
func StartCPUProfiling(cpuProfilePath string) (func(), error) {
	if cpuProfilePath == "" {
		// Return a no-op cleanup function if profiling is disabled
		return func() {}, nil
	}

	cpuFile, err := os.Create(cpuProfilePath)
	if err != nil {
		return nil, fmt.Errorf("error creating CPU profile file: %v", err)
	}

	err = pprof.StartCPUProfile(cpuFile)
	if err != nil {
		cpuFile.Close()
		return nil, fmt.Errorf("error starting CPU profile: %v", err)
	}

	// Return cleanup function that stops profiling and closes the file
	cleanup := func() {
		pprof.StopCPUProfile()
		cpuFile.Close()
	}

	return cleanup, nil
}
