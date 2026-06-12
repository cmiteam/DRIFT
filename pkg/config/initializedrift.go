package config

import (
	"drift/pkg/core"
	"fmt"
	"os"
	"path/filepath"
	"runtime/pprof"
)

// SetupDirectories creates necessary directories and handles setup tasks
func SetupDirectories(resultsPath string) error {
	// Create the results directory if it doesn't exist
	if _, err := os.Stat(resultsPath); os.IsNotExist(err) {
		err = os.MkdirAll(resultsPath, 0755)
		if err != nil {
			return fmt.Errorf("error creating results directory: %v", err)
		}
	}
	return nil
}

// SetupUserDirectories creates user-specific directory structure
// Structure: users/{username}/{model_name}/
func SetupUserDirectories(model *core.Model) error {
	if model.Username == "" {
		return nil // No username specified, skip user directory setup
	}

	// Create users directory
	usersDir := "users"
	if _, err := os.Stat(usersDir); os.IsNotExist(err) {
		if err := os.Mkdir(usersDir, 0755); err != nil {
			return fmt.Errorf("error creating users directory: %v", err)
		}
	}

	// Create user directory
	userDir := filepath.Join(usersDir, model.Username)
	if _, err := os.Stat(userDir); os.IsNotExist(err) {
		if err := os.Mkdir(userDir, 0755); err != nil {
			return fmt.Errorf("error creating user directory: %v", err)
		}
	}

	// Create model directory
	modelDir := filepath.Join(userDir, model.ModelName)
	if _, err := os.Stat(modelDir); os.IsNotExist(err) {
		if err := os.MkdirAll(modelDir, 0755); err != nil {
			return fmt.Errorf("error creating model directory: %v", err)
		}
	}

	// Update model's ResultsDir to point to user's model directory
	model.ResultsDir = modelDir

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
