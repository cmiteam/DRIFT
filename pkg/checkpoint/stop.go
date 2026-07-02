package checkpoint

import (
	"os"
	"path/filepath"
	"strings"
)

// stopFile is a sentinel dropped into a run's output directory to request a
// graceful stop-and-save. It is checked by the running simulation each year and
// removed once handled. A file-based signal is used because cross-process
// signals are unreliable on Windows.
const stopFile = ".stop"

// RequestStop writes a graceful-stop request into a run's output directory. The
// running simulation will save its state under saveName (see the sim's year
// loop) and then finish the run cleanly. An empty saveName lets the sim choose a
// default. outputDir is the directory passed to the run as -output-dir.
func RequestStop(outputDir, saveName string) error {
	return os.WriteFile(filepath.Join(outputDir, stopFile), []byte(saveName), 0644)
}

// StopRequested reports whether a graceful stop was requested for the run whose
// output directory is dir, and if so the save name to use (may be empty).
func StopRequested(dir string) (saveName string, ok bool) {
	data, err := os.ReadFile(filepath.Join(dir, stopFile))
	if err != nil {
		return "", false
	}
	return strings.TrimSpace(string(data)), true
}

// ClearStop removes a stop request (best effort).
func ClearStop(dir string) {
	_ = os.Remove(filepath.Join(dir, stopFile))
}
