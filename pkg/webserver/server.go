package webserver

import (
	"encoding/json"
	"log"
	"net/http"
	"os"
	"path/filepath"
)

// ServerConfig holds web server configuration
type ServerConfig struct {
	Port         string
	StaticDir    string
	DriftExePath string
}

// Start initializes and starts the web server
func Start() error {
	config := &ServerConfig{
		Port:         getEnv("PORT", "8080"),
		StaticDir:    "web/static",
		DriftExePath: getDriftExecutablePath(),
	}

	// Initialize simulation queue
	initQueue()

	// Setup routes
	mux := http.NewServeMux()

	// Serve static files (HTML, CSS, JS). no-cache forces the browser to REVALIDATE
	// every time rather than serving app.js from memory cache — it may still get a 304,
	// so this costs nothing, but an edit to app.js now shows up on an ordinary refresh.
	// Without it a stale app.js survives reloads and silently hides newly-added form
	// fields, which reads as "the parameters disappeared".
	fs := http.FileServer(http.Dir(config.StaticDir))
	mux.Handle("/", noCache(fs))

	// API routes
	mux.HandleFunc("/api/login", handleLogin)
	mux.HandleFunc("/api/register", handleRegister)
	mux.HandleFunc("/api/simulation/start", handleSimulationStart)
	mux.HandleFunc("/api/simulation/status", handleSimulationStatus)
	mux.HandleFunc("/api/simulation/cancel", handleSimulationCancel)
	mux.HandleFunc("/api/simulation/stop-save", handleSimulationStopSave)
	mux.HandleFunc("/api/progress/", handleProgress)
	mux.HandleFunc("/api/results/", handleResults)
	mux.HandleFunc("/api/plot-data/", handlePlotData)
	mux.HandleFunc("/api/users", handleGetUsers)
	mux.HandleFunc("/api/basemodels/list", handleBaseModelsList)
	mux.HandleFunc("/api/basemodels/parameters", handleBaseModelParameters)
	mux.HandleFunc("/api/modules/list", handlePhaseModules)
	mux.HandleFunc("/api/saves/list", handleSavesList)
	mux.HandleFunc("/api/saves/delete", handleSavesDelete)
	mux.HandleFunc("/api/maps/list", handleMapsList)

	// Model management routes
	mux.HandleFunc("/api/models/list", handleUserModelsList)
	mux.HandleFunc("/api/models/load", handleModelLoad)
	mux.HandleFunc("/api/models/save", handleModelSave)
	mux.HandleFunc("/api/models/reset", handleModelReset)

	// Admin routes
	mux.HandleFunc("/api/admin/queue", handleAdminQueue)
	mux.HandleFunc("/api/admin/users", handleAdminUsers)

	addr := ":" + config.Port
	log.Printf("Starting Drift web server on http://localhost%s", addr)
	log.Printf("Drift executable (config): %s", config.DriftExePath)
	log.Printf("Drift executable (detected): %s", getDriftExecutablePath())
	log.Printf("Current working directory: %s", func() string { wd, _ := os.Getwd(); return wd }())
	log.Printf("Serving static files from: %s", config.StaticDir)

	return http.ListenAndServe(addr, mux)
}

// getEnv gets environment variable or returns default value
func getEnv(key, defaultValue string) string {
	if value := os.Getenv(key); value != "" {
		return value
	}
	return defaultValue
}

// getDriftExecutablePath returns the path to the drift executable
// During development, returns "go" to use with "go run drift.go"
func getDriftExecutablePath() string {
	// Get current working directory
	wd, err := os.Getwd()
	if err != nil {
		wd = "."
	}

	// For development: use "go run drift.go" to always use latest code
	// Check if drift.go exists in current directory
	driftGoPath := filepath.Join(wd, "drift.go")
	if _, err := os.Stat(driftGoPath); err == nil {
		return "go"
	}

	// Production mode: use compiled executable
	// Check if drift.exe exists in current directory
	driftExePath := filepath.Join(wd, "drift.exe")
	if _, err := os.Stat(driftExePath); err == nil {
		absPath, _ := filepath.Abs(driftExePath)
		return absPath
	}

	// Check for drift (Unix-style)
	driftPath := filepath.Join(wd, "drift")
	if _, err := os.Stat(driftPath); err == nil {
		absPath, _ := filepath.Abs(driftPath)
		return absPath
	}

	// Fallback: return "go" and hope drift.go is in the current directory
	// This is better than returning a non-existent drift.exe
	return "go"
}

// Response helper functions
func respondJSON(w http.ResponseWriter, data interface{}) {
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(data)
}

func respondError(w http.ResponseWriter, message string, statusCode int) {
	w.Header().Set("Content-Type", "application/json")
	w.WriteHeader(statusCode)
	json.NewEncoder(w).Encode(map[string]interface{}{
		"success": false,
		"error":   message,
	})
}

// noCache asks the browser to revalidate a static asset on every request. ETag and
// Last-Modified still apply, so an unchanged file costs one 304 rather than a transfer.
func noCache(h http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Cache-Control", "no-cache, must-revalidate")
		h.ServeHTTP(w, r)
	})
}
