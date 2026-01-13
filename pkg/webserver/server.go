package webserver

import (
	"fmt"
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

	// Serve static files (HTML, CSS, JS)
	fs := http.FileServer(http.Dir(config.StaticDir))
	mux.Handle("/", fs)

	// API routes
	mux.HandleFunc("/api/login", handleLogin)
	mux.HandleFunc("/api/simulation/start", handleSimulationStart)
	mux.HandleFunc("/api/simulation/status", handleSimulationStatus)
	mux.HandleFunc("/api/simulation/cancel", handleSimulationCancel)
	mux.HandleFunc("/api/progress/", handleProgress)
	mux.HandleFunc("/api/results/", handleResults)
	mux.HandleFunc("/api/users", handleGetUsers)

	// Admin routes
	mux.HandleFunc("/api/admin/queue", handleAdminQueue)
	mux.HandleFunc("/api/admin/users", handleAdminUsers)

	addr := ":" + config.Port
	log.Printf("Starting Drift web server on http://localhost%s", addr)
	log.Printf("Drift executable: %s", config.DriftExePath)
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
func getDriftExecutablePath() string {
	// Check if drift.exe exists in current directory
	if _, err := os.Stat("drift.exe"); err == nil {
		absPath, _ := filepath.Abs("drift.exe")
		return absPath
	}

	// Check for drift (Unix-style)
	if _, err := os.Stat("drift"); err == nil {
		absPath, _ := filepath.Abs("drift")
		return absPath
	}

	// Fallback to assuming drift.exe in current directory
	absPath, _ := filepath.Abs("drift.exe")
	return absPath
}

// Response helper functions
func respondJSON(w http.ResponseWriter, data interface{}) {
	w.Header().Set("Content-Type", "application/json")
	fmt.Fprintf(w, "%v", data)
}

func respondError(w http.ResponseWriter, message string, statusCode int) {
	w.Header().Set("Content-Type", "application/json")
	w.WriteHeader(statusCode)
	fmt.Fprintf(w, `{"error": "%s"}`, message)
}
