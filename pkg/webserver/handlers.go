package webserver

import (
	"encoding/json"
	"fmt"
	"net/http"
	"os"
	"path/filepath"
	"strings"
)

// LoginRequest represents a login request
type LoginRequest struct {
	Username string `json:"username"`
	Password string `json:"password"`
}

// LoginResponse represents a login response
type LoginResponse struct {
	Success  bool   `json:"success"`
	Username string `json:"username,omitempty"`
	Message  string `json:"message,omitempty"`
}

// handleLogin handles user login
func handleLogin(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		respondError(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	var req LoginRequest
	err := json.NewDecoder(r.Body).Decode(&req)
	if err != nil {
		respondError(w, "Invalid request body", http.StatusBadRequest)
		return
	}

	// Validate user
	valid := usersDB.ValidateUser(req.Username, req.Password)
	if !valid {
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(LoginResponse{
			Success: false,
			Message: "Invalid username or password",
		})
		return
	}

	// Success
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(LoginResponse{
		Success:  true,
		Username: req.Username,
	})
}

// StartSimulationRequest represents a request to start a simulation
type StartSimulationRequest struct {
	Username  string `json:"username"`
	ModelName string `json:"model_name"`
}

// handleSimulationStart handles starting a new simulation
func handleSimulationStart(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		respondError(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	var req StartSimulationRequest
	err := json.NewDecoder(r.Body).Decode(&req)
	if err != nil {
		respondError(w, "Invalid request body", http.StatusBadRequest)
		return
	}

	// Validate user exists
	_, err = usersDB.GetUser(req.Username)
	if err != nil {
		respondError(w, "User not found", http.StatusNotFound)
		return
	}

	// Add job to queue
	job, err := queue.AddJob(req.Username, req.ModelName)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to queue simulation: %v", err), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(job)
}

// handleSimulationStatus returns status of user's simulations
func handleSimulationStatus(w http.ResponseWriter, r *http.Request) {
	username := r.URL.Query().Get("username")
	if username == "" {
		respondError(w, "Username parameter required", http.StatusBadRequest)
		return
	}

	jobs := queue.GetUserJobs(username)

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"jobs": jobs,
	})
}

// handleSimulationCancel cancels a simulation
func handleSimulationCancel(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		respondError(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	var req struct {
		JobID string `json:"job_id"`
	}
	err := json.NewDecoder(r.Body).Decode(&req)
	if err != nil {
		respondError(w, "Invalid request body", http.StatusBadRequest)
		return
	}

	err = queue.CancelJob(req.JobID)
	if err != nil {
		respondError(w, err.Error(), http.StatusBadRequest)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"success": true,
		"message": "Job cancelled",
	})
}

// handleProgress returns the progress.json for a specific simulation
func handleProgress(w http.ResponseWriter, r *http.Request) {
	// URL format: /api/progress/{username}/{modelname}
	parts := strings.Split(strings.TrimPrefix(r.URL.Path, "/api/progress/"), "/")
	if len(parts) < 2 {
		respondError(w, "Invalid URL format. Expected /api/progress/{username}/{modelname}", http.StatusBadRequest)
		return
	}

	username := parts[0]
	modelName := parts[1]

	progressPath := filepath.Join("users", username, modelName, "progress.json")

	data, err := os.ReadFile(progressPath)
	if err != nil {
		if os.IsNotExist(err) {
			respondError(w, "Progress file not found", http.StatusNotFound)
			return
		}
		respondError(w, fmt.Sprintf("Failed to read progress: %v", err), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	w.Write(data)
}

// handleResults serves result files for a user's simulation
func handleResults(w http.ResponseWriter, r *http.Request) {
	// URL format: /api/results/{username}/{modelname}/{filename}
	parts := strings.Split(strings.TrimPrefix(r.URL.Path, "/api/results/"), "/")
	if len(parts) < 2 {
		respondError(w, "Invalid URL format", http.StatusBadRequest)
		return
	}

	username := parts[0]
	modelName := parts[1]

	resultsDir := filepath.Join("users", username, modelName)

	// If no filename specified, list files
	if len(parts) == 2 {
		files, err := os.ReadDir(resultsDir)
		if err != nil {
			respondError(w, "Results directory not found", http.StatusNotFound)
			return
		}

		fileList := make([]map[string]interface{}, 0)
		for _, file := range files {
			if !file.IsDir() {
				info, _ := file.Info()
				fileList = append(fileList, map[string]interface{}{
					"name":    file.Name(),
					"size":    info.Size(),
					"modtime": info.ModTime(),
				})
			}
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"files": fileList,
		})
		return
	}

	// Serve specific file
	filename := parts[2]
	filePath := filepath.Join(resultsDir, filename)

	// Security check: ensure file is within results directory
	absFilePath, _ := filepath.Abs(filePath)
	absResultsDir, _ := filepath.Abs(resultsDir)
	if !strings.HasPrefix(absFilePath, absResultsDir) {
		respondError(w, "Invalid file path", http.StatusForbidden)
		return
	}

	// Check if file exists
	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		respondError(w, "File not found", http.StatusNotFound)
		return
	}

	// Serve the file
	http.ServeFile(w, r, filePath)
}

// handleGetUsers returns list of all users (for login dropdown)
func handleGetUsers(w http.ResponseWriter, r *http.Request) {
	usernames := make([]string, 0)
	for _, user := range usersDB.Users {
		usernames = append(usernames, user.Username)
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"users": usernames,
	})
}

// Admin handlers

// handleAdminQueue returns the full queue status (admin only)
func handleAdminQueue(w http.ResponseWriter, r *http.Request) {
	status := queue.GetQueueStatus()

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(status)
}

// handleAdminUsers returns all user information (admin only)
func handleAdminUsers(w http.ResponseWriter, r *http.Request) {
	if r.Method == http.MethodGet {
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(usersDB)
		return
	}

	if r.Method == http.MethodPost {
		// Add new user
		var req struct {
			Username string `json:"username"`
			Password string `json:"password"`
		}
		err := json.NewDecoder(r.Body).Decode(&req)
		if err != nil {
			respondError(w, "Invalid request body", http.StatusBadRequest)
			return
		}

		err = usersDB.AddUser(req.Username, req.Password)
		if err != nil {
			respondError(w, err.Error(), http.StatusBadRequest)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"success": true,
			"message": "User created",
		})
		return
	}

	respondError(w, "Method not allowed", http.StatusMethodNotAllowed)
}
