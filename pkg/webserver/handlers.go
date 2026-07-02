package webserver

import (
	"drift/pkg/checkpoint"
	"drift/pkg/config"
	"drift/pkg/modules"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"log"
	"net/http"
	"os"
	"path/filepath"
	"strconv"
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

// RegisterRequest represents a registration request
type RegisterRequest struct {
	Username string `json:"username"`
	Password string `json:"password"`
}

// RegisterResponse represents a registration response
type RegisterResponse struct {
	Success          bool   `json:"success"`
	Message          string `json:"message,omitempty"`
	RequiresApproval bool   `json:"requires_approval"`
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

// handleRegister handles user registration
func handleRegister(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		respondError(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	var req RegisterRequest
	err := json.NewDecoder(r.Body).Decode(&req)
	if err != nil {
		respondError(w, "Invalid request body", http.StatusBadRequest)
		return
	}

	// Validate input
	if req.Username == "" {
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(RegisterResponse{
			Success: false,
			Message: "Username is required",
		})
		return
	}

	// TODO: Implement proper password hashing (bcrypt)
	// For now, storing password as-is (NOT SECURE - placeholder only)
	// Empty password is allowed
	passwordHash := req.Password

	// Add user to database
	err = usersDB.AddUser(req.Username, passwordHash)
	if err != nil {
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(RegisterResponse{
			Success: false,
			Message: err.Error(),
		})
		return
	}

	// Initialize user models (copy all base models to user directory)
	err = config.GetOrCreateUser(req.Username)
	if err != nil {
		log.Printf("Warning: Failed to initialize models for user %s: %v", req.Username, err)
		// Don't fail registration - models can be initialized later
	} else {
		log.Printf("Initialized models for new user: %s", req.Username)
	}

	// Success - currently no approval required
	// TODO: Add admin approval system here when needed
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(RegisterResponse{
		Success:          true,
		RequiresApproval: false, // Placeholder: Set to true when admin approval is implemented
	})
}

// StartSimulationRequest represents a request to start a simulation
type StartSimulationRequest struct {
	Username       string  `json:"username"`
	BaseModel      string  `json:"base_model,omitempty"`
	ModelName      string  `json:"model_name"`
	EndYear        int     `json:"end_year,omitempty"`
	StartPopSize   int     `json:"start_pop_size,omitempty"`
	MaxPopSize     int     `json:"max_pop_size,omitempty"`
	SaveInterval   int     `json:"save_interval,omitempty"`
	NumRuns        int     `json:"num_runs,omitempty"`
	TrackDNA       int     `json:"track_DNA,omitempty"`
	TrackMap       int     `json:"track_map,omitempty"`
	TrackDrift     int     `json:"track_drift,omitempty"`
	TrackMutations int     `json:"track_mutations,omitempty"`
	MapName        string  `json:"map_name,omitempty"`
	Scaling        float64 `json:"scaling,omitempty"`
	SaveState      string  `json:"save_state,omitempty"`
	LoadState      string  `json:"load_state,omitempty"`
	ForkSeed       int     `json:"fork_seed,omitempty"`
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

	// Validate user exists and ensure models are initialized
	_, err = usersDB.GetUser(req.Username)
	if err != nil {
		respondError(w, "User not found", http.StatusNotFound)
		return
	}

	// Ensure user has models initialized
	err = config.GetOrCreateUser(req.Username)
	if err != nil {
		log.Printf("Warning: Could not initialize user models: %v", err)
	}

	// Add job to queue with parameters
	job, err := queue.AddJob(req.Username, req.ModelName, req)
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

// handleSimulationStopSave requests a graceful stop of the running simulation:
// it saves its current state under a name and finishes the run cleanly.
func handleSimulationStopSave(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		respondError(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}
	var req struct {
		JobID    string `json:"job_id"`
		SaveName string `json:"save_name"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		respondError(w, "Invalid request body", http.StatusBadRequest)
		return
	}
	name := strings.TrimSpace(req.SaveName)
	if name == "" {
		name = "stopped"
	}
	if err := queue.StopSaveJob(req.JobID, name); err != nil {
		respondError(w, err.Error(), http.StatusBadRequest)
		return
	}
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"success": true,
		"message": fmt.Sprintf("Stopping and saving as %q; the run will finish the current year and exit.", name),
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

	progressPath := filepath.Join("users", username, "models", modelName, "results", "progress.json")

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

// handlePlotData serves incremental plot data from CSV results
func handlePlotData(w http.ResponseWriter, r *http.Request) {
	// URL format: /api/plot-data/{username}/{modelname}?from_row=N
	parts := strings.Split(strings.TrimPrefix(r.URL.Path, "/api/plot-data/"), "/")
	if len(parts) < 2 {
		respondError(w, "Invalid URL format. Expected /api/plot-data/{username}/{modelname}", http.StatusBadRequest)
		return
	}

	username := parts[0]
	modelName := parts[1]

	// Build path to results CSV
	csvPath := filepath.Join("users", username, "models", modelName, "results", modelName+"_results.csv")

	// Check if file exists
	fileInfo, err := os.Stat(csvPath)
	if err != nil {
		if os.IsNotExist(err) {
			respondError(w, "Results file not found", http.StatusNotFound)
			return
		}
		respondError(w, fmt.Sprintf("Failed to access results: %v", err), http.StatusInternalServerError)
		return
	}

	// Check If-Modified-Since header for caching
	ifModifiedSince := r.Header.Get("If-Modified-Since")
	if ifModifiedSince != "" {
		// Parse the time
		t, err := http.ParseTime(ifModifiedSince)
		if err == nil && !fileInfo.ModTime().After(t) {
			w.WriteHeader(http.StatusNotModified)
			return
		}
	}

	// Get from_row parameter (default 0)
	fromRow := 0
	if fromRowStr := r.URL.Query().Get("from_row"); fromRowStr != "" {
		fromRow, _ = strconv.Atoi(fromRowStr)
	}

	// Read CSV file
	file, err := os.Open(csvPath)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to open results file: %v", err), http.StatusInternalServerError)
		return
	}
	defer file.Close()

	reader := csv.NewReader(file)

	// Read header
	headers, err := reader.Read()
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to read CSV header: %v", err), http.StatusInternalServerError)
		return
	}

	// Read all rows
	allRows, err := reader.ReadAll()
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to read CSV data: %v", err), http.StatusInternalServerError)
		return
	}

	// Get rows starting from fromRow
	var rows [][]string
	if fromRow < len(allRows) {
		rows = allRows[fromRow:]
	} else {
		rows = [][]string{}
	}

	// Convert to JSON-friendly format
	data := make([]map[string]string, 0, len(rows))
	for _, row := range rows {
		rowData := make(map[string]string)
		for i, header := range headers {
			if i < len(row) {
				rowData[header] = row[i]
			}
		}
		data = append(data, rowData)
	}

	// Send response
	w.Header().Set("Content-Type", "application/json")
	w.Header().Set("Last-Modified", fileInfo.ModTime().UTC().Format(http.TimeFormat))
	json.NewEncoder(w).Encode(map[string]interface{}{
		"headers":    headers,
		"rows":       data,
		"total_rows": len(allRows),
		"from_row":   fromRow,
		"new_rows":   len(rows),
		"modified":   fileInfo.ModTime().UTC().Format(http.TimeFormat),
	})
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

	resultsDir := filepath.Join("users", username, "models", modelName, "results")

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

// handleBaseModelsList returns list of available base models
func handleBaseModelsList(w http.ResponseWriter, r *http.Request) {
	models, err := config.ListBaseModels()
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to load base models: %v", err), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"models": models,
	})
}

// handleSavesList returns the named saved run-states for a model, with metadata.
func handleSavesList(w http.ResponseWriter, r *http.Request) {
	username := r.URL.Query().Get("username")
	modelName := r.URL.Query().Get("model_name")
	if username == "" || modelName == "" {
		respondError(w, "username and model_name required", http.StatusBadRequest)
		return
	}
	saves, err := checkpoint.List(checkpoint.ModelFor(username, modelName))
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to list saves: %v", err), http.StatusInternalServerError)
		return
	}
	if saves == nil {
		saves = []checkpoint.SaveInfo{}
	}
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"success": true,
		"saves":   saves,
	})
}

// handleSavesDelete removes a named saved run-state.
func handleSavesDelete(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		respondError(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}
	var req struct {
		Username  string `json:"username"`
		ModelName string `json:"model_name"`
		Name      string `json:"name"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		respondError(w, "Invalid request body", http.StatusBadRequest)
		return
	}
	if req.Username == "" || req.ModelName == "" || req.Name == "" {
		respondError(w, "username, model_name and name required", http.StatusBadRequest)
		return
	}
	if err := checkpoint.DeleteNamed(checkpoint.ModelFor(req.Username, req.ModelName), req.Name); err != nil {
		respondError(w, fmt.Sprintf("Failed to delete save: %v", err), http.StatusInternalServerError)
		return
	}
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{"success": true})
}

// handlePhaseModules returns the available module names per phase parameter
// (birth_style, death_style, mating_style, seed_style, setup_style) so the GUI
// can render a dropdown for each. Populated from the compiled-in registry.
func handlePhaseModules(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"styles": modules.AllStyles(),
	})
}

// handleBaseModelParameters returns parameters for a specific base model
func handleBaseModelParameters(w http.ResponseWriter, r *http.Request) {
	modelID := r.URL.Query().Get("id")
	if modelID == "" {
		respondError(w, "Missing model ID parameter", http.StatusBadRequest)
		return
	}

	// Load parameters from base model CSV file
	paramsPath := filepath.Join("static", "basemodels", modelID, "parameters.csv")
	file, err := os.Open(paramsPath)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to open base model parameters: %v", err), http.StatusNotFound)
		return
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to read parameters CSV: %v", err), http.StatusInternalServerError)
		return
	}

	// Convert CSV to map (skip header row)
	parameters := make(map[string]interface{})
	for i, record := range records {
		if i == 0 || len(record) < 2 {
			continue // Skip header
		}
		paramName := strings.TrimSpace(record[0])
		paramValue := strings.TrimSpace(record[1])

		// Try to convert to appropriate type
		if intVal, err := strconv.Atoi(paramValue); err == nil {
			parameters[paramName] = intVal
		} else if floatVal, err := strconv.ParseFloat(paramValue, 64); err == nil {
			parameters[paramName] = floatVal
		} else {
			parameters[paramName] = paramValue
		}
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"success":    true,
		"parameters": parameters,
		"base_model": modelID,
	})
}

// handleMapsList returns list of available maps
func handleMapsList(w http.ResponseWriter, r *http.Request) {
	mapsDir := "maps"

	// Check if maps directory exists
	if _, err := os.Stat(mapsDir); os.IsNotExist(err) {
		// Return default map if directory doesn't exist
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"maps": []string{"sandbox"},
		})
		return
	}

	// Read directory contents
	entries, err := os.ReadDir(mapsDir)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to read maps directory: %v", err), http.StatusInternalServerError)
		return
	}

	// Collect map names (directories and .csv files without extension)
	mapNames := make([]string, 0)
	seenMaps := make(map[string]bool)

	for _, entry := range entries {
		var mapName string

		if entry.IsDir() {
			mapName = entry.Name()
		} else if strings.HasSuffix(entry.Name(), ".csv") {
			mapName = strings.TrimSuffix(entry.Name(), ".csv")
		} else {
			continue
		}

		// Avoid duplicates
		if !seenMaps[mapName] {
			mapNames = append(mapNames, mapName)
			seenMaps[mapName] = true
		}
	}

	// If no maps found, return default
	if len(mapNames) == 0 {
		mapNames = []string{"sandbox"}
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"maps": mapNames,
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
