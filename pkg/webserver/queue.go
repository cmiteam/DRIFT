package webserver

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"sync"
	"time"
)

// SimulationStatus represents the status of a simulation
type SimulationStatus string

const (
	StatusQueued    SimulationStatus = "queued"
	StatusRunning   SimulationStatus = "running"
	StatusCompleted SimulationStatus = "completed"
	StatusFailed    SimulationStatus = "failed"
	StatusCancelled SimulationStatus = "cancelled"
)

// SimulationJob represents a simulation in the queue
type SimulationJob struct {
	ID             string           `json:"id"`
	Username       string           `json:"username"`
	BaseModel      string           `json:"base_model,omitempty"`
	ModelName      string           `json:"model_name"`
	Status         SimulationStatus `json:"status"`
	QueuedAt       time.Time        `json:"queued_at"`
	StartedAt      *time.Time       `json:"started_at,omitempty"`
	CompletedAt    *time.Time       `json:"completed_at,omitempty"`
	Error          string           `json:"error,omitempty"`
	Position       int              `json:"position"`
	OutputDir      string           `json:"output_dir"`
	EndYear        int              `json:"end_year,omitempty"`
	StartPopSize   int              `json:"start_pop_size,omitempty"`
	MaxPopSize     int              `json:"max_pop_size,omitempty"`
	SaveInterval   int              `json:"save_interval,omitempty"`
	NumRuns        int              `json:"num_runs,omitempty"`
	TrackDNA       int              `json:"track_DNA,omitempty"`
	TrackMap       int              `json:"track_map,omitempty"`
	TrackDrift     int              `json:"track_drift,omitempty"`
	TrackMutations int              `json:"track_mutations,omitempty"`
	MapName        string           `json:"map_name,omitempty"`
	Scaling        float64          `json:"scaling,omitempty"`
	cmd            *exec.Cmd
}

// SimulationQueue manages the queue of simulations
type SimulationQueue struct {
	jobs      []*SimulationJob
	running   *SimulationJob
	completed []*SimulationJob
	mu        sync.RWMutex
	jobChan   chan *SimulationJob
}

var queue *SimulationQueue

// initQueue initializes the simulation queue
func initQueue() {
	queue = &SimulationQueue{
		jobs:      make([]*SimulationJob, 0),
		completed: make([]*SimulationJob, 0),
		jobChan:   make(chan *SimulationJob, 100),
	}

	// Start the queue processor
	go queue.processor()
}

// AddJob adds a new simulation job to the queue
func (q *SimulationQueue) AddJob(username, modelName string, params interface{}) (*SimulationJob, error) {
	q.mu.Lock()
	defer q.mu.Unlock()

	// Generate unique job ID
	jobID := fmt.Sprintf("%s_%s_%d", username, modelName, time.Now().Unix())

	// Create output directory path
	outputDir := fmt.Sprintf("users/%s/models/%s/results", username, modelName)

	job := &SimulationJob{
		ID:        jobID,
		Username:  username,
		ModelName: modelName,
		Status:    StatusQueued,
		QueuedAt:  time.Now(),
		Position:  len(q.jobs) + 1,
		OutputDir: outputDir,
	}

	// Copy parameters if provided
	if req, ok := params.(StartSimulationRequest); ok {
		job.BaseModel = req.BaseModel
		job.EndYear = req.EndYear
		job.StartPopSize = req.StartPopSize
		job.MaxPopSize = req.MaxPopSize
		job.SaveInterval = req.SaveInterval
		job.NumRuns = req.NumRuns
		job.TrackDNA = req.TrackDNA
		job.TrackMap = req.TrackMap
		job.TrackDrift = req.TrackDrift
		job.TrackMutations = req.TrackMutations
		job.MapName = req.MapName
		job.Scaling = req.Scaling
	}

	q.jobs = append(q.jobs, job)
	q.jobChan <- job

	return job, nil
}

// GetJob returns a job by ID
func (q *SimulationQueue) GetJob(jobID string) (*SimulationJob, error) {
	q.mu.RLock()
	defer q.mu.RUnlock()

	// Check running job
	if q.running != nil && q.running.ID == jobID {
		return q.running, nil
	}

	// Check queued jobs
	for _, job := range q.jobs {
		if job.ID == jobID {
			return job, nil
		}
	}

	// Check completed jobs
	for _, job := range q.completed {
		if job.ID == jobID {
			return job, nil
		}
	}

	return nil, fmt.Errorf("job not found: %s", jobID)
}

// GetUserJobs returns all jobs for a specific user
func (q *SimulationQueue) GetUserJobs(username string) []*SimulationJob {
	q.mu.RLock()
	defer q.mu.RUnlock()

	userJobs := make([]*SimulationJob, 0)

	// Check running job
	if q.running != nil && q.running.Username == username {
		userJobs = append(userJobs, q.running)
	}

	// Check queued jobs
	for _, job := range q.jobs {
		if job.Username == username {
			userJobs = append(userJobs, job)
		}
	}

	// Check completed jobs (most recent first)
	for i := len(q.completed) - 1; i >= 0; i-- {
		if q.completed[i].Username == username {
			userJobs = append(userJobs, q.completed[i])
		}
	}

	return userJobs
}

// CancelJob cancels a queued or running job
func (q *SimulationQueue) CancelJob(jobID string) error {
	q.mu.Lock()
	defer q.mu.Unlock()

	// Check if it's the running job
	if q.running != nil && q.running.ID == jobID {
		if q.running.cmd != nil && q.running.cmd.Process != nil {
			err := q.running.cmd.Process.Kill()
			if err != nil {
				return fmt.Errorf("failed to kill process: %v", err)
			}
		}
		q.running.Status = StatusCancelled
		now := time.Now()
		q.running.CompletedAt = &now
		// Move to completed list
		q.completed = append(q.completed, q.running)
		q.running = nil
		return nil
	}

	// Check queued jobs
	for i, job := range q.jobs {
		if job.ID == jobID {
			job.Status = StatusCancelled
			now := time.Now()
			job.CompletedAt = &now
			// Remove from queue and add to completed
			q.jobs = append(q.jobs[:i], q.jobs[i+1:]...)
			q.completed = append(q.completed, job)
			return nil
		}
	}

	return fmt.Errorf("job not found or already completed: %s", jobID)
}

// processor processes jobs from the queue
func (q *SimulationQueue) processor() {
	for job := range q.jobChan {
		if job.Status == StatusCancelled {
			continue
		}

		q.mu.Lock()
		q.running = job
		job.Status = StatusRunning
		now := time.Now()
		job.StartedAt = &now
		q.mu.Unlock()

		// Remove from queue
		q.mu.Lock()
		for i, j := range q.jobs {
			if j.ID == job.ID {
				q.jobs = append(q.jobs[:i], q.jobs[i+1:]...)
				break
			}
		}
		q.mu.Unlock()

		// Run the simulation
		err := q.runSimulation(job)

		q.mu.Lock()
		if err != nil {
			job.Status = StatusFailed
			job.Error = err.Error()
		} else {
			job.Status = StatusCompleted
		}
		completed := time.Now()
		job.CompletedAt = &completed
		// Move to completed list
		q.completed = append(q.completed, job)
		q.running = nil
		q.mu.Unlock()
	}
}

// runSimulation executes the drift simulation
func (q *SimulationQueue) runSimulation(job *SimulationJob) error {
	// Create output directory if it doesn't exist
	err := os.MkdirAll(job.OutputDir, 0755)
	if err != nil {
		return fmt.Errorf("failed to create output directory: %v", err)
	}

	// Update user model parameters directly
	err = q.updateUserModelParameters(job)
	if err != nil {
		return fmt.Errorf("failed to update model parameters: %v", err)
	}

	// Get drift executable path from server config
	driftPath := getDriftExecutablePath()
	log.Printf("Drift executable path: %s", driftPath)

	// Build command arguments - use the user model system
	args := []string{
		"-username", job.Username,
		"-model", job.ModelName,
		"-output-dir", job.OutputDir,
		"-map-root", "maps",
	}

	var cmd *exec.Cmd

	// If using "go" command for development, prepend "run drift.go"
	if driftPath == "go" {
		// Development mode: go run drift.go [args...]
		goArgs := []string{"run", "drift.go"}
		goArgs = append(goArgs, args...)
		cmd = exec.Command(driftPath, goArgs...)
	} else {
		// Production mode: drift.exe [args...]
		cmd = exec.Command(driftPath, args...)
	}

	// Set working directory to the project root where drift.go is located
	// This is important for "go run drift.go" to find the source file
	wd, _ := os.Getwd()
	cmd.Dir = wd

	job.cmd = cmd

	log.Printf("Running simulation: %s %v in directory: %s", driftPath, args, cmd.Dir)

	// Run the command
	output, err := cmd.CombinedOutput()
	if err != nil {
		errorMsg := fmt.Sprintf("simulation failed: %v\nOutput: %s", err, string(output))
		log.Printf("Simulation error for job %s: %s", job.ID, errorMsg)
		return fmt.Errorf(errorMsg)
	}

	log.Printf("Simulation completed successfully for job %s", job.ID)
	return nil
}

// updateUserModelParameters updates the user's model parameters.csv file directly
func (q *SimulationQueue) updateUserModelParameters(job *SimulationJob) error {
	// Path to user's model parameters file
	paramsPath := filepath.Join("users", job.Username, "models", job.ModelName, "config", "parameters.csv")

	log.Printf("Updating user model parameters at %s", paramsPath)

	// Read existing parameters
	file, err := os.Open(paramsPath)
	if err != nil {
		return fmt.Errorf("failed to open parameters file: %v", err)
	}

	reader := csv.NewReader(file)
	reader.FieldsPerRecord = -1 // tolerate inconsistent widths from legacy/corrupt files
	records, err := reader.ReadAll()
	file.Close()
	if err != nil {
		return fmt.Errorf("failed to read parameters file: %v", err)
	}

	// Detect schema from header: 6-col (parameter,label,entry,format,value,group) or 2-col (parameter,value)
	headerWidth := 2
	if len(records) > 0 && len(records[0]) >= 6 {
		headerWidth = 6
	}
	valueIdx := 1
	if headerWidth == 6 {
		valueIdx = 4
	}

	// Build a map of parameter updates from job
	updates := make(map[string]string)
	if job.EndYear > 0 {
		updates["end_year"] = strconv.Itoa(job.EndYear)
	}
	if job.StartPopSize > 0 {
		updates["start_pop_size"] = strconv.Itoa(job.StartPopSize)
	}
	if job.MaxPopSize > 0 {
		updates["max_pop_size"] = strconv.Itoa(job.MaxPopSize)
	}
	if job.SaveInterval > 0 {
		updates["save_interval"] = strconv.Itoa(job.SaveInterval)
	}
	if job.NumRuns > 0 {
		updates["num_runs"] = strconv.Itoa(job.NumRuns)
	}
	updates["track_DNA"] = strconv.Itoa(job.TrackDNA)
	updates["track_map"] = strconv.Itoa(job.TrackMap)
	updates["track_drift"] = strconv.Itoa(job.TrackDrift)
	updates["track_mutations"] = strconv.Itoa(job.TrackMutations)
	if job.MapName != "" {
		updates["map_name"] = job.MapName
	}
	if job.Scaling > 0 {
		updates["scaling"] = strconv.FormatFloat(job.Scaling, 'f', -1, 64)
	}

	// Track which parameters exist in file
	existingParams := make(map[string]bool)

	// Update existing records
	for i, record := range records {
		if i == 0 || len(record) < 2 {
			continue // Skip header or invalid rows
		}
		paramName := record[0]
		existingParams[paramName] = true
		if newValue, ok := updates[paramName]; ok {
			// Ensure the row is wide enough to hold the value at valueIdx (pad if a legacy 2-col row snuck into a 6-col file)
			for len(record) <= valueIdx {
				record = append(record, "")
			}
			record[valueIdx] = newValue
			records[i] = record
			log.Printf("Updated %s = %s", paramName, newValue)
		}
	}

	// Add any missing parameters with rows matching the file's schema width
	for param, value := range updates {
		if !existingParams[param] {
			var row []string
			if headerWidth == 6 {
				row = []string{param, param, "Text", "float", value, "Main"}
			} else {
				row = []string{param, value}
			}
			records = append(records, row)
			log.Printf("Added new parameter %s = %s", param, value)
		}
	}

	// Write back to file
	outFile, err := os.Create(paramsPath)
	if err != nil {
		return fmt.Errorf("failed to create parameters file: %v", err)
	}
	defer outFile.Close()

	writer := csv.NewWriter(outFile)
	err = writer.WriteAll(records)
	if err != nil {
		return fmt.Errorf("failed to write parameters file: %v", err)
	}

	return nil
}

// createCustomConfig creates a custom config directory with user's parameter overrides
func (q *SimulationQueue) createCustomConfig(job *SimulationJob, configRoot string) error {
	// Create config directory
	err := os.MkdirAll(configRoot, 0755)
	if err != nil {
		return err
	}

	// Copy all static config files to custom config directory
	staticFiles := []string{
		"actuarial_table.csv",
		"chromosome_data.csv",
		"Pisum_chromosome_data.csv",
		"plot_defaults.csv",
		"basemodels.csv",
	}

	// Also copy base model parameter files if they exist
	baseModelFiles := []string{
		"basemodel_standard.csv",
		"basemodel_flood.csv",
		"basemodel_out_of_africa.csv",
	}

	for _, file := range staticFiles {
		src := filepath.Join("static", file)
		dst := filepath.Join(configRoot, file)
		err := copyFile(src, dst)
		if err != nil {
			return fmt.Errorf("failed to copy %s: %v", file, err)
		}
	}

	// Copy base model files (ignore errors if they don't exist)
	for _, file := range baseModelFiles {
		src := filepath.Join("static", file)
		dst := filepath.Join(configRoot, file)
		copyFile(src, dst) // Ignore errors - file may not exist
	}

	// Create custom parameter_defaults.csv with user's overrides
	err = q.createParameterFile(job, configRoot)
	if err != nil {
		return fmt.Errorf("failed to create parameter file: %v", err)
	}

	return nil
}

// createParameterFile creates parameter_defaults.csv with user's custom values
func (q *SimulationQueue) createParameterFile(job *SimulationJob, configRoot string) error {
	log.Printf("Creating parameter file for job %s with SaveInterval=%d", job.ID, job.SaveInterval)

	// Read default parameters
	defaultPath := filepath.Join("static", "parameter_defaults.csv")
	file, err := os.Open(defaultPath)
	if err != nil {
		return err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		return err
	}

	// Update values based on job parameters
	for i, record := range records {
		if i == 0 {
			continue // Skip header
		}
		if len(record) < 5 {
			continue
		}

		param := record[0]
		switch param {
		case "end_year":
			if job.EndYear > 0 {
				record[4] = strconv.Itoa(job.EndYear)
			}
		case "start_pop_size":
			if job.StartPopSize > 0 {
				record[4] = strconv.Itoa(job.StartPopSize)
			}
		case "max_pop_size":
			if job.MaxPopSize > 0 {
				record[4] = strconv.Itoa(job.MaxPopSize)
			}
		case "save_interval":
			if job.SaveInterval > 0 {
				log.Printf("Updating save_interval from %s to %d", record[4], job.SaveInterval)
				record[4] = strconv.Itoa(job.SaveInterval)
			}
		case "num_runs":
			if job.NumRuns > 0 {
				record[4] = strconv.Itoa(job.NumRuns)
			}
		case "track_DNA":
			record[4] = strconv.Itoa(job.TrackDNA)
		case "track_map":
			record[4] = strconv.Itoa(job.TrackMap)
		case "map_name":
			if job.MapName != "" {
				record[4] = job.MapName
			}
		case "scaling":
			if job.Scaling > 0 {
				record[4] = strconv.FormatFloat(job.Scaling, 'f', -1, 64)
			}
		case "track_drift":
			record[4] = strconv.Itoa(job.TrackDrift)
		case "track_mutations":
			record[4] = strconv.Itoa(job.TrackMutations)
		case "username":
			record[4] = job.Username
		case "model_name":
			record[4] = job.ModelName
		case "base_model_id":
			if job.BaseModel != "" {
				record[4] = job.BaseModel
			}
		}
	}

	// Write updated parameters to custom config directory
	customPath := filepath.Join(configRoot, "parameter_defaults.csv")
	outFile, err := os.Create(customPath)
	if err != nil {
		return err
	}
	defer outFile.Close()

	writer := csv.NewWriter(outFile)
	err = writer.WriteAll(records)
	if err != nil {
		return err
	}

	return nil
}

// copyFile copies a file from src to dst
func copyFile(src, dst string) error {
	data, err := os.ReadFile(src)
	if err != nil {
		return err
	}
	return os.WriteFile(dst, data, 0644)
}

// GetQueueStatus returns current queue status
func (q *SimulationQueue) GetQueueStatus() map[string]interface{} {
	q.mu.RLock()
	defer q.mu.RUnlock()

	status := map[string]interface{}{
		"queue_length":    len(q.jobs),
		"running":         q.running,
		"queued_jobs":     q.jobs,
		"completed_jobs":  q.completed,
		"completed_count": len(q.completed),
	}

	return status
}
