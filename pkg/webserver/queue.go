package webserver

import (
	"fmt"
	"os/exec"
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
	ID          string           `json:"id"`
	Username    string           `json:"username"`
	ModelName   string           `json:"model_name"`
	Status      SimulationStatus `json:"status"`
	QueuedAt    time.Time        `json:"queued_at"`
	StartedAt   *time.Time       `json:"started_at,omitempty"`
	CompletedAt *time.Time       `json:"completed_at,omitempty"`
	Error       string           `json:"error,omitempty"`
	Position    int              `json:"position"`
	OutputDir   string           `json:"output_dir"`
	cmd         *exec.Cmd
}

// SimulationQueue manages the queue of simulations
type SimulationQueue struct {
	jobs    []*SimulationJob
	running *SimulationJob
	mu      sync.RWMutex
	jobChan chan *SimulationJob
}

var queue *SimulationQueue

// initQueue initializes the simulation queue
func initQueue() {
	queue = &SimulationQueue{
		jobs:    make([]*SimulationJob, 0),
		jobChan: make(chan *SimulationJob, 100),
	}

	// Start the queue processor
	go queue.processor()
}

// AddJob adds a new simulation job to the queue
func (q *SimulationQueue) AddJob(username, modelName string) (*SimulationJob, error) {
	q.mu.Lock()
	defer q.mu.Unlock()

	// Generate unique job ID
	jobID := fmt.Sprintf("%s_%s_%d", username, modelName, time.Now().Unix())

	// Create output directory path
	outputDir := fmt.Sprintf("users/%s/%s", username, modelName)

	job := &SimulationJob{
		ID:        jobID,
		Username:  username,
		ModelName: modelName,
		Status:    StatusQueued,
		QueuedAt:  time.Now(),
		Position:  len(q.jobs) + 1,
		OutputDir: outputDir,
	}

	q.jobs = append(q.jobs, job)
	q.jobChan <- job

	return job, nil
}

// GetJob returns a job by ID
func (q *SimulationQueue) GetJob(jobID string) (*SimulationJob, error) {
	q.mu.RLock()
	defer q.mu.RUnlock()

	for _, job := range q.jobs {
		if job.ID == jobID {
			return job, nil
		}
	}

	if q.running != nil && q.running.ID == jobID {
		return q.running, nil
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
		return nil
	}

	// Check queued jobs
	for i, job := range q.jobs {
		if job.ID == jobID {
			job.Status = StatusCancelled
			now := time.Now()
			job.CompletedAt = &now
			// Remove from queue
			q.jobs = append(q.jobs[:i], q.jobs[i+1:]...)
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
		q.running = nil
		q.mu.Unlock()
	}
}

// runSimulation executes the drift simulation
func (q *SimulationQueue) runSimulation(job *SimulationJob) error {
	// Get drift executable path from server config
	driftPath := getDriftExecutablePath()

	// Build command
	cmd := exec.Command(driftPath,
		"--output-dir", job.OutputDir,
		"--config-root", "static", // Use default for now
	)

	job.cmd = cmd

	// Run the command
	output, err := cmd.CombinedOutput()
	if err != nil {
		return fmt.Errorf("simulation failed: %v\nOutput: %s", err, string(output))
	}

	return nil
}

// GetQueueStatus returns current queue status
func (q *SimulationQueue) GetQueueStatus() map[string]interface{} {
	q.mu.RLock()
	defer q.mu.RUnlock()

	status := map[string]interface{}{
		"queue_length": len(q.jobs),
		"running":      q.running,
		"queued_jobs":  q.jobs,
	}

	return status
}
