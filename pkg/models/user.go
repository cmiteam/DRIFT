// File: pkg/models/user.go
package models

import (
	"fmt"
	"os"
	"path/filepath"
	"time"
)

// InitializeUserModels creates a new user's models directory and copies all base models
func (m *ModelManager) InitializeUserModels(username string) error {
	userModelsPath := filepath.Join(m.UserModelsPath, username, "models")

	// Check if already initialized
	if _, err := os.Stat(userModelsPath); err == nil {
		return fmt.Errorf("user %s already has models initialized", username)
	}

	// Create models directory
	err := os.MkdirAll(userModelsPath, 0755)
	if err != nil {
		return fmt.Errorf("failed to create user models directory: %w", err)
	}

	// Get all base models
	baseModels, err := m.ListBaseModels()
	if err != nil {
		return fmt.Errorf("failed to list base models: %w", err)
	}

	// Copy each base model to user's models directory
	for _, baseModel := range baseModels {
		err = m.CreateModel(username, baseModel.ID, baseModel.ID, baseModel.Description)
		if err != nil {
			return fmt.Errorf("failed to create model %s for user: %w", baseModel.ID, err)
		}
	}

	return nil
}

// EnsureUserModelsExist checks if user has models, and initializes them if not
func (m *ModelManager) EnsureUserModelsExist(username string) error {
	userModelsPath := filepath.Join(m.UserModelsPath, username, "models")

	// Check if models directory exists
	if _, err := os.Stat(userModelsPath); os.IsNotExist(err) {
		// Initialize user models
		return m.InitializeUserModels(username)
	}

	return nil
}

// User represents a user in the system
type User struct {
	Username    string    `json:"username"`
	DisplayName string    `json:"display_name"`
	Created     time.Time `json:"created"`
}

// GetOrCreateUser ensures a user exists and has models initialized
func (m *ModelManager) GetOrCreateUser(username string) error {
	userPath := filepath.Join(m.UserModelsPath, username)

	// Create user directory if it doesn't exist
	err := os.MkdirAll(userPath, 0755)
	if err != nil {
		return fmt.Errorf("failed to create user directory: %w", err)
	}

	// Ensure models are initialized
	return m.EnsureUserModelsExist(username)
}
