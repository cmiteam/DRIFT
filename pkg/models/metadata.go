// File: pkg/models/metadata.go
package models

import (
	"encoding/json"
	"fmt"
	"os"
	"time"
)

// ModelMetadata contains information about a user model
type ModelMetadata struct {
	Name         string    `json:"name"`
	Description  string    `json:"description"`
	BaseModel    string    `json:"base_model"`     // "standard", "flood", etc.
	Created      time.Time `json:"created"`
	LastModified time.Time `json:"last_modified"`
	Author       string    `json:"author"`
}

// LoadMetadata loads metadata from a JSON file
func LoadMetadata(filepath string) (*ModelMetadata, error) {
	data, err := os.ReadFile(filepath)
	if err != nil {
		return nil, fmt.Errorf("failed to read metadata file: %w", err)
	}

	var metadata ModelMetadata
	err = json.Unmarshal(data, &metadata)
	if err != nil {
		return nil, fmt.Errorf("failed to parse metadata JSON: %w", err)
	}

	return &metadata, nil
}

// SaveMetadata saves metadata to a JSON file
func (m *ModelMetadata) Save(filepath string) error {
	data, err := json.MarshalIndent(m, "", "  ")
	if err != nil {
		return fmt.Errorf("failed to marshal metadata: %w", err)
	}

	err = os.WriteFile(filepath, data, 0644)
	if err != nil {
		return fmt.Errorf("failed to write metadata file: %w", err)
	}

	return nil
}

// BaseModelInfo contains information from the base model registry
type BaseModelInfo struct {
	ID          string `json:"id"`
	Name        string `json:"name"`
	Description string `json:"description"`
	Scenario    string `json:"scenario"`
}
