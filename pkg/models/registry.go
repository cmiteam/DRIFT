// File: pkg/models/registry.go
package models

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strings"
)

// Registry holds base model information
type Registry struct {
	models map[string]*BaseModelInfo
}

// NewRegistry creates a new registry
func NewRegistry() *Registry {
	return &Registry{
		models: make(map[string]*BaseModelInfo),
	}
}

// LoadFromFile loads registry from a CSV file
func (r *Registry) LoadFromFile(filepath string) error {
	file, err := os.Open(filepath)
	if err != nil {
		return fmt.Errorf("failed to open registry file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)

	// Skip header
	_, err = reader.Read()
	if err != nil {
		return err
	}

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		if len(record) < 4 {
			continue
		}

		info := &BaseModelInfo{
			ID:          strings.TrimSpace(record[0]),
			Name:        strings.TrimSpace(record[1]),
			Description: strings.TrimSpace(record[2]),
			Scenario:    strings.TrimSpace(record[3]),
		}

		r.models[info.ID] = info
	}

	return nil
}

// Get returns a base model by ID
func (r *Registry) Get(id string) (*BaseModelInfo, error) {
	info, exists := r.models[id]
	if !exists {
		return nil, fmt.Errorf("base model '%s' not found in registry", id)
	}
	return info, nil
}

// List returns all base models
func (r *Registry) List() []BaseModelInfo {
	var list []BaseModelInfo
	for _, info := range r.models {
		list = append(list, *info)
	}
	return list
}
