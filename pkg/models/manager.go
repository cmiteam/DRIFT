// File: pkg/models/manager.go
package models

import (
	"drift/pkg/core"
	"drift/pkg/utils"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"time"
)

// ModelManager handles loading, saving, and managing models
type ModelManager struct {
	BaseModelsPath string // e.g., "static/basemodels"
	UserModelsPath string // e.g., "users"
}

// NewModelManager creates a new ModelManager
func NewModelManager(baseModelsPath, userModelsPath string) *ModelManager {
	return &ModelManager{
		BaseModelsPath: baseModelsPath,
		UserModelsPath: userModelsPath,
	}
}

// LoadModel loads a user model with inheritance from base model
func (m *ModelManager) LoadModel(username, modelName string) (*core.Model, error) {
	userModelPath := filepath.Join(m.UserModelsPath, username, "models", modelName)
	metadataPath := filepath.Join(userModelPath, "metadata.json")

	// Load metadata to get base model
	metadata, err := LoadMetadata(metadataPath)
	if err != nil {
		return nil, fmt.Errorf("failed to load model metadata: %w", err)
	}

	// Create empty model
	model := &core.Model{
		Parameters:      make(map[string]float64),
		PlotFlags:       make(map[string]bool),
		ChromosomeArms:  make(map[int]map[int][]int),
		DeathRisk:       make(map[int]float64),
		CumulativeProb:  make(map[int]float64),
		FreeParameters:  make(map[string]int),
		Map:             nil,
		PrevFrequencies: make(map[int]float64),
		HetHistory:      make([]float64, 0),
		TimeHistory:     make([]int, 0),
		ResultsDir:      filepath.Join(userModelPath, "results"),
		BaseModelID:     metadata.BaseModel,
		Username:        username,
		ModelName:       modelName,
	}

	// Load base model info to get scenario
	baseInfo, err := m.GetBaseModelInfo(metadata.BaseModel)
	if err != nil {
		return nil, err
	}
	model.Scenario = baseInfo.Scenario

	// Load parameters with inheritance
	configPath := filepath.Join(userModelPath, "config")
	err = m.loadParametersWithInheritance(model, configPath, metadata.BaseModel)
	if err != nil {
		return nil, err
	}

	// Load actuarial table (try user override first, then Default base model)
	actuarialPath := filepath.Join(configPath, "actuarial_table.csv")
	if _, err := os.Stat(actuarialPath); os.IsNotExist(err) {
		actuarialPath = filepath.Join(m.BaseModelsPath, "Default", "actuarial_table.csv")
	}
	err = utils.LoadActuarialTableFromPath(model, actuarialPath)
	if err != nil {
		return nil, err
	}

	// Load chromosomes (try user override first, then Default base model)
	chromosomePath := filepath.Join(configPath, "chromosome_data.csv")
	if _, err := os.Stat(chromosomePath); os.IsNotExist(err) {
		chromosomePath = filepath.Join(m.BaseModelsPath, "Default", "chromosome_data.csv")
	}
	err = utils.LoadChromosomesFromPath(model, chromosomePath)
	if err != nil {
		return nil, err
	}

	return model, nil
}

// loadParametersWithInheritance loads parameters with inheritance support
func (m *ModelManager) loadParametersWithInheritance(model *core.Model, configPath, baseModelID string) error {
	// First, load parameter_defaults.csv (full set of parameters)
	_, err := utils.LoadParameterDefaults(model)
	if err != nil {
		return fmt.Errorf("failed to load parameter defaults: %w", err)
	}

	// Then, load user overrides on top
	userParamsPath := filepath.Join(configPath, "parameters.csv")
	err = utils.LoadParameterOverrides(model, userParamsPath)
	if err != nil {
		return fmt.Errorf("failed to load user parameters: %w", err)
	}

	return nil
}

// LoadParameterFile loads a single parameter CSV file
func (m *ModelManager) LoadParameterFile(model *core.Model, filePath string, processInherit bool) error {
	file, err := os.Open(filePath)
	if err != nil {
		return err
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

		if len(record) < 2 {
			continue
		}

		paramName := strings.TrimSpace(record[0])
		paramValue := strings.TrimSpace(record[1])

		// Handle inheritance
		if paramName == "_inherit" && processInherit {
			inheritFrom := paramValue
			// Resolve relative to current file's directory
			dir := filePath[:strings.LastIndex(filePath, string(os.PathSeparator))]
			parentDir := dir[:strings.LastIndex(dir, string(os.PathSeparator))]
			inheritPath := filepath.Join(parentDir, inheritFrom, "parameters.csv")

			err = m.LoadParameterFile(model, inheritPath, true)
			if err != nil {
				return fmt.Errorf("failed to load inherited parameters from %s: %w", inheritFrom, err)
			}
			continue
		}

		// Parse and set parameter value
		err = utils.SetParameter(model, paramName, paramValue)
		if err != nil {
			return err
		}
	}

	return nil
}

// SaveModel saves the current model state to user directory
func (m *ModelManager) SaveModel(username, modelName string, model *core.Model) error {
	userModelPath := filepath.Join(m.UserModelsPath, username, "models", modelName)
	configPath := filepath.Join(userModelPath, "config")
	metadataPath := filepath.Join(userModelPath, "metadata.json")

	// Load existing metadata to preserve creation date
	metadata, err := LoadMetadata(metadataPath)
	if err != nil {
		return fmt.Errorf("failed to load metadata: %w", err)
	}

	// Update last modified
	metadata.LastModified = time.Now()

	// Save metadata
	err = metadata.Save(metadataPath)
	if err != nil {
		return err
	}

	// Save parameters
	paramsPath := filepath.Join(configPath, "parameters.csv")
	err = utils.SaveParameters(model, paramsPath)
	if err != nil {
		return fmt.Errorf("failed to save parameters: %w", err)
	}

	return nil
}

// ResetModel resets a user model to its base model defaults
func (m *ModelManager) ResetModel(username, modelName string) error {
	userModelPath := filepath.Join(m.UserModelsPath, username, "models", modelName)
	metadataPath := filepath.Join(userModelPath, "metadata.json")

	// Load metadata to get base model
	metadata, err := LoadMetadata(metadataPath)
	if err != nil {
		return fmt.Errorf("failed to load metadata: %w", err)
	}

	userConfigPath := filepath.Join(userModelPath, "config")
	userParamsPath := filepath.Join(userConfigPath, "parameters.csv")

	// Create a temporary model to load defaults + base model overrides
	tempModel := &core.Model{
		Parameters: make(map[string]float64),
	}

	// Load parameter_defaults.csv first
	_, err = utils.LoadParameterDefaults(tempModel)
	if err != nil {
		return fmt.Errorf("failed to load parameter defaults: %w", err)
	}

	// Apply base model overrides on top
	baseParamsPath := filepath.Join(m.BaseModelsPath, metadata.BaseModel, "parameters.csv")
	err = utils.LoadParameterOverrides(tempModel, baseParamsPath)
	if err != nil {
		return fmt.Errorf("failed to load base model overrides: %w", err)
	}

	// Save as user overrides (will only save what differs from defaults)
	err = utils.SaveParameters(tempModel, userParamsPath)
	if err != nil {
		return fmt.Errorf("failed to save user parameters: %w", err)
	}

	// Update metadata
	metadata.LastModified = time.Now()
	err = metadata.Save(metadataPath)
	if err != nil {
		return err
	}

	return nil
}

// CreateModel creates a new user model based on a base model
func (m *ModelManager) CreateModel(username, modelName, baseModelID, description string) error {
	userModelPath := filepath.Join(m.UserModelsPath, username, "models", modelName)
	configPath := filepath.Join(userModelPath, "config")
	resultsPath := filepath.Join(userModelPath, "results")

	// Create directories
	err := os.MkdirAll(configPath, 0755)
	if err != nil {
		return fmt.Errorf("failed to create config directory: %w", err)
	}
	err = os.MkdirAll(resultsPath, 0755)
	if err != nil {
		return fmt.Errorf("failed to create results directory: %w", err)
	}

	// Create a temporary model to load defaults + base model overrides
	tempModel := &core.Model{
		Parameters: make(map[string]float64),
	}

	// Load parameter_defaults.csv first
	_, err = utils.LoadParameterDefaults(tempModel)
	if err != nil {
		return fmt.Errorf("failed to load parameter defaults: %w", err)
	}

	// Apply base model overrides on top
	baseParamsPath := filepath.Join(m.BaseModelsPath, baseModelID, "parameters.csv")
	err = utils.LoadParameterOverrides(tempModel, baseParamsPath)
	if err != nil {
		return fmt.Errorf("failed to load base model overrides: %w", err)
	}

	// Save as user overrides (will only save what differs from defaults)
	userParamsPath := filepath.Join(configPath, "parameters.csv")
	err = utils.SaveParameters(tempModel, userParamsPath)
	if err != nil {
		return fmt.Errorf("failed to save user parameters: %w", err)
	}

	// Create metadata
	metadata := &ModelMetadata{
		Name:         modelName,
		Description:  description,
		BaseModel:    baseModelID,
		Created:      time.Now(),
		LastModified: time.Now(),
		Author:       username,
	}

	metadataPath := filepath.Join(userModelPath, "metadata.json")
	err = metadata.Save(metadataPath)
	if err != nil {
		return err
	}

	return nil
}

// CloneModel creates a copy of an existing user model
func (m *ModelManager) CloneModel(username, sourceModel, newModel, description string) error {
	sourcePath := filepath.Join(m.UserModelsPath, username, "models", sourceModel)
	destPath := filepath.Join(m.UserModelsPath, username, "models", newModel)

	// Copy entire directory
	err := m.copyDir(sourcePath, destPath)
	if err != nil {
		return fmt.Errorf("failed to copy model: %w", err)
	}

	// Load and update metadata
	metadataPath := filepath.Join(destPath, "metadata.json")
	metadata, err := LoadMetadata(metadataPath)
	if err != nil {
		return err
	}

	metadata.Name = newModel
	if description != "" {
		metadata.Description = description
	}
	metadata.Created = time.Now()
	metadata.LastModified = time.Now()

	err = metadata.Save(metadataPath)
	if err != nil {
		return err
	}

	return nil
}

// ListUserModels returns all models for a user
func (m *ModelManager) ListUserModels(username string) ([]ModelMetadata, error) {
	userPath := filepath.Join(m.UserModelsPath, username, "models")

	entries, err := os.ReadDir(userPath)
	if err != nil {
		if os.IsNotExist(err) {
			return []ModelMetadata{}, nil
		}
		return nil, err
	}

	var models []ModelMetadata
	for _, entry := range entries {
		if !entry.IsDir() {
			continue
		}

		metadataPath := filepath.Join(userPath, entry.Name(), "metadata.json")
		metadata, err := LoadMetadata(metadataPath)
		if err != nil {
			// Skip if no valid metadata
			continue
		}

		models = append(models, *metadata)
	}

	return models, nil
}

// ListBaseModels returns all available base models
func (m *ModelManager) ListBaseModels() ([]BaseModelInfo, error) {
	registryPath := filepath.Join(m.BaseModelsPath, "registry.csv")

	file, err := os.Open(registryPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open registry: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)

	// Skip header
	_, err = reader.Read()
	if err != nil {
		return nil, err
	}

	var baseModels []BaseModelInfo
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}

		if len(record) < 4 {
			continue
		}

		baseModels = append(baseModels, BaseModelInfo{
			ID:          strings.TrimSpace(record[0]),
			Name:        strings.TrimSpace(record[1]),
			Description: strings.TrimSpace(record[2]),
			Scenario:    strings.TrimSpace(record[3]),
		})
	}

	return baseModels, nil
}

// GetBaseModelInfo returns info for a specific base model
func (m *ModelManager) GetBaseModelInfo(baseModelID string) (*BaseModelInfo, error) {
	baseModels, err := m.ListBaseModels()
	if err != nil {
		return nil, err
	}

	for _, bm := range baseModels {
		if bm.ID == baseModelID {
			return &bm, nil
		}
	}

	return nil, fmt.Errorf("base model '%s' not found", baseModelID)
}

// Helper functions

func (m *ModelManager) copyFile(src, dst string) error {
	data, err := os.ReadFile(src)
	if err != nil {
		return err
	}
	return os.WriteFile(dst, data, 0644)
}

func (m *ModelManager) copyDir(src, dst string) error {
	return filepath.Walk(src, func(path string, info os.FileInfo, err error) error {
		if err != nil {
			return err
		}

		// Get relative path
		relPath, err := filepath.Rel(src, path)
		if err != nil {
			return err
		}

		dstPath := filepath.Join(dst, relPath)

		if info.IsDir() {
			return os.MkdirAll(dstPath, info.Mode())
		}

		return m.copyFile(path, dstPath)
	})
}
