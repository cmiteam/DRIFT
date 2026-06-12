// File: pkg/config/modelmanager.go
package config

import (
	"drift/pkg/core"
	"drift/pkg/models"
)

var globalModelManager *models.ModelManager

// InitializeModelManager initializes the global model manager
func InitializeModelManager(baseModelsPath, userModelsPath string) {
	globalModelManager = models.NewModelManager(baseModelsPath, userModelsPath)
}

// GetModelManager returns the global model manager
func GetModelManager() *models.ModelManager {
	if globalModelManager == nil {
		// Use default paths
		globalModelManager = models.NewModelManager("static/basemodels", "users")
	}
	return globalModelManager
}

// LoadUserModel loads a user model using the ModelManager
func LoadUserModel(username, modelName string) (*core.Model, error) {
	manager := GetModelManager()
	return manager.LoadModel(username, modelName)
}

// CreateUserModel creates a new user model from a base model
func CreateUserModel(username, modelName, baseModelID, description string) error {
	manager := GetModelManager()
	return manager.CreateModel(username, modelName, baseModelID, description)
}

// SaveUserModel saves a user model
func SaveUserModel(username, modelName string, model *core.Model) error {
	manager := GetModelManager()
	return manager.SaveModel(username, modelName, model)
}

// ResetUserModel resets a user model to base defaults
func ResetUserModel(username, modelName string) error {
	manager := GetModelManager()
	return manager.ResetModel(username, modelName)
}

// CloneUserModel clones a user model
func CloneUserModel(username, sourceModel, newModel, description string) error {
	manager := GetModelManager()
	return manager.CloneModel(username, sourceModel, newModel, description)
}

// ListUserModels lists all models for a user
func ListUserModels(username string) ([]models.ModelMetadata, error) {
	manager := GetModelManager()
	return manager.ListUserModels(username)
}

// ListBaseModels lists all available base models
func ListBaseModels() ([]models.BaseModelInfo, error) {
	manager := GetModelManager()
	return manager.ListBaseModels()
}

// InitializeModelFromBaseModel creates a model from base model without saving
// This is useful for temporary models or testing
func InitializeModelFromBaseModel(baseModelID string) (*core.Model, error) {
	manager := GetModelManager()

	// Load base model using temp path as user model
	// This creates a temporary in-memory model without user modifications
	tempModel, err := manager.LoadModel("_temp", "_temp")
	if err != nil {
		// If temp doesn't exist, we need to load directly from base
		// For now, return error - this function should be used with existing user models
		return nil, err
	}

	return tempModel, nil
}

// GetOrCreateUser ensures a user exists with all base models initialized
func GetOrCreateUser(username string) error {
	manager := GetModelManager()
	return manager.GetOrCreateUser(username)
}

// InitializeUserModels initializes all base models for a new user
func InitializeUserModels(username string) error {
	manager := GetModelManager()
	return manager.InitializeUserModels(username)
}
