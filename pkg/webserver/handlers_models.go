// File: pkg/webserver/handlers_models.go
package webserver

import (
	"drift/pkg/config"
	"encoding/json"
	"fmt"
	"net/http"
)

// handleUserModelsList returns list of user's models
func handleUserModelsList(w http.ResponseWriter, r *http.Request) {
	username := r.URL.Query().Get("username")
	if username == "" {
		respondError(w, "Username required", http.StatusBadRequest)
		return
	}

	// Ensure user models exist
	err := config.GetOrCreateUser(username)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to initialize user: %v", err), http.StatusInternalServerError)
		return
	}

	models, err := config.ListUserModels(username)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to list models: %v", err), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"success": true,
		"models":  models,
	})
}

// handleModelSave saves current model parameters
func handleModelSave(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		respondError(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	var req struct {
		Username   string                 `json:"username"`
		ModelName  string                 `json:"model_name"`
		Parameters map[string]interface{} `json:"parameters"`
	}

	err := json.NewDecoder(r.Body).Decode(&req)
	if err != nil {
		respondError(w, "Invalid request body", http.StatusBadRequest)
		return
	}

	// Load current model
	model, err := config.LoadUserModel(req.Username, req.ModelName)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to load model: %v", err), http.StatusInternalServerError)
		return
	}

	// Update parameters
	for key, value := range req.Parameters {
		if floatVal, ok := value.(float64); ok {
			model.Parameters[key] = floatVal
		}
	}

	// Save model
	err = config.SaveUserModel(req.Username, req.ModelName, model)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to save model: %v", err), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"success": true,
		"message": "Model saved successfully",
	})
}

// handleModelReset resets model to base defaults
func handleModelReset(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		respondError(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	var req struct {
		Username  string `json:"username"`
		ModelName string `json:"model_name"`
	}

	err := json.NewDecoder(r.Body).Decode(&req)
	if err != nil {
		respondError(w, "Invalid request body", http.StatusBadRequest)
		return
	}

	err = config.ResetUserModel(req.Username, req.ModelName)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to reset model: %v", err), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"success": true,
		"message": "Model reset to defaults",
	})
}

// handleModelLoad loads model parameters
func handleModelLoad(w http.ResponseWriter, r *http.Request) {
	username := r.URL.Query().Get("username")
	modelName := r.URL.Query().Get("model_name")

	if username == "" || modelName == "" {
		respondError(w, "Username and model_name required", http.StatusBadRequest)
		return
	}

	model, err := config.LoadUserModel(username, modelName)
	if err != nil {
		respondError(w, fmt.Sprintf("Failed to load model: %v", err), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"success":    true,
		"parameters": model.Parameters,
		"scenario":   model.Scenario,
		"base_model": model.BaseModelID,
	})
}
