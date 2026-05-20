// File: pkg/utils/parameters.go
package utils

import (
	"drift/pkg/core"
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// ParameterRecord holds the full metadata for a parameter from parameter_defaults.csv
type ParameterRecord struct {
	Name   string
	Label  string
	Entry  string
	Format string
	Value  string
	Group  string
}

// LoadParameterDefaults loads parameter_defaults.csv and returns both the records and applies values to model
func LoadParameterDefaults(model *core.Model) ([]ParameterRecord, error) {
	defaultsPath := filepath.Join("static", "parameter_defaults.csv")
	return LoadParameterDefaultsFromPath(model, defaultsPath)
}

// LoadParameterDefaultsFromPath loads parameter defaults from a specific path
func LoadParameterDefaultsFromPath(model *core.Model, path string) ([]ParameterRecord, error) {
	records, err := LoadCSV(path)
	if err != nil {
		return nil, fmt.Errorf("failed to load parameter defaults: %w", err)
	}

	records = StripHeader(records)
	var paramRecords []ParameterRecord

	for _, row := range records {
		if len(row) < 5 {
			continue
		}

		rec := ParameterRecord{
			Name:   strings.TrimSpace(row[0]),
			Label:  row[1],
			Entry:  row[2],
			Format: row[3],
			Value:  row[4],
		}
		if len(row) >= 6 {
			rec.Group = row[5]
		}

		paramRecords = append(paramRecords, rec)

		// Apply value to model
		err := SetParameter(model, rec.Name, rec.Value)
		if err != nil {
			// Log but don't fail - some parameters might be special
			continue
		}
	}

	return paramRecords, nil
}

// LoadParameterOverrides loads an override file and applies values on top of existing model parameters.
// Accepts either the 2-column schema (parameter,value) used by base models or the 6-column schema
// (parameter,label,entry,format,value,group) used by parameter_defaults.csv and saved user models.
func LoadParameterOverrides(model *core.Model, path string) error {
	if _, err := os.Stat(path); os.IsNotExist(err) {
		return nil // No overrides file is OK
	}

	records, err := LoadCSV(path)
	if err != nil {
		return fmt.Errorf("failed to load parameter overrides: %w", err)
	}

	records = StripHeader(records)

	for _, row := range records {
		if len(row) < 2 {
			continue
		}

		paramName := strings.TrimSpace(row[0])

		// Value column depends on schema: 6-col (parameter_defaults.csv) puts value at index 4,
		// 2-col (base model parameters.csv) puts it at index 1.
		var paramValue string
		if len(row) >= 5 {
			paramValue = strings.TrimSpace(row[4])
		} else {
			paramValue = strings.TrimSpace(row[1])
		}

		err := SetParameter(model, paramName, paramValue)
		if err != nil {
			continue
		}
	}

	return nil
}

// GetDefaultValues returns a map of parameter name -> default value from parameter_defaults.csv
func GetDefaultValues() (map[string]string, error) {
	defaultsPath := filepath.Join("static", "parameter_defaults.csv")
	records, err := LoadCSV(defaultsPath)
	if err != nil {
		return nil, err
	}

	records = StripHeader(records)
	defaults := make(map[string]string)

	for _, row := range records {
		if len(row) >= 5 {
			defaults[strings.TrimSpace(row[0])] = strings.TrimSpace(row[4])
		}
	}

	return defaults, nil
}

// SetParameter sets a single parameter on the model
func SetParameter(model *core.Model, paramName, paramValue string) error {
	paramName = strings.TrimSpace(paramName)
	paramValue = strings.TrimSpace(paramValue)

	// Handle special string parameters
	if paramName == "username" {
		model.Username = paramValue
		return nil
	}
	if paramName == "model_name" {
		model.ModelName = paramValue
		return nil
	}
	if paramName == "scenario" {
		// Only set scenario if not already set
		if model.Scenario == "" {
			model.Scenario = paramValue
		}
		return nil
	}
	if paramName == "map_name" {
		model.MapName = paramValue
		return nil
	}

	// Handle numeric parameters
	value, err := strconv.ParseFloat(paramValue, 64)
	if err != nil {
		return fmt.Errorf("parameter %s: invalid float '%s': %w", paramName, paramValue, err)
	}
	model.Parameters[paramName] = value

	return nil
}

// SaveParameters saves model parameters to a CSV file (only overrides that differ from defaults)
func SaveParameters(model *core.Model, destPath string) error {
	// Load defaults to compare against
	defaultsPath := filepath.Join("static", "parameter_defaults.csv")
	defaultRecords, err := LoadCSV(defaultsPath)
	if err != nil {
		return fmt.Errorf("failed to load parameter defaults: %w", err)
	}
	defaultRecords = StripHeader(defaultRecords)

	// Build map of default values and full records
	defaults := make(map[string]string)
	recordMap := make(map[string][]string)
	for _, row := range defaultRecords {
		if len(row) >= 5 {
			name := strings.TrimSpace(row[0])
			defaults[name] = strings.TrimSpace(row[4])
			recordMap[name] = row
		}
	}

	// Create output file
	file, err := os.Create(destPath)
	if err != nil {
		return fmt.Errorf("failed to create parameters file: %w", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write header (same format as parameter_defaults.csv)
	err = writer.Write([]string{"parameter", "label", "entry", "format", "value", "group"})
	if err != nil {
		return err
	}

	// Check string parameters for differences
	stringParams := map[string]string{
		"username":   model.Username,
		"model_name": model.ModelName,
		"scenario":   model.Scenario,
		"map_name":   model.MapName,
	}

	for paramName, paramValue := range stringParams {
		if paramValue == "" {
			continue
		}
		defaultVal, hasDefault := defaults[paramName]
		if !hasDefault || paramValue != defaultVal {
			// Write override with full metadata from defaults
			if rec, exists := recordMap[paramName]; exists {
				row := make([]string, 6)
				copy(row, rec)
				row[4] = paramValue // Update value
				writer.Write(row)
			} else {
				// Parameter not in defaults, write minimal record
				writer.Write([]string{paramName, paramName, "Text", "string", paramValue, "Main"})
			}
		}
	}

	// Check numeric parameters for differences
	for paramName, paramValue := range model.Parameters {
		// Format current value
		var valueStr string
		if paramValue == float64(int64(paramValue)) {
			valueStr = fmt.Sprintf("%d", int64(paramValue))
		} else {
			valueStr = fmt.Sprintf("%g", paramValue)
		}

		// Compare against default
		defaultVal, hasDefault := defaults[paramName]
		if !hasDefault || valueStr != defaultVal {
			// Write override with full metadata from defaults
			if rec, exists := recordMap[paramName]; exists {
				row := make([]string, 6)
				copy(row, rec)
				if len(row) < 6 {
					row = append(row, "Main")
				}
				row[4] = valueStr // Update value
				writer.Write(row)
			} else {
				// Parameter not in defaults, write minimal record
				writer.Write([]string{paramName, paramName, "Text", "float", valueStr, "Main"})
			}
		}
	}

	return nil
}

// LoadActuarialTableFromPath loads actuarial table from a specific path
func LoadActuarialTableFromPath(model *core.Model, filepath string) error {
	records, err := LoadCSV(filepath)
	if err != nil {
		return err
	}

	records = StripHeader(records)
	model.DeathRisk = make(map[int]float64)

	for i, row := range records {
		if len(row) < 2 {
			continue
		}

		age, err := ParseInt(row[0], i+1, 0)
		if err != nil {
			return err
		}

		risk, err := ParseFloat(row[1], i+1, 1)
		if err != nil {
			return err
		}

		model.DeathRisk[age] = risk / model.Parameters["death_risk_adj"]
	}

	return nil
}

// LoadChromosomesFromPath loads chromosomes from a specific path
func LoadChromosomesFromPath(model *core.Model, filepath string) error {
	records, err := LoadCSV(filepath)
	if err != nil {
		return err
	}

	records = StripHeader(records)
	model.ChromosomeArms = make(map[int]map[int][]int)

	for i, row := range records {
		if len(row) < 5 {
			continue
		}

		chrom, err := ParseInt(row[0], i+1, 0)
		if err != nil {
			return err
		}

		pStart, err := ParseInt(row[1], i+1, 1)
		if err != nil {
			return err
		}

		pLen, err := ParseInt(row[2], i+1, 2)
		if err != nil {
			return err
		}

		qStart, err := ParseInt(row[3], i+1, 3)
		if err != nil {
			return err
		}

		qLen, err := ParseInt(row[4], i+1, 4)
		if err != nil {
			return err
		}

		if model.ChromosomeArms[chrom] == nil {
			model.ChromosomeArms[chrom] = make(map[int][]int)
		}

		model.ChromosomeArms[chrom][0] = []int{pStart, pLen}
		model.ChromosomeArms[chrom][1] = []int{qStart, qLen}
	}

	return nil
}
