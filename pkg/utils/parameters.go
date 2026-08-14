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

// stringParams is the set of parameters whose VALUE is genuinely a string — a
// spec, table, list, or module-selector name — even when that value happens to be
// all-numeric (e.g. movement_barriers=3, a bare terrain code the §3 barriers doc
// explicitly supports). These are ALWAYS stored in model.StringParams and read
// back via Model.StringParam().
//
// The bug this fixes: routing by strconv.ParseFloat success (the old behavior)
// silently dropped an all-numeric value into the numeric Parameters map, where
// StringParam() could never see it — so the feature was inert (movement_barriers=3
// disabled, while movement_barriers=3:block worked). This holds for every string
// param shipped across §1–§3 (mutation_classes, mutation_rate_map, migration_matrix,
// habitat_suitability, the *_demes slots, …).
//
// This set is exactly the format=string rows of parameter_defaults.csv EXCEPT the
// two legacy numeric-coded dropdowns — mating_style and selection — which store an
// integer code in Parameters and resolve through the legacy-int mapping in
// pkg/modules (moving them would break the default distance-mating run). username /
// model_name / scenario / map_name are handled specially in SetParameter/LoadParameters
// before routing and are not listed here. setup_style is a module selector consumed
// via StringParam but not present in parameter_defaults.csv; it is included defensively.
//
// TestStringParamAllowlist derives the expected set from parameter_defaults.csv and
// asserts this map stays in sync, so a future string param cannot be silently forgotten.
var stringParams = map[string]bool{
	"environmental_events": true,
	"selection_mode":       true,
	"fitness_model":        true,
	"mutation_classes":     true,
	"mutation_rate_map":    true,
	"linkage_model":        true,
	"mutation_count_model": true,
	"dfe_model":            true,
	"recombination_model":  true,
	"migration_matrix":     true,
	"geographic_fst_grid":  true,
	"fstats_demes":         true,
	"joint_sfs_demes":      true,
	"dating_mu_factors":    true,
	"dating_gen_times":     true,
	"haplostats_demes":     true,
	"deme_inheritance":     true,
	"arg_loci":             true,
	"founder_allele_model": true,
	"habitat_suitability":  true,
	"movement_barriers":    true,
	"birth_style":          true,
	"death_style":          true,
	"setup_style":          true,
}

// IsStringParam reports whether paramName is a genuinely-string parameter that must
// be routed to model.StringParams regardless of whether its value parses as a number.
func IsStringParam(paramName string) bool {
	return stringParams[paramName]
}

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

	// Genuinely-string parameters are routed by declared type, not by ParseFloat
	// success, so an all-numeric spec value (movement_barriers=3, migration_matrix
	// edges, etc.) is stored where StringParam() can read it instead of being
	// silently dropped into the numeric Parameters map.
	if stringParams[paramName] {
		if model.StringParams == nil {
			model.StringParams = make(map[string]string)
		}
		model.StringParams[paramName] = paramValue
		return nil
	}

	// Remaining parameters are numeric; anything non-numeric that slips through
	// (e.g. an ad-hoc module selector not in the string set) is stored as a string.
	value, err := strconv.ParseFloat(paramValue, 64)
	if err != nil {
		if model.StringParams == nil {
			model.StringParams = make(map[string]string)
		}
		model.StringParams[paramName] = paramValue
		return nil
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

	// Check string parameters (e.g. module selectors like birth_style) for differences
	for paramName, paramValue := range model.StringParams {
		defaultVal, hasDefault := defaults[paramName]
		if !hasDefault || paramValue != defaultVal {
			if rec, exists := recordMap[paramName]; exists {
				row := make([]string, 6)
				copy(row, rec)
				if len(row) < 6 {
					row = append(row, "Main")
				}
				row[4] = paramValue
				writer.Write(row)
			} else {
				writer.Write([]string{paramName, paramName, "Dropdown", "string", paramValue, "Modules"})
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

// LoadChromosomesFromPath loads chromosomes from a specific path, parsing the
// same one-row-per-arm format as LoadChromosomes (delegates to a shared parser).
func LoadChromosomesFromPath(model *core.Model, path string) error {
	records, err := LoadCSV(path)
	if err != nil {
		return err
	}
	return parseChromosomeRecords(model, records)
}
