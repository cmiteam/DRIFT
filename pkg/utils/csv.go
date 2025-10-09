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

// Generic CSV utilities
func LoadCSV(filename string) ([][]string, error) {
    file, err := os.Open(filename)
    if err != nil {
        return nil, fmt.Errorf("failed to open %s: %w", filename, err)
    }
    defer file.Close()

    reader := csv.NewReader(file)
    records, err := reader.ReadAll()
    if err != nil {
        return nil, fmt.Errorf("failed to read %s: %w", filename, err)
    }

    return records, nil
}

func StripHeader(records [][]string) [][]string {
    if len(records) <= 1 {
        return [][]string{}
    }
    return records[1:]
}

func ParseFloat(value string, row, col int) (float64, error) {
    value = strings.TrimSpace(value)
    f, err := strconv.ParseFloat(value, 64)
    if err != nil {
        return 0, fmt.Errorf("row %d, col %d: invalid float '%s'", row, col, value)
    }
    return f, nil
}

func ParseInt(value string, row, col int) (int, error) {
    value = strings.TrimSpace(value)
    i, err := strconv.Atoi(value)
    if err != nil {
        return 0, fmt.Errorf("row %d, col %d: invalid int '%s'", row, col, value)
    }
    return i, nil
}

// LoadParameters loads parameter_defaults.csv into model.Parameters
func LoadParameters(model *core.Model, configRoot string) error {
    filename := filepath.Join(configRoot, "parameter_defaults.csv")
    records, err := LoadCSV(filename)
    if err != nil {
        return err
    }

    // Skip header row
    records = StripHeader(records)
    
    // Initialize maps
    model.Parameters = make(map[string]float64)
    model.PlotFlags = make(map[string]bool)

    for i, record := range records {
        // Ensure record has at least 5 fields
        if len(record) < 5 {
            return fmt.Errorf("row %d: expected at least 5 fields, got %d", i+2, len(record))
        }

        paramName := strings.TrimSpace(record[0])
        paramValue := strings.TrimSpace(record[4])
        
        // Check if there's a group field (column 5)
        var group string
        if len(record) > 5 {
            group = strings.TrimSpace(record[5])
        }

        // Handle special string parameters
        if paramName == "model_name" {
            model.ModelName = paramValue
            continue
        }
        
        if paramName == "map_name" {
            model.MapName = paramValue
            continue
        }

        // Handle plot flags
        if group == "Plot" {
            boolValue, err := strconv.ParseBool(paramValue)
            if err != nil {
                return fmt.Errorf("row %d: parameter %s: invalid bool '%s'", i+2, paramName, paramValue)
            }
            model.PlotFlags[paramName] = boolValue
            continue
        }

        // Handle numeric parameters
        value, err := strconv.ParseFloat(paramValue, 64)
        if err != nil {
            return fmt.Errorf("row %d: parameter %s: invalid float '%s'", i+2, paramName, paramValue)
        }
        model.Parameters[paramName] = value
    }

    return nil
}

// LoadChromosomes loads chromosome_data.csv into model.ChromosomeArms
func LoadChromosomes(model *core.Model, configRoot string) error {
    filename := filepath.Join(configRoot, "chromosome_data.csv")
    records, err := LoadCSV(filename)
    if err != nil {
        return err
    }

    records = StripHeader(records)
    model.ChromosomeArms = make(map[int]map[int][]int)

    for i, row := range records {
        if len(row) < 4 {
            continue
        }

        chromNum, err := ParseInt(row[0], i+1, 0)
        if err != nil {
            return err
        }

        arm, err := ParseInt(row[1], i+1, 1)  // 0 = p arm, 1 = q arm
        if err != nil {
            return err
        }

        start, err := ParseInt(row[2], i+1, 2)
        if err != nil {
            return err
        }

        length, err := ParseInt(row[3], i+1, 3)
        if err != nil {
            return err
        }

        // Initialize chromosome if needed
        if model.ChromosomeArms[chromNum] == nil {
            model.ChromosomeArms[chromNum] = make(map[int][]int)
        }

        // Store: model.ChromosomeArms[chromosome][arm] = [start, length]
        model.ChromosomeArms[chromNum][arm] = []int{start, length}
    }

    return nil
}

// LoadActuarialTable loads actuarial_table.csv into model.DeathRisk
func LoadActuarialTable(model *core.Model, configRoot string) error {
    filename := filepath.Join(configRoot, "actuarial_table.csv")
    records, err := LoadCSV(filename)
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

        model.DeathRisk[age] = risk
    }

    return nil
}