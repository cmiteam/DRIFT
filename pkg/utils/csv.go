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
	reader.FieldsPerRecord = -1 // tolerate rows with varying widths (e.g. mixed 2-col/6-col parameters.csv)
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
	model.StringParams = make(map[string]string)

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
		if paramName == "username" {
			model.Username = paramValue
			continue
		}
		if paramName == "model_name" {
			model.ModelName = paramValue
			continue
		}
		if paramName == "scenario" {
			model.Scenario = paramValue
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

		// Numeric parameters go in Parameters; non-numeric values (e.g. a module
		// selector like birth_style=standard) are stored as string parameters.
		value, err := strconv.ParseFloat(paramValue, 64)
		if err != nil {
			model.StringParams[paramName] = paramValue
			continue
		}
		model.Parameters[paramName] = value
	}

	return nil
}

// LoadChromosomes loads chromosome_data.csv into model.ChromosomeArms
func LoadChromosomes(model *core.Model, configRoot string) error {
	return LoadChromosomesFromPath(model, filepath.Join(configRoot, "chromosome_data.csv"))
}

// parseChromosomeRecords parses the chromosome CSV format shared by every load
// path: one row per arm — Chromosome, Arm (0=p, 1=q), Start, Length. Previously
// LoadChromosomesFromPath used a divergent 5-column-per-chromosome parser that
// silently skipped every row of these 4-column files, leaving ChromosomeArms
// empty and genome_bits = 0 on the user-model path.
func parseChromosomeRecords(model *core.Model, records [][]string) error {
	records = StripHeader(records)
	model.ChromosomeArms = make(map[int]map[int][]int)
	model.ArmCM = make(map[int]map[int]float64)

	for i, row := range records {
		if len(row) < 4 {
			continue
		}

		chromNum, err := ParseInt(row[0], i+1, 0)
		if err != nil {
			return err
		}
		arm, err := ParseInt(row[1], i+1, 1) // 0 = p arm, 1 = q arm
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

		if model.ChromosomeArms[chromNum] == nil {
			model.ChromosomeArms[chromNum] = make(map[int][]int)
		}
		model.ChromosomeArms[chromNum][arm] = []int{start, length}

		// Optional 5th column: the arm's genetic length in centiMorgans (roadmap
		// §6a). Present only in map-enabled chromosome files; a blank/absent cell
		// leaves this arm out of ArmCM so the map falls back to the uniform rate for
		// it. A malformed value is a hard error (fail-fast at load, like other params).
		if len(row) >= 5 && strings.TrimSpace(row[4]) != "" {
			cm, err := strconv.ParseFloat(strings.TrimSpace(row[4]), 64)
			if err != nil {
				return fmt.Errorf("chromosome_data row %d col 5 (cM): %w", i+1, err)
			}
			if model.ArmCM[chromNum] == nil {
				model.ArmCM[chromNum] = make(map[int]float64)
			}
			model.ArmCM[chromNum][arm] = cm
		}
	}

	SetGenomeBits(model)
	return nil
}

// SetGenomeBits computes genome_bits — the highest occupied bit position across
// all chromosome arms — from ChromosomeArms and stores it in FreeParameters. It
// must be called after chromosomes are loaded; both LoadChromosomes and
// LoadChromosomesFromPath call it so every load path sets genome_bits. (Before
// this, the user-model load path never set it, silently disabling all DNA /
// mutation tracking.) Returns the computed value.
func SetGenomeBits(model *core.Model) int {
	totalBits := 0
	for chrom := range model.ChromosomeArms {
		if arm0, ok := model.ChromosomeArms[chrom][0]; ok && len(arm0) >= 2 {
			if end := arm0[0] + arm0[1]; end > totalBits {
				totalBits = end
			}
		}
		if arm1, ok := model.ChromosomeArms[chrom][1]; ok && len(arm1) >= 2 {
			if end := arm1[0] + arm1[1]; end > totalBits {
				totalBits = end
			}
		}
	}
	if model.FreeParameters == nil {
		model.FreeParameters = make(map[string]int)
	}
	model.FreeParameters["genome_bits"] = totalBits
	return totalBits
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
