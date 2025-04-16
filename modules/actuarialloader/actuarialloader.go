package actuarialloader

import (
	"drift/types"
	"fmt"
	"strconv"
)

func LoadActuarialTable(model *types.Model,
	records [][]string) error {

	// Ensure the CSV file has at least one record
	if len(records) == 0 {
		return fmt.Errorf("No records found in the CSV file")
	}

	// Skip the header row and process each record
	var cumulative float64
	for _, row := range records[1:] { // Skip the header row

		// Ensure the record has at least 3 fields
		if len(row) < 3 {
			return fmt.Errorf("Invalid record: %v", row)
		}

		age, err := strconv.Atoi(row[0])
		if err != nil {
			return err
		}

		risk, err := strconv.ParseFloat(row[1], 64)
		if err != nil {
			return err
		}
		model.DeathRisk[age] = risk

		popProb, err := strconv.ParseFloat(row[2], 64)
		if err != nil {
			return err
		}
		cumulative += popProb * 4.4
		model.CumulativeProb[age] = cumulative
	}
	return nil
}
