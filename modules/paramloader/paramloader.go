package paramloader

import (
	"drift/types"
	"fmt"
	"strconv"
)

func LoadParameters(model *types.Model,
	records [][]string) error {

	// Ensure the CSV file has at least one record
	if len(records) == 0 {
		return fmt.Errorf("No records found in the CSV file")
	}
	// Skip the header row and process each record
	for _, record := range records[1:] {

		// Ensure the record has at least 5 fields
		if len(record) < 5 {
			return fmt.Errorf("Invalid record: %v", record)
		}
		// Convert the value based on the group
		if record[0] == "model_name" {
			model.ModelName = record[0]
		} else if record[5] == "Plot" {
			boolValue, err := strconv.ParseBool(record[4])
			if err != nil {
				return err
			}
			model.PlotFlags[record[0]] = boolValue
		} else {
			value, err := strconv.ParseFloat(record[4], 64)
			if err != nil {
				return err
			}
			model.Parameters[record[0]] = value
		}
	}

	return nil
}
