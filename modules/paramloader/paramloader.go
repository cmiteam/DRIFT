package paramloader

import (
	"drift/modules/csvutils"
	"drift/types"
)

const myFileName = "parameter_defaults.csv"

// Load the parameters from a CSV file and populate the model's Parameters map.
func LoadParameters(model *types.Model, configRoot string) error {
	// Load the CSV file
	csvLoader := csvutils.CSVLoader{
		FileName:   myFileName,
		Dir:        configRoot,
		MinRecords: 2,
	}
	records, err := csvLoader.LoadCSV()
	if err != nil {
		return err
	}

	// Skip the header row and process each record
	for _, record := range records[1:] {
		// Ensure the record has at least 5 fields
		err := csvLoader.CheckRecord(record, 5)
		if err != nil {
			return err
		}

		// Convert the value based on the group
		if record[0] == "model_name" {
			model.ModelName = record[4]
		} else if record[0] == "map_name" {
			model.MapName = record[4]
		} else if record[0] == "mating_style" {
			model.MatingStyle = record[4]
		} else if record[5] == "Plot" {
			boolValue, err := csvLoader.ParseBool(record, 4)
			if err != nil {
				return err
			}
			model.PlotFlags[record[0]] = boolValue
		} else {
			value, err := csvLoader.ParseFloat64(record, 4)
			if err != nil {
				return err
			}
			model.Parameters[record[0]] = value
		}
	}
	return nil
}
