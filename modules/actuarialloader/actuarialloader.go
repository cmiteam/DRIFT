package actuarialloader

import (
	"drift/pkg/core"
	"drift/pkg/utils"
)

// Name of the CSV file containing the actuarial table.
const myFileName = "actuarial_table.csv"

// Load the actuarial table from a CSV file and populate the model's DeathRisk and CumulativeProb maps.
func LoadActuarialTable(model *core.Model, configRoot string) error {
	// Load the CSV file
	csvLoader := utils.CSVLoader{
		FileName:   myFileName,
		Dir:        configRoot,
		MinRecords: 2,
	}
	records, err := csvLoader.LoadCSV()
	if err != nil {
		return err
	}

	// Skip the header row and process each record
	var cumulative float64
	for _, record := range records[1:] { // Skip the header row
		// Ensure the record has at least 3 fields
		err := csvLoader.CheckRecord(record, 3)
		if err != nil {
			return err
		}

		age, err := csvLoader.Atoi(record, 0)
		if err != nil {
			return err
		}

		risk, err := csvLoader.ParseFloat64(record, 1)
		if err != nil {
			return err
		}
		model.DeathRisk[age] = risk

		popProb, err := csvLoader.ParseFloat64(record, 2)
		if err != nil {
			return err
		}
		cumulative += popProb * 4.4
		model.CumulativeProb[age] = cumulative
	}
	return nil
}
