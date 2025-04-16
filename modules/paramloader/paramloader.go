package paramloader

import (
	"drift/types"
	"encoding/csv"
	"os"
	"strconv"
)

func LoadParameters(model *types.Model) {

	file, err := os.Open("C:/Go/Programs/Drift/static/parameter_defaults.csv")
	if err != nil {
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
	}

	for _, record := range records[1:] {
		// Convert the value based on the group
		if record[0] == "model_name" {
			model.ModelName = record[0]
		} else if record[5] == "Plot" {
			boolValue, err := strconv.ParseBool(record[4])
			if err != nil {
			}
			model.PlotFlags[record[0]] = boolValue
		} else {
			value, err := strconv.ParseFloat(record[4], 64)
			if err != nil {
			}
			model.Parameters[record[0]] = value
		}
	}
}
