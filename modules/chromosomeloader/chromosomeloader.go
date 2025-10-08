package chromosomeloader

import (
	"drift/pkg/utils"
	"drift/pkg/core"
)

// Name of the CSV file containing the chromosome arms.
const myFileName = "chromosome_data.csv"

// Load the chromosome arms from a CSV file and populate the model's ChromosomeArms map.
func LoadChromosomeArms(model *core.Model, configRoot string) error {
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
	var maxStart, maxStartLength int
	for _, record := range records[1:] {
		// Ensure the record has at least 4 fields
		err := csvLoader.CheckRecord(record, 4)
		if err != nil {
			return err
		}

		chromosome, err := csvLoader.Atoi(record, 0)
		if err != nil {
			return err
		}

		arm, err := csvLoader.Atoi(record, 1)
		if err != nil {
			return err
		}

		start, err := csvLoader.Atoi(record, 2)
		if err != nil {
			return err
		}

		length, err := csvLoader.Atoi(record, 3)
		if err != nil {
			return err
		}

		// Check if this is the highest start position so we can see how many bits we need
		if start > maxStart {
			maxStart = start
			maxStartLength = length
		}

		if model.ChromosomeArms[chromosome] == nil {
			model.ChromosomeArms[chromosome] = make(map[int][]int)
		}
		if model.ChromosomeArms[chromosome][arm] == nil {
			model.ChromosomeArms[chromosome][arm] = make([]int, 2)
		}
		model.ChromosomeArms[chromosome][arm][0] = start * int(model.Parameters["multiplier"])
		model.ChromosomeArms[chromosome][arm][1] = length * int(model.Parameters["multiplier"])
	}

	model.FreeParameters["genome_bits"] = maxStart + maxStartLength

	return nil
}
