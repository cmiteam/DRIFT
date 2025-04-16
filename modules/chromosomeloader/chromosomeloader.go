package chromosomeloader

import (
	"drift/types"
	"fmt"
	"strconv"
)

func LoadChromosomeArmsFromCSV(model *types.Model,
	records [][]string) error {

	// Ensure the CSV file has at least one record
	if len(records) == 0 {
		return fmt.Errorf("No records found in the CSV file")
	}

	// Skip the header row and process each record
	var totallen int
	for _, record := range records[1:] {

		// Ensure the record has at least 4 fields
		if len(record) < 4 {
			return fmt.Errorf("Invalid record: %v", record)
		}

		chromosome, err := strconv.Atoi(record[0])
		if err != nil {
			return err
		}

		arm, err := strconv.Atoi(record[1])
		if err != nil {
			return err
		}

		start, err := strconv.Atoi(record[2])
		if err != nil {
			return err
		}

		length, err := strconv.Atoi(record[3])
		if err != nil {
			return err
		}

		if model.ChromosomeArms[chromosome] == nil {
			model.ChromosomeArms[chromosome] = make(map[int][]int)
		}
		if model.ChromosomeArms[chromosome][arm] == nil {
			model.ChromosomeArms[chromosome][arm] = make([]int, 2)
		}
		model.ChromosomeArms[chromosome][arm][0] = start * int(model.Parameters["multiplier"])
		model.ChromosomeArms[chromosome][arm][1] = length * int(model.Parameters["multiplier"])
		totallen += model.ChromosomeArms[chromosome][arm][1]
	}

	model.FreeParameters["numbits"] = totallen

	return nil
}
