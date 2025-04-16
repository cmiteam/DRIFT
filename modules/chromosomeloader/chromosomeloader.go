package chromosomeloader

import (
	"drift/types"
	"encoding/csv"
	"os"
	"strconv"
)

func LoadChromosomeArmsFromCSV(model *types.Model) {

	file, err := os.Open("C:/Go/Programs/Drift/static/chromosome_data.csv")
	if err != nil {
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
	}

	var totallen int

	for _, record := range records[1:] {
		chromosome, err := strconv.Atoi(record[0])
		if err != nil {
		}

		arm, err := strconv.Atoi(record[1])
		if err != nil {
		}

		start, err := strconv.Atoi(record[2])
		if err != nil {
		}

		length, err := strconv.Atoi(record[3])
		if err != nil {
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

}
