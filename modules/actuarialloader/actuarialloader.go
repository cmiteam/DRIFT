package actuarialloader

import (
	"drift/types"
	"encoding/csv"
	"os"
	"strconv"
)

func LoadActuarialTable(model *types.Model) {

	file, err := os.Open("C:/Go/Programs/Drift/static/actuarial_table.csv")
	if err != nil {
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
	}

	var cumulative float64

	for _, row := range records[1:] { // Skip the header row

		age, err := strconv.Atoi(row[0])
		if err != nil {
		}

		risk, err := strconv.ParseFloat(row[1], 64)
		if err != nil {
		}
		model.DeathRisk[age] = risk

		popProb, err := strconv.ParseFloat(row[2], 64)
		if err != nil {
		}
		cumulative += popProb * 4.4
		model.CumulativeProb[age] = cumulative
	}

}
