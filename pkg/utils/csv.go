package utils

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
)

// Load CSV files
type CSVLoader struct {
	FileName   string
	Dir        string
	MinRecords int
}

// Error type for invalid records.
type ErrInvalidRecord struct {
	CSVLoader CSVLoader
	Record    []string
	Message   string
}

// Error type for invalid fields.
type ErrInvalidField struct {
	CSVLoader CSVLoader
	Record    []string
	Field     int
	Message   string
}

// Attempt to load a CSV file and return its contents as a slice of slices of strings.
func (csvLoader CSVLoader) LoadCSV() ([][]string, error) {
	path := filepath.Join(csvLoader.Dir, csvLoader.FileName)
	file, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("Error loading %s: %v", path, err)
	}
	defer file.Close()

	// Read all records from the file.
	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("Error reading records from %s: %v", path, err)
	}

	// Ensure the CSV file has enough records
	if len(records) < csvLoader.MinRecords {
		return nil, fmt.Errorf("Error: %s has fewer than %d records", path, csvLoader.MinRecords)
	}

	// yay
	return records, nil
}

// Check a record to ensure it has at least a certain number of fields.
func (csvLoader *CSVLoader) CheckRecord(record []string, minFields int) error {
	if len(record) < minFields {
		return ErrInvalidRecord{
			CSVLoader: *csvLoader,
			Record:    record,
			Message:   fmt.Sprintf("Less than %d fields", minFields),
		}
	}
	return nil
}

// Convert a field in a record to an integer.
// You could do it yourself but this will give you a nice error message
// including the file name, record, and field index.
func (csvLoader *CSVLoader) Atoi(record []string, field int) (int, error) {
	value, err := strconv.Atoi(record[field])
	if err != nil {
		return 0, ErrInvalidField{
			CSVLoader: *csvLoader,
			Record:    record,
			Field:     field,
			Message:   fmt.Sprintf("%v is not an integer", record[field]),
		}
	}
	return value, nil
}

// Convert a field in a record to an float64.
// You could do it yourself but this will give you a nice error message
// including the file name, record, and field index.
func (csvLoader *CSVLoader) ParseFloat64(record []string, field int) (float64, error) {
	value, err := strconv.ParseFloat(record[field], 64)
	if err != nil {
		return 0, ErrInvalidField{
			CSVLoader: *csvLoader,
			Record:    record,
			Field:     field,
			Message:   fmt.Sprintf("%v is not a float", record[field]),
		}
	}
	return value, nil
}

// Convert a field in a record to a boolean value.
// You could do it yourself but this will give you a nice error message
// including the file name, record, and field index.
func (csvLoader *CSVLoader) ParseBool(record []string, field int) (bool, error) {
	value, err := strconv.ParseBool(record[field])
	if err != nil {
		return false, ErrInvalidField{
			CSVLoader: *csvLoader,
			Record:    record,
			Field:     field,
			Message:   fmt.Sprintf("%v is not a boolean", record[field]),
		}
	}
	return value, nil
}

// Error() method for ErrInvalidRecord.
func (e ErrInvalidRecord) Error() string {
	path := filepath.Join(e.CSVLoader.Dir, e.CSVLoader.FileName)
	return fmt.Sprintf("Invalid record in %v: %v. %s", path, e.Record, e.Message)
}

// Error() method for ErrInvalidField.
func (e ErrInvalidField) Error() string {
	path := filepath.Join(e.CSVLoader.Dir, e.CSVLoader.FileName)
	return fmt.Sprintf("Invalid field in %v: %v. At index %d: %s", path, e.Record, e.Field, e.Message)
}
