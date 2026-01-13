package config

import (
	"drift/pkg/core"
	"drift/pkg/simulation"
	"drift/pkg/utils"
	"fmt"
	"math"
)

// validateBasicParameters checks that critical parameters exist and are within valid ranges
func validateBasicParameters(params map[string]float64) error {
	// Check required parameters exist and are in valid ranges
	requiredParams := map[string]struct {
		min float64
		max float64
	}{
		"start_pop_size":  {1, 1e6},
		"max_pop_size":    {1, 1e6},
		"maturity":        {0, 100},
		"lifespan":        {1, 1000},
		"min_lifespan":    {1, 1000},
		"end_year":        {1, 1e6},
		"max_growth_rate": {0, 10},
		"seed_year":       {0, 1e6},
		"spacing":         {1, 100},
		"birth_prob":      {0, 100},
		"menopause":       {0, 1},
	}

	// Validate each required parameter
	for param, bounds := range requiredParams {
		value, exists := params[param]
		if !exists {
			return fmt.Errorf("required parameter '%s' is missing", param)
		}
		if value < bounds.min || value > bounds.max {
			return fmt.Errorf("parameter '%s' value %.2f is out of valid range [%.0f, %.0f]",
				param, value, bounds.min, bounds.max)
		}
	}

	// Logical validations - check parameter relationships
	if params["start_pop_size"] > params["max_pop_size"] {
		return fmt.Errorf("start_pop_size (%.0f) cannot exceed max_pop_size (%.0f)",
			params["start_pop_size"], params["max_pop_size"])
	}

	if params["min_lifespan"] > params["lifespan"] {
		return fmt.Errorf("min_lifespan (%.0f) cannot exceed lifespan (%.0f)",
			params["min_lifespan"], params["lifespan"])
	}

	if params["maturity"] >= params["lifespan"] {
		return fmt.Errorf("maturity (%.0f) must be less than lifespan (%.0f)",
			params["maturity"], params["lifespan"])
	}

	if params["seed_year"] > params["end_year"] {
		return fmt.Errorf("seed_year (%.0f) cannot be after end_year (%.0f)",
			params["seed_year"], params["end_year"])
	}

	// Validate bottleneck parameters if enabled
	bottleneckStart := params["bottleneck_start"]
	bottleneckEnd := params["bottleneck_end"]

	// bottleneck_start of -1 means bottleneck is disabled
	if bottleneckStart >= 0 {
		if bottleneckEnd <= bottleneckStart {
			return fmt.Errorf("bottleneck_end (%.0f) must be greater than bottleneck_start (%.0f)",
				bottleneckEnd, bottleneckStart)
		}
		if params["bottleneck_size"] < 2 {
			return fmt.Errorf("bottleneck_size (%.0f) must be at least 2", params["bottleneck_size"])
		}
		if params["bottleneck_size"] > params["max_pop_size"] {
			return fmt.Errorf("bottleneck_size (%.0f) cannot exceed max_pop_size (%.0f)",
				params["bottleneck_size"], params["max_pop_size"])
		}
	}

	// Validate mating style
	matingStyle := params["mating_style"]
	if matingStyle < 0 || matingStyle > 2 {
		return fmt.Errorf("mating_style (%.0f) must be 0 (random), 1 (distance), or 2 (age_distance)",
			matingStyle)
	}

	// Validate seed style
	seedStyle := params["seed_style"]
	if seedStyle < 0 || seedStyle > 2 {
		return fmt.Errorf("seed_style (%.0f) must be 0 (single), 1 (population), or 2 (max_het)",
			seedStyle)
	}

	return nil
}

// Initializes the model based on the configuration files

func InitializeModel(configRoot string) (*core.Model, error) {
	model := &core.Model{
		Parameters:      make(map[string]float64),
		PlotFlags:       make(map[string]bool),
		ChromosomeArms:  make(map[int]map[int][]int),
		DeathRisk:       make(map[int]float64),
		CumulativeProb:  make(map[int]float64),
		FreeParameters:  make(map[string]int),
		Map:             nil, // this will be set later, once the dimensions are known
		PrevFrequencies: make(map[int]float64),
		HetHistory:      make([]float64, 0),
		TimeHistory:     make([]int, 0),
		ResultsDir:      "results", // Default results directory
	}

	// Attempt to load each config file. Failure will be fatal.
	err := utils.LoadParameters(model, configRoot)
	if err != nil {
		return nil, err
	}

	// Validate parameters before proceeding
	err = validateBasicParameters(model.Parameters)
	if err != nil {
		return nil, fmt.Errorf("parameter validation failed: %w", err)
	}

	err = utils.LoadChromosomes(model, configRoot)
	if err != nil {
		return nil, err
	}
	err = utils.LoadActuarialTable(model, configRoot)
	if err != nil {
		return nil, err
	}

	// Calculate derived values
	model.Parameters["mu_sig_figs"] = math.Pow(1, model.Parameters["mu_sig_figs"])

	// Calculate total genome bits from chromosome data
	totalBits := 0
	for chrom := range model.ChromosomeArms {
		if arm0, ok := model.ChromosomeArms[chrom][0]; ok {
			totalBits = max(totalBits, arm0[0]+arm0[1]) // pstart + plen
		}
		if arm1, ok := model.ChromosomeArms[chrom][1]; ok {
			totalBits = max(totalBits, arm1[0]+arm1[1]) // qstart + qlen
		}
	}

	model.FreeParameters["genome_bits"] = totalBits
	// Prepare output files
	simulation.SaveHeaders(model.ModelName, model.ResultsDir)

	// Initialize free parameters
	model.FreeParameters["indID"] = 0         // Starting ID for individuals
	model.FreeParameters["seed"] = -1         // No seed initially
	model.FreeParameters["last_pop_size"] = 0 // Required for growth rate calculations
	model.FreeParameters["mutID"] = 0         // Starting ID for mutations

	return model, nil
}
