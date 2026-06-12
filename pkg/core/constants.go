package core

// Mating style constants
// These correspond to the mating_style parameter in parameter_defaults.csv
const (
	MatingRandom      = 0 // Random mating without distance considerations
	MatingDistance    = 1 // Distance-based mating with geographic constraints
	MatingAgeDistance = 2 // Combined age and distance-based mating
)

// Seeding style constants
// These correspond to the seed_style parameter in parameter_defaults.csv
const (
	SeedingSingle     = 0 // Seed from a single individual
	SeedingPopulation = 1 // Seed from entire population
	SeedingMaxHet     = 2 // Seed to maximize heterozygosity
)
