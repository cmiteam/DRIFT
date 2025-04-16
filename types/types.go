package types

type Model struct {
	Parameters     map[string]float64
	FreeParameters map[string]int
	PlotFlags      map[string]bool
	ChromosomeArms map[int]map[int][]int
	DeathRisk      map[int]float64
	CumulativeProb map[int]float64
	ModelName      string
}

type Pop struct {
	IndData       map[int]map[string]int // Individual data
	Chromosomes   map[int][][]uint64     // Genetic data
	Centromeres   map[int][]uint64       // Centromere information
	IndMutations  map[int]map[int][]int  // Mutations per individual
	MutationPool  map[int]Mutation       // Global pool of mutations
	MutationHist  map[int]int            // Mutation history/statistics
	MutationCount int
	Tracking      map[string]int
}

type Mutation struct {
	Id        int     // Unique mutation identifier
	Position  int     // Position in genome (base-pair level)
	Effect    float64 // Mutation effect value
	Origin    int     // Original individual
	Dominance int     // 0, 100, or somewhere in between
	Count     int     // Number of instances in circulation
}

type Chromosomes struct {
	data [][]uint64
}
