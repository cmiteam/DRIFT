package core

import (
	"image"
)

type Model struct {
	Parameters      map[string]float64
	FreeParameters  map[string]int
	PlotFlags       map[string]bool
	ChromosomeArms  map[int]map[int][]int
	DeathRisk       map[int]float64
	CumulativeProb  map[int]float64
	Map             [][]int
	ModelName       string
	MapName         string
	MatingStyle     string
	PrevFrequencies map[int]float64
	HetHistory      []float64
	TimeHistory     []int
}

type Pop struct {
	IndData       map[int][]int         // Individual data
	Chromosomes   map[int][][]uint64    // Genetic data
	Centromeres   map[int][2][]uint64   // Centromere information
	IndMutations  map[int]map[int][]int // Mutations per individual
	MutationPool  map[int]Mutation      // Global pool of mutations
	MutationHist  map[int]int           // Mutation history/statistics
	MutationCount int
	Tracking      map[string]int
	AlleleFreqs   map[int][]int16
	SFSHistory    map[int][]int
	MaleDB        map[int]Ancestor
	FemaleDB      map[int]Ancestor
}

type IndData int
const (
	Sex               IndData = iota // 0 = male, 1 = female
	Dad                              // ID of dad
	Mom                              // ID of mom
	BirthYear                        // year of birth
	MarriageState                    // ID of spouse
	Lifespan                         // potential lifespan
	Lat                              // latitude
	Lon                              // longitude
	YGens                            // generations from male seed
	MtGens                           // generations from female seed
	MinGenealoGens                   // shortest path on family tree to seed
	MaxGenealoGens                   // longest path on family tree to seed
	AlleleCount                      // tracking descent from seed individual(s)
	NumBlocks                        // number of contiguous blocks of set bits
	NumCentromeres                   // number of centromeres
	Fitness                          // used for survival calculations
	NumMutations                     // number of mutations
	NumBirths                        // number of children
	LastBirthYear                    // year of last birth
	Sons                             // number of sons born to this individual
	Daughters                        // number of daughters born to this individual
	IndDataFieldCount                // number of fields in IndData
)
const EmptyField = -999

type Mutation struct {
	Id        int     // Unique mutation identifier
	Position  int     // Position in genome (base-pair level)
	Effect    float64 // Mutation effect value
	Origin    int     // Original individual
	Dominance int     // 0, 100, or somewhere in between
	Count     int     // Number of instances in circulation
}

type AlleleCounts struct {
	counts  []int16 // allele count * 10000
	popSize int16   // population size * 2
}

// AnimationsContainer manages multiple animations with a shared base map
type AnimationsContainer struct {
	BaseMap    *image.RGBA
	Collection map[string]*AnimationWriter
	TileSize   int // Calculated once during initialization
}

// AnimationWriter tracks and creates a single animation
type AnimationWriter struct {
	Frames []*image.Paletted
	Delays []int
	Width  int
	Height int
}

type DriftStats struct {
	NumSegSites  int     // Allelic richness (number of segregating sites)
	AvgHet       float64 // Average heterozygosity
	AllelesLost  int     // Sites that went to 0% this generation
	AllelesFixed int     // Sites that went to 100% this generation
	AvgFreq      float64 // Average allele frequency
	FreqVariance float64 // Variance in allele frequencies
	EffectiveNe  float64 // Effective population size estimate
}

// SFSResult contains the site frequency spectrum data and derived statistics
type SFSResult struct {
	FoldedSFS   []int   // Frequency bins for folded SFS (0 to n/2)
	UnfoldedSFS []int   // Frequency bins for unfolded SFS (0 to n) - if ancestral state known
	TotalSites  int     // Total number of segregating sites
	Singletons  int     // Number of singleton variants
	Doubletons  int     // Number of doubleton variants
	TajimasD    float64 // Tajima's D statistic
	FuLiD       float64 // Fu & Li's D* statistic
	FayWuH      float64 // Fay & Wu's H statistic
	ThetaW      float64 // Watterson's theta (from number of segregating sites)
	ThetaPi     float64 // Tajima's pi (from pairwise differences)
	NumSamples  int     // Number of chromosomes sampled
}

// LDResult contains linkage disequilibrium statistics between two sites
type LDResult struct {
	Site1    int     // Position of first site
	Site2    int     // Position of second site
	Distance int     // Physical distance between sites
	RSquared float64 // r-squared - correlation coefficient
	DPrime   float64 // D' - normalized LD coefficient
	ChiSq    float64 // Chi-square test statistic
	PValue   float64 // P-value for independence test
}

// CoalescenceResult stores results from pairwise coalescence analysis
type CoalescenceResult struct {
	Individual1     int  // ID of first individual
	Individual2     int  // ID of second individual
	CoalescenceTime int  // Generations to MRCA
	CommonAncestor  int  // ID of MRCA (if found)
	PathLength1     int  // Generations from ind1 to MRCA
	PathLength2     int  // Generations from ind2 to MRCA
	FoundMRCA       bool // Whether MRCA was found in tracked genealogy
}

type TMRCAResult struct {
	MeanTMRCA     float64 // Average coalescence time across samples
	MaxTMRCA      float64 // Maximum coalescence time found
	VarianceTMRCA float64 // Variance in coalescence times
	SampleSize    int     // Number of successful coalescence events
	EstimatedNe   float64 // Rough estimate of effective population size
	ActualPopSize int     // Current census population size
}

// YAdamResult stores Y-chromosome Adam analysis results
type YAdamResult struct {
	YAdamID         int  // ID of the Y-chromosome Adam (if found)
	GenerationsBack int  // Generations back to Y-Adam
	YAdamBirthYear  int  // Birth year of Y-Adam
	NumMales        int  // Number of males analyzed
	FoundYAdam      bool // Whether Y-Adam was identified
}

type MtEveResult struct {
	MtEveID         int  // ID of mitochondrial Eve (if found)
	GenerationsBack int  // Generations back to Mt-Eve
	MtEveBirthYear  int  // Birth year of Mt-Eve
	FoundMtEve      bool // Whether Mt-Eve was identified
}

type Ancestor struct {
	ID        int
	BirthYear int
}
