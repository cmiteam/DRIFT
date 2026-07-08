package core

import (
	"image"
)

type Model struct {
	Parameters      map[string]float64
	StringParams    map[string]string // non-numeric params, e.g. birth_style=standard
	Scenario        string
	FreeParameters  map[string]int
	PlotFlags       map[string]bool
	ChromosomeArms  map[int]map[int][]int
	DeathRisk       map[int]float64
	CumulativeProb  map[int]float64
	Map             [][]int
	ModelName       string
	MapName         string
	MatingStyle     string
	Username        string // User running the simulation (for GUI mode)
	BaseModelID     string // Selected base model (e.g., "standard", "flood")
	PrevFrequencies map[int]float64
	HetHistory      []float64
	TimeHistory     []int
	ResultsDir      string // Directory for saving results (supports GUI mode)

	// Demographic scenario (roadmap §6c). When a model ships a demography.json
	// (e.g. the canonical out-of-Africa model), the loader attaches a scheduler
	// here; the run loop calls Apply once per year to advance the schedule
	// (per-deme sizes, split events, and the migration matrix). Nil for ordinary
	// runs, in which case the engine keeps its existing single-cap / scalar-
	// migration behavior unchanged.
	DemographyScheduler DemographyScheduler
	// DemeCaps is the per-deme census target for the CURRENT year, written by the
	// scheduler's Apply. When non-empty, Death culls each deme to its cap instead
	// of applying the global max_pop_size / max_growth_rate limits. Nil/empty when
	// no scenario is active (strict no-op: existing runs are byte-identical).
	DemeCaps map[int]int
}

// DemographyScheduler advances a time-varying demographic scenario. Defined here
// (rather than in pkg/demography) so core.Model can hold one without importing the
// implementation — pkg/demography imports core, not the reverse. Apply is called
// once per simulated year, after model.FreeParameters["year"] is set and before
// Birth: it performs any split due this year, migrates individuals per the current
// epoch's migration matrix, and writes the per-deme caps into model.DemeCaps.
type DemographyScheduler interface {
	Apply(model *Model, pop *Pop)
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

// StringParam returns a string-valued parameter, or def if unset/empty.
func (m *Model) StringParam(key, def string) string {
	if m.StringParams != nil {
		if v, ok := m.StringParams[key]; ok && v != "" {
			return v
		}
	}
	return def
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
	Deme                             // deme (subpopulation) membership; 0 = single-population default
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
	Class     int     // Mutation-class index (roadmap §1): 0 = point, 1.. per the
	// configured mutation_classes spectrum. Lets downstream stats filter/partition
	// mutations by class (point / indel / CNV / large deletion). Default 0 (point),
	// which is also the gob zero value, so old checkpoints decode unchanged.
	Size int // Target size in genome bits recorded for this mutation's class
	// (1 for point substitutions; larger for indel / CNV / large-deletion classes).
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
