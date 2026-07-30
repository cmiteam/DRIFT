package core

import (
	"image"
)

type Model struct {
	Parameters     map[string]float64
	StringParams   map[string]string // non-numeric params, e.g. birth_style=standard
	Scenario       string
	FreeParameters map[string]int
	PlotFlags      map[string]bool
	ChromosomeArms map[int]map[int][]int
	// ArmCM holds an optional per-arm genetic length in centiMorgans (roadmap §6a),
	// parsed from a 5th column of chromosome_data.csv (Chromosome,Arm,Start,Length,cM).
	// Empty when the file has no cM column (legacy 4-column files) — the map then
	// falls back to a uniform recomb_cM_per_bit rate. Keyed [chrom][arm].
	ArmCM map[int]map[int]float64
	// GeneticMap is the lazily-built recombination map (roadmap §6a) derived from
	// ChromosomeArms + ArmCM. Nil until first needed and only built when a map is
	// actually used (a cM column is present, or recombination_model="map"); a nil map
	// means the legacy single-interior-segment createMask is in effect and the §6d LD
	// leg falls back to the ne_cM_per_bit scalar. Built once (no RNG), then cached.
	GeneticMap      *GeneticMap
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

	// Scriptable environmental events (roadmap §2). When a model sets the
	// `environmental_events` param the loader parses it into a schedule and
	// attaches it here; the run loop calls Apply once per year, after the
	// demographic scheduler and before Seed/Birth, to fire famines, migration
	// pulses, etc. Nil for runs with no scheduled events (strict no-op).
	EventScheduler EventScheduler
	// EventMortalityFactor is a transient per-year multiplier on actuarial death
	// risk, rewritten by the event scheduler every year (a famine raises it above
	// 1 during its window; it is reset to 1 when no mortality event is active).
	// Read by Death only when EventScheduler is non-nil, so runs without events
	// take the byte-identical original death path.
	EventMortalityFactor float64

	// Habitat suitability table (roadmap §3). Lazily built (see ensureHabitat) from
	// the `habitat_suitability` param, which maps terrain codes to a per-cell
	// suitability multiplier. Nil until first consulted; once built, Active reports
	// whether any multiplier was configured. When inactive the cell accessors fall
	// back to today's binary land/water semantics, so an unset param leaves every
	// existing path byte-identical. Rebuilt from the param, not serialized (like the
	// schedulers) — no checkpoint schema bump.
	Habitat *HabitatTable

	// Movement barriers table (roadmap §3). Lazily built (see ensureBarriers) from
	// the `movement_barriers` param, which lists terrain codes that block *crossing*
	// (as opposed to habitat_suitability, which governs where an individual may live).
	// Nil until first consulted; when inactive (unset param) CanTraverse returns true
	// unconditionally, so Wander's move-accept predicate is byte-identical. Rebuilt
	// from the param, not serialized (like the schedulers) — no checkpoint schema bump.
	Barriers *BarrierTable

	// Deme migration matrix (roadmap §3). Lazily built (see ensureMigration) from the
	// `migration_matrix` param, which lists directional per-pair per-year migration
	// rates between demes (an asymmetric generalization of the scalar deme_migration_rate
	// island model). Nil until first consulted; when inactive (unset param) the island
	// mating path uses the scalar deme_migration_rate unchanged, so an unset param leaves
	// every existing path byte-identical. Rebuilt from the param, not serialized (like the
	// schedulers / habitat / barriers) — no checkpoint schema bump.
	Migration *MigrationMatrix
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

// EventScheduler advances scriptable environmental events (roadmap §2) — famine,
// migration pulse, and so on. Defined here (like DemographyScheduler) so core.Model
// can hold one without importing pkg/events (pkg/events imports core, not the
// reverse). Apply is called once per simulated year, after
// model.FreeParameters["year"] is set and after the demographic scheduler, before
// Seed/Birth: it resets the per-year event modifiers, then fires any events due
// this year (mortality windows, deme relabels). Implementations consume RNG only
// in sorted-id order — and none at all in years no event acts — so a fixed
// rng_seed stays byte-reproducible and years before the first event match a
// no-event baseline exactly.
type EventScheduler interface {
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

	// Ne-through-time state (roadmap §6d). NeHistory holds one NeSnapshot per
	// captured window, keyed by simulated year; SaveNeTimeSeries writes it at
	// end-of-run. NePrevSnap is the previous window's allele-frequency sample,
	// held transiently so the temporal (Waples) estimator can compare successive
	// windows. Both are populated only when track_Ne is set (default off ⇒ nil ⇒
	// no capture, no RNG, byte-identical). Exported so they survive gob
	// checkpointing; nil is a valid zero value on resume.
	NeHistory  map[int]*NeSnapshot
	NePrevSnap *FreqSnapshot

	// Genetic-load-through-time state (roadmap §6f). One LoadSnapshot per captured
	// window, keyed by simulated year; SaveLoadTimeSeries writes it at end-of-run.
	// Populated only when track_load is set (default off ⇒ nil ⇒ no capture, no
	// RNG, byte-identical). A map so gob checkpointing tolerates it; nil is a valid
	// zero value on resume.
	LoadHistory map[int]*LoadSnapshot
}

// LoadSnapshot holds the genetic-load statistics for the living population in one
// captured window (roadmap §6f: DFE & genetic load / mutational-meltdown). It is a
// READ-ONLY characterization of the de-novo mutation pool — a deterministic full-
// population scan that consumes no RNG, so a track_load run is byte-identical to a
// track_load-off run. Exported fields so it survives gob checkpointing.
//
// Load decomposition: TotalLoad = 1 - MeanFitness is the realized load under the
// active fitness_model (so it reflects dominance + epistasis exactly). FixedLoad
// is the load an individual carrying ONLY the fixed (frequency-1) derived variants
// would bear (identical for everyone, since fixed lineages are homozygous in all);
// SegLoad = TotalLoad - FixedLoad is the residual from segregating variation
// (exact under additive fitness; an approximation under multiplicative/synergistic,
// where loci interact). NumFixedDel is the Muller's-ratchet counter.
type LoadSnapshot struct {
	Year          int
	CensusN       int
	MeanFitness   float64         // mean relative fitness over living (per fitness_model)
	TotalLoad     float64         // 1 - MeanFitness
	SegLoad       float64         // load attributable to segregating (0<f<1) lineages
	FixedLoad     float64         // load from fixed (f==1) derived lineages
	NumSeg        int             // distinct segregating lineages among living
	NumFixed      int             // distinct fixed lineages among living
	NumFixedDel   int             // fixed deleterious lineages (Muller's-ratchet counter)
	NumFixedBen   int             // fixed beneficial lineages
	MeanMutPerInd float64         // mean mutation copies carried per individual
	MeanDelPerInd float64         // mean deleterious mutation copies per individual
	DemeLoad      map[int]float64 // per-deme total load (nil when a single deme)
	ClassLoad     map[int]float64 // per mutation-class additive load (nil when only class 0)
}

// FreqSnapshot is one temporal sample of allele frequencies at segregating sites,
// grouped by population unit (group key -1 = whole population pooled; key >= 0 = a
// deme id). Held transiently on Pop to feed the temporal Ne estimator (§6d): each
// window is compared against the previous snapshot. Frequencies are derived-allele
// frequencies over the Chromosomes bitfield. Exported fields so it survives gob
// checkpointing.
type FreqSnapshot struct {
	Year     int
	Groups   map[int]map[int]float64 // group -> genome-bit site -> derived-allele frequency
	HaploidN map[int]int             // group -> haploid sample size (2 * sampled diploids)
}

// NeSnapshot holds the effective-population-size estimates for one captured window
// (§6d). TemporalNe is the pooled temporal (Waples 1989 plan-II) estimate over the
// interval since the previous window; DemeNe carries the same per deme when the
// island model is active. LDNe is the caveated LD-r² recent-Ne companion (0 when
// disabled or unresolved; approximate pending a real recombination map, §6a). A
// non-finite Ne means drift was below sampling noise (unresolved) that window.
type NeSnapshot struct {
	Year             int
	CensusN          int             // living individuals at this window
	IntervalGen      float64         // t: generations since the previous window
	TemporalNe       float64         // pooled temporal Ne over the interval
	TemporalNumSites int             // segregating sites used for the pooled estimate
	LDNe             float64         // LD-based recent Ne (approx; 0 if disabled/unresolved)
	LDNumPairs       int             // site pairs used for the LD estimate
	DemeNe           map[int]float64 // per-deme temporal Ne (nil when single deme)
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
