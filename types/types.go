package types

import (
	"image"
)

type Model struct {
	Parameters     map[string]float64
	FreeParameters map[string]int
	PlotFlags      map[string]bool
	ChromosomeArms map[int]map[int][]int
	DeathRisk      map[int]float64
	CumulativeProb map[int]float64
	Map            map[float64]map[float64]int
	ModelName      string
	MapName        string
	MatingStyle    string
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
}

type Mutation struct {
	Id        int     // Unique mutation identifier
	Position  int     // Position in genome (base-pair level)
	Effect    float64 // Mutation effect value
	Origin    int     // Original individual
	Dominance int     // 0, 100, or somewhere in between
	Count     int     // Number of instances in circulation
}

// AnimationsContainer manages multiple animations with a shared base map
type AnimationsContainer struct {
	BaseMap     *image.RGBA
	Collection  map[string]*AnimationWriter
	MinLat      float64
	MaxLat      float64
	MinLon      float64
	MaxLon      float64
	MapWidth    int
	MapHeight   int
	latScale    float64
	lonScale    float64
	ScaleFactor float64
}

// AnimationWriter tracks and creates a single animation
type AnimationWriter struct {
	Frames []*image.Paletted
	Delays []int
	Width  int
	Height int
}
