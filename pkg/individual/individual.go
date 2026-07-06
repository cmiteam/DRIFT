package individual

const EmptyField = -1

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

// MakeIndData creates an individual data array with all fields initialized to EmptyField

func MakeIndData() []int {
	data := make([]int, IndDataFieldCount)
	// Initialize fields that need non-zero defaults
	data[MarriageState] = -1
	data[Dad] = -1
	data[Mom] = -1
	return data
}

func GetStringField(field IndData) string {
	switch field {
	case Sex:
		return "Sex"
	case Dad:
		return "Dad"
	case Mom:
		return "Mom"
	case BirthYear:
		return "BirthYear"
	case MarriageState:
		return "MarriageState"
	case Lifespan:
		return "Lifespan"
	case Lat:
		return "Lat"
	case Lon:
		return "Lon"
	case YGens:
		return "YGens"
	case MtGens:
		return "MtGens"
	case MinGenealoGens:
		return "MinGenealoGens"
	case MaxGenealoGens:
		return "MaxGenealoGens"
	case AlleleCount:
		return "AlleleCount"
	case NumBlocks:
		return "NumBlocks"
	case NumCentromeres:
		return "NumCentromeres"
	case Fitness:
		return "Fitness"
	case NumMutations:
		return "NumMutations"
	case NumBirths:
		return "NumBirths"
	case LastBirthYear:
		return "LastBirthYear"
	case Sons:
		return "Sons"
	case Daughters:
		return "Daughters"
	case Deme:
		return "Deme"
	default:
		return "UnknownField"
	}
}
