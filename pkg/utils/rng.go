package utils

import (
	"encoding/binary"
	"fmt"
	"math"
	"math/rand"
)

// simRNG is the single source of randomness for a simulation run. It is seeded
// once per run via SeedRNG so that runs are reproducible.
//
// Simulation runs execute as separate OS processes — a direct CLI invocation, or
// one exec.Command per web-queue job (see pkg/webserver/queue.go) — so no two
// runs share this generator within a process. A package-level, unsynchronized
// generator is therefore safe here and avoids the deprecated global rand.Seed.
//
// The source is a small custom generator (xoshiro256** seeded via splitmix64)
// rather than math/rand's default source, because the default source is not
// serializable in this Go version. A capturable state is required for
// checkpoint/resume (see RNGState / RestoreRNG and pkg/checkpoint).
var (
	simSource = newXoshiro(1)
	simRNG    = rand.New(simSource)
	currentSeed int64 = 1
)

// SeedRNG (re)initializes the shared generator for a run. Call once at the start
// of each run. Passing the same seed reproduces the run byte-for-byte.
func SeedRNG(seed int64) {
	currentSeed = seed
	simSource.seed(seed)
	simRNG = rand.New(simSource)
}

// CurrentSeed returns the seed the shared generator was last initialized with.
func CurrentSeed() int64 { return currentSeed }

// RNGState returns a snapshot of the generator's exact position, for resuming a
// checkpoint on the identical random stream. Pair with RestoreRNG.
func RNGState() []byte {
	b := make([]byte, 32)
	for i := 0; i < 4; i++ {
		binary.LittleEndian.PutUint64(b[i*8:], simSource.s[i])
	}
	return b
}

// RestoreRNG restores generator state previously captured by RNGState.
func RestoreRNG(state []byte) error {
	if len(state) != 32 {
		return fmt.Errorf("invalid RNG state length %d (want 32)", len(state))
	}
	for i := 0; i < 4; i++ {
		simSource.s[i] = binary.LittleEndian.Uint64(state[i*8:])
	}
	simRNG = rand.New(simSource)
	return nil
}

// Rand exposes the shared generator for callers that need a *rand.Rand directly.
func Rand() *rand.Rand { return simRNG }

// RandIntn returns a pseudo-random int in [0, n) from the shared generator.
func RandIntn(n int) int { return simRNG.Intn(n) }

// RandFloat64 returns a pseudo-random float64 in [0.0, 1.0) from the shared generator.
func RandFloat64() float64 { return simRNG.Float64() }

// RandShuffle pseudo-randomizes the order of n elements using the shared generator.
func RandShuffle(n int, swap func(i, j int)) { simRNG.Shuffle(n, swap) }

// RandPoisson draws from a Poisson distribution with the given mean using the
// shared generator. It replaces the gonum distuv.Poisson draw so that mutation
// counts are reproducible under a fixed seed. Knuth's algorithm is used for
// small means; a normal approximation is used for large means to avoid exp()
// underflow and long loops.
func RandPoisson(lambda float64) int {
	if lambda <= 0 {
		return 0
	}
	if lambda > 30 {
		v := int(math.Round(lambda + math.Sqrt(lambda)*simRNG.NormFloat64()))
		if v < 0 {
			v = 0
		}
		return v
	}
	l := math.Exp(-lambda)
	k := 0
	p := 1.0
	for {
		k++
		p *= simRNG.Float64()
		if p <= l {
			return k - 1
		}
	}
}

// xoshiro256** — a fast, high-quality PRNG with 256 bits of state. Its state is
// four uint64s, which makes it trivially serializable (unlike math/rand's default
// source). Implements math/rand.Source64.
type xoshiro struct {
	s [4]uint64
}

func newXoshiro(seed int64) *xoshiro {
	x := &xoshiro{}
	x.seed(seed)
	return x
}

// seed initializes the 256-bit state from a 64-bit seed using splitmix64, as
// recommended by the xoshiro authors.
func (x *xoshiro) seed(seed int64) {
	sm := uint64(seed)
	for i := 0; i < 4; i++ {
		sm += 0x9E3779B97F4A7C15
		z := sm
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
		z = (z ^ (z >> 27)) * 0x94D049BB133111EB
		x.s[i] = z ^ (z >> 31)
	}
}

func rotl(x uint64, k int) uint64 { return (x << k) | (x >> (64 - k)) }

func (x *xoshiro) Uint64() uint64 {
	result := rotl(x.s[1]*5, 7) * 9
	t := x.s[1] << 17
	x.s[2] ^= x.s[0]
	x.s[3] ^= x.s[1]
	x.s[1] ^= x.s[2]
	x.s[0] ^= x.s[3]
	x.s[2] ^= t
	x.s[3] = rotl(x.s[3], 45)
	return result
}

func (x *xoshiro) Int63() int64 { return int64(x.Uint64() >> 1) }

// Seed satisfies rand.Source; SeedRNG is the intended entry point.
func (x *xoshiro) Seed(seed int64) { x.seed(seed) }
