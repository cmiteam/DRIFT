package utils

import (
	"math"
)

func Min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func Max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func Abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func WeibullRandom(shape, scale float64) float64 {
	u := RandFloat64()
	return scale * math.Pow(-math.Log(u), 1/shape)
}
