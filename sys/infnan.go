package sys

import (
	"math"
)

// isNaN reports whether f is an IEEE 754 ``not-a-number'' value.
func isNaN(f float64) (is bool) {
	return math.IsNaN(f)
}

// isFinite reports whether f is neither NaN nor an infinity.
func isFinite(f float64) bool {
	return !isNaN(f) && !isInf(f)
}

// isInf reports whether f is an infinity.
func isInf(f float64) bool {
	return math.IsInf(f, 1) || math.IsInf(f, -1)
}

func NaN() float64 {
	return Fdiv(0.0, 0.0)
}

func Posinf() float64 {
	return Fdiv(+1.0, 0.0)
}

func Neginf() float64 {
	return Fdiv(-1.0, 0.0)
}

func IsNaN(x float64) int {
	if isNaN(x) {
		return 1
	}

	return 0
}

func IsInf(x float64) int {

	if math.IsInf(x, 1) {
		return 1
	} else if math.IsInf(x, -1) {
		return -1
	}

	return 0
}

func Finite(x float64) int {
	if isFinite(x) {
		return 1
	}

	return 0
}
