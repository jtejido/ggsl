package ggsl

import (
	"math"
)

const (
	E           = 2.71828182845904523536028747135
	Ln2E        = 1.44269504088896340735992468100 /* log_2 (e) */
	Ln10E       = 0.43429448190325182765112891892 /* log_10 (e) */
	Sqrt2       = 1.41421356237309504880168872421
	SqrtHalf    = 0.70710678118654752440084436210 // sqrt(1/2)
	Sqrt3       = 1.73205080756887729352744634151 /* sqrt(3) */
	Pi          = 3.14159265358979323846264338328
	TwoPi       = 6.28318530717958647692528676656 // 2 * pi
	PiOver2     = 1.57079632679489661923132169164 // pi/2
	PiOver4     = 0.78539816339744830961566084582 // pi/4
	SqrtPi      = 1.77245385090551602729816748334
	OneOverPi   = 0.31830988618379067153776752675 // 1/pi
	TwoOverPi   = 0.63661977236758134307553505349 // 2/pi
	TwoOverSqPi = 1.12837916709551257389615890312
	Ln10        = 2.30258509299404568401799145468
	Ln2         = 0.69314718055994530941723212146
	Ln3         = 1.09861228866810969139524523692
	LnPi        = 1.14472988584940017414342735135 // log(pi)
	Euler       = 0.57721566490153286060651209008
	PiSq        = 9.86960440108935799230494012590
	Apery       = 1.2020569031595942853997381615114499907649 //OEIS: A002117
	Catalan     = 0.91596559417721901505460351493
)

// f(a) = x
type Function interface {
	Evaluate(x float64) float64
}

func Evaluate(f Function, x float64) float64 {
	return f.Evaluate(x)
}

func Sign(x float64) float64 {
	if (x) >= 0.0 {
		return 1.
	}

	return -1.
}

func IsOdd(n int) bool {
	return n&1 == 1
}

func IsEven(n int) bool {
	return !IsOdd(n)
}

func IsReal(x float64) bool {
	return !math.IsNaN(x) && !math.IsInf(x, 0)
}

func AbsInt(n int) int {
	if n < 0 {
		return -n
	}
	return n
}

func Max(in ...interface{}) float64 {
	var max_i int
	flin := make([]float64, len(in))
	for i, v := range in {
		switch v.(type) {
		case int:
			flin[i] = float64(v.(int))
		case float64:
			flin[i] = v.(float64)
		}
		if flin[i] > flin[max_i] {
			max_i = i
		}
	}
	return flin[max_i]
}

func Min(in ...interface{}) float64 {
	var max_i int
	flin := make([]float64, len(in))
	for i, v := range in {
		switch v.(type) {
		case int:
			flin[i] = float64(v.(int))
		case float64:
			flin[i] = v.(float64)
		}
		if flin[i] < flin[max_i] {
			max_i = i
		}
	}
	return flin[max_i]
}
