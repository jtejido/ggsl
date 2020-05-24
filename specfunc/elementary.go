package specfunc

import (
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"math"
)

func Multiply_e(x, y float64, result *Result) err.GSLError {
	ax := math.Abs(x)
	ay := math.Abs(y)

	if x == 0.0 || y == 0.0 {
		/* It is necessary to eliminate this immediately.
		 */
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if (ax <= 1.0 && ay >= 1.0) || (ay <= 1.0 && ax >= 1.0) {
		/* Straddling 1.0 is always safe.
		 */
		result.val = x * y
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		f := 1.0 - 2.0*gsl.Float64Eps
		min := gsl.Min(math.Abs(x), math.Abs(y))
		max := gsl.Max(math.Abs(x), math.Abs(y))
		if max < 0.9*gsl.SqrtMaxFloat64 || min < (f*math.MaxFloat64)/max {
			result.val = x * y
			result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
			CheckUnderflow(result)
			return nil
		} else {
			return OverflowError(result)
		}
	}
}

func Multiply_err_e(x, dx, y, dy float64, result *Result) err.GSLError {
	status := Multiply_e(x, y, result)
	result.err += math.Abs(dx*y) + math.Abs(dy*x)
	return status
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

func Multiply(x, y float64) float64 {
	result := new(Result)
	status := Multiply_e(x, y, result)
	return EvalResult(result, status)
}
