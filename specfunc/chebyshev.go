package specfunc

import (
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"math"
)

type chebyshevSeries struct {
	c        []float64 // coefficients
	order    int       // degree of expansion
	a, b     float64   // lower and upper interval
	order_sp int       /* effective single precision order */
}

func (cs *chebyshevSeries) evaluate(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	d := 0.0
	dd := 0.0

	y := (2.0*x - cs.a - cs.b) / (cs.b - cs.a)
	y2 := 2.0 * y

	e := 0.0

	eval_order := cs.order

	if mode != gsl.PREC_DOUBLE {
		eval_order = cs.order_sp
	}

	for j := eval_order; j >= 1; j-- {
		temp := d
		d = y2*d - dd + cs.c[j]
		e += math.Abs(y2*temp) + math.Abs(dd) + math.Abs(cs.c[j])
		dd = temp
	}

	{
		temp := d
		d = y*d - dd + 0.5*cs.c[0]
		e += math.Abs(y*temp) + math.Abs(dd) + 0.5*math.Abs(cs.c[0])
	}

	result.val = d
	result.err = gsl.Float64Eps*e + math.Abs(cs.c[eval_order])

	return nil
}

func (cs *chebyshevSeries) Evaluate(x float64, result *Result) err.GSLError {
	return cs.evaluate(x, gsl.PREC_DOUBLE, result)
}

func (cs *chebyshevSeries) EvaluateMode(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	return cs.evaluate(x, mode, result)
}
