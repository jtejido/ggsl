/* specfunc/atanint.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
package specfunc

import (
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"math"
)

var (
	atanint_data = []float64{
		1.91040361296235937512,
		-0.4176351437656746940e-01,
		0.275392550786367434e-02,
		-0.25051809526248881e-03,
		0.2666981285121171e-04,
		-0.311890514107001e-05,
		0.38833853132249e-06,
		-0.5057274584964e-07,
		0.681225282949e-08,
		-0.94212561654e-09,
		0.13307878816e-09,
		-0.1912678075e-10,
		0.278912620e-11,
		-0.41174820e-12,
		0.6142987e-13,
		-0.924929e-14,
		0.140387e-14,
		-0.21460e-15,
		0.3301e-16,
		-0.511e-17,
		0.79e-18,
	}

	atanint_cs = &chebyshevSeries{
		atanint_data,
		20,
		-1, 1,
		10,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*/
func Atanint_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)
	sgn := gsl.Sign(x)

	if ax == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if ax < 0.5*gsl.SqrtFloat64Eps {
		result.val = x
		result.err = 0.0
		return nil
	} else if ax <= 1.0 {
		t := 2.0 * (x*x - 0.5)
		result_c := new(Result)
		atanint_cs.Evaluate(t, result_c)
		result.val = x * result_c.val
		result.err = x * result_c.err
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if ax < 1.0/gsl.SqrtFloat64Eps {
		t := 2.0 * (1.0/(x*x) - 0.5)
		result_c := new(Result)
		atanint_cs.Evaluate(t, result_c)
		result.val = sgn * (0.5*gsl.Pi*math.Log(ax) + result_c.val/ax)
		result.err = result_c.err/ax + math.Abs(result.val*gsl.Float64Eps)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = sgn * (0.5*gsl.Pi*math.Log(ax) + 1.0/ax)
		result.err = 2.0 * math.Abs(result.val*gsl.Float64Eps)
		return nil
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Atanint(x float64) float64 {
	result := new(Result)
	status := Atanint_e(x, result)
	return EvalResult(result, status)
}
