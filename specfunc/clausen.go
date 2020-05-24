/* specfunc/clausen.c
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
	aclaus_data = []float64{
		2.142694363766688447e+00,
		0.723324281221257925e-01,
		0.101642475021151164e-02,
		0.3245250328531645e-04,
		0.133315187571472e-05,
		0.6213240591653e-07,
		0.313004135337e-08,
		0.16635723056e-09,
		0.919659293e-11,
		0.52400462e-12,
		0.3058040e-13,
		0.18197e-14,
		0.1100e-15,
		0.68e-17,
		0.4e-18,
	}
	aclaus_cs = &chebyshevSeries{
		aclaus_data,
		14,
		-1, 1,
		8, /* FIXME:  this is a guess, correct value needed here BJG */
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Clausen_e(x float64, result *Result) err.GSLError {
	x_cut := gsl.Pi * gsl.SqrtFloat64Eps

	sgn := 1.0
	var status_red err.GSLError

	if x < 0.0 {
		x = -x
		sgn = -1.0
	}

	/* Argument reduction to [0, 2pi) */
	status_red = Angle_restrict_pos_e(&x)

	/* Further reduction to [0,pi) */
	if x > gsl.Pi {
		/* simulated extra precision: 2PI = p0 + p1 */
		p0 := 6.28125
		p1 := 0.19353071795864769253e-02
		x = (p0 - x) + p1
		sgn = -sgn
	}

	if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
	} else if x < x_cut {
		result.val = x * (1.0 - math.Log(x))
		result.err = x * gsl.Float64Eps
	} else {
		t := 2.0 * (x*x/(gsl.Pi*gsl.Pi) - 0.5)
		var result_c Result
		aclaus_cs.Evaluate(t, &result_c)
		result.val = x * (result_c.val - math.Log(x))
		result.err = x * (result_c.err + gsl.Float64Eps)
	}

	result.val *= sgn

	return status_red
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Clausen(x float64) float64 {
	result := new(Result)
	status := Clausen_e(x, result)
	return EvalResult(result, status)
}
