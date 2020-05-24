/* specfunc/bessel_J1.c
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

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besj1, 1983 version, w. fullerton */

/* chebyshev expansions

series for bj1        on the interval  0.          to  1.60000d+01
                                       with weighted error   4.48e-17
                                        log weighted error  16.35
                              significant figures required  15.77
                                   decimal places required  16.89

*/
var (
	bj1_data = []float64{
		-0.11726141513332787,
		-0.25361521830790640,
		0.050127080984469569,
		-0.004631514809625081,
		0.000247996229415914,
		-0.000008678948686278,
		0.000000214293917143,
		-0.000000003936093079,
		0.000000000055911823,
		-0.000000000000632761,
		0.000000000000005840,
		-0.000000000000000044,
	}

	bj1_cs = &chebyshevSeries{
		bj1_data,
		11,
		-1, 1,
		8,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_J1_e(x float64, result *Result) err.GSLError {
	y := math.Abs(x)

	if y == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if y < 2.0*gsl.MinFloat64 {
		return UnderflowError(result)
	} else if y < ROOT_EIGHT*gsl.SqrtFloat64Eps {
		result.val = 0.5 * x
		result.err = 0.0
		return nil
	} else if y < 4.0 {
		c := new(Result)
		bj1_cs.Evaluate(0.125*y*y-1.0, c)
		result.val = x * (0.25 + c.val)
		result.err = math.Abs(x * c.err)
		return nil
	} else {
		/* Because the leading term in the phase is y,
		 * which we assume is exactly known, the error
		 * in the cos() evaluation is bounded.
		 */
		z := 32.0/(y*y) - 1.0
		ca, ct, sp := new(Result), new(Result), new(Result)
		stat_ca := bm1.Evaluate(z, ca)
		stat_ct := bth1.Evaluate(z, ct)
		stat_sp := bessel_sin_pi4_e(y, ct.val/y, sp)
		sqrty := math.Sqrt(y)
		ampl := (0.75 + ca.val) / sqrty
		result.val = sp.val
		if x < 0.0 {
			result.val *= -ampl
		} else {
			result.val *= ampl
		}
		result.err = math.Abs(sp.val)*ca.err/sqrty + math.Abs(ampl)*sp.err
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_ca, stat_ct, stat_sp)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_J1(x float64) float64 {
	result := new(Result)
	status := Bessel_J1_e(x, result)
	return EvalResult(result, status)
}
