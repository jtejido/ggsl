/* specfunc/bessel_Y0.c
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
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"math"
)

/* based on SLATEC besy0, 1980 version, w. fullerton */

/* chebyshev expansions

series for by0        on the interval  0.          to  1.60000d+01
                                       with weighted error   1.20e-17
                                        log weighted error  16.92
                              significant figures required  16.15
                                   decimal places required  17.48
*/

var (
	by0_data = []float64{
		-0.011277839392865573,
		-0.128345237560420350,
		-0.104378847997942490,
		0.023662749183969695,
		-0.002090391647700486,
		0.000103975453939057,
		-0.000003369747162423,
		0.000000077293842676,
		-0.000000001324976772,
		0.000000000017648232,
		-0.000000000000188105,
		0.000000000000001641,
		-0.000000000000000011,
	}

	by0_cs = &chebyshevSeries{
		by0_data,
		12,
		-1, 1,
		8,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_Y0_e(x float64, result *Result) err.GSLError {
	two_over_pi := 2.0 / gsl.Pi
	xmax := 1.0 / gsl.Float64Eps

	/* CHECK_POINTER(result) */

	if x <= 0.0 {
		return DomainError(result)
	} else if x < 4.0 {
		J0, c := new(Result), new(Result)
		stat_J0 := Bessel_J0_e(x, J0)
		by0_cs.Evaluate(0.125*x*x-1.0, c)
		result.val = two_over_pi*(-math.Ln2+math.Log(x))*J0.val + 0.375 + c.val
		result.err = 2.0*gsl.Float64Eps*math.Abs(result.val) + c.err
		return stat_J0
	} else if x < xmax {
		/* Leading behaviour of phase is x, which is exact,
		 * so the error is bounded.
		 */
		z := 32.0/(x*x) - 1.0
		c1, c2, sp := new(Result), new(Result), new(Result)
		stat_c1 := bm0.Evaluate(z, c1)
		stat_c2 := bth0.Evaluate(z, c2)
		stat_sp := bessel_sin_pi4_e(x, c2.val/x, sp)
		sqrtx := math.Sqrt(x)
		ampl := (0.75 + c1.val) / sqrtx
		result.val = ampl * sp.val
		result.err = math.Abs(sp.val)*c1.err/sqrtx + math.Abs(ampl)*sp.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_sp, stat_c1, stat_c2)
	} else {
		return UnderflowError(result)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_Y0(x float64) float64 {
	result := new(Result)
	status := Bessel_Y0_e(x, result)
	return EvalResult(result, status)

}
