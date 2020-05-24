/* specfunc/bessel_J0.c
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

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC besj0, 1977 version, w. fullerton */

/* chebyshev expansions for Bessel functions

series for bj0        on the interval  0.          to  1.60000d+01
                                       with weighted error   7.47e-18
                                        log weighted error  17.13
                              significant figures required  16.98
                                   decimal places required  17.68

*/

var (
	bj0_data = []float64{
		0.100254161968939137,
		-0.665223007764405132,
		0.248983703498281314,
		-0.0332527231700357697,
		0.0023114179304694015,
		-0.0000991127741995080,
		0.0000028916708643998,
		-0.0000000612108586630,
		0.0000000009838650793,
		-0.0000000000124235515,
		0.0000000000001265433,
		-0.0000000000000010619,
		0.0000000000000000074,
	}

	bj0_cs = &chebyshevSeries{
		bj0_data,
		12,
		-1, 1,
		9,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_J0_e(x float64, result *Result) err.GSLError {
	y := math.Abs(x)

	if y < 2.0*gsl.SqrtFloat64Eps {

		result.val = 1.0
		result.err = y * y
		return nil
	} else if y <= 4.0 {
		return bj0_cs.Evaluate(0.125*y*y-1.0, result)
	} else {
		z := 32.0/(y*y) - 1.0
		ca, ct, cp := new(Result), new(Result), new(Result)
		stat_ca := bm0.Evaluate(z, ca)
		stat_ct := bth0.Evaluate(z, ct)
		stat_cp := bessel_cos_pi4_e(y, ct.val/y, cp)
		sqrty := math.Sqrt(y)
		ampl := (0.75 + ca.val) / sqrty
		result.val = ampl * cp.val
		result.err = math.Abs(cp.val)*ca.err/sqrty + math.Abs(ampl)*cp.err
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_ca, stat_ct, stat_cp)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_J0(x float64) float64 {
	result := new(Result)
	status := Bessel_J0_e(x, result)
	return EvalResult(result, status)
}
