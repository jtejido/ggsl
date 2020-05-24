/* specfunc/bessel_Y1.c
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

/* based on SLATEC besy1, 1977 version, w. fullerton */

/* chebyshev expansions

series for by1        on the interval  0.          to  1.60000d+01
                                       with weighted error   1.87e-18
                                        log weighted error  17.73
                              significant figures required  17.83
                                   decimal places required  18.30
*/
var (
	by1_data = []float64{
		0.03208047100611908629,
		1.262707897433500450,
		0.00649996189992317500,
		-0.08936164528860504117,
		0.01325088122175709545,
		-0.00089790591196483523,
		0.00003647361487958306,
		-0.00000100137438166600,
		0.00000001994539657390,
		-0.00000000030230656018,
		0.00000000000360987815,
		-0.00000000000003487488,
		0.00000000000000027838,
		-0.00000000000000000186,
	}

	by1_cs = &chebyshevSeries{
		by1_data,
		13,
		-1, 1,
		10,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_Y1_e(x float64, result *Result) err.GSLError {
	two_over_pi := 2.0 / gsl.Pi
	xmin := 1.571 * gsl.MinFloat64 /*exp ( amax1(alog(r1mach(1)), -alog(r1mach(2)))+.01) */
	x_small := 2.0 * gsl.SqrtFloat64Eps
	xmax := 1.0 / gsl.Float64Eps

	if x <= 0.0 {
		return DomainError(result)
	} else if x < xmin {
		return OverflowError(result)
	} else if x < x_small {
		lnterm := math.Log(0.5 * x)
		J1, c := new(Result), new(Result)
		status := Bessel_J1_e(x, J1)
		by1_cs.Evaluate(-1.0, c)
		result.val = two_over_pi*lnterm*J1.val + (0.5+c.val)/x
		result.err = math.Abs(lnterm)*(math.Abs(gsl.Float64Eps*J1.val)+J1.err) + c.err/x
		return status
	} else if x < 4.0 {
		lnterm := math.Log(0.5 * x)
		J1, c := new(Result), new(Result)
		by1_cs.Evaluate(0.125*x*x-1.0, c)
		status := Bessel_J1_e(x, J1)
		result.val = two_over_pi*lnterm*J1.val + (0.5+c.val)/x
		result.err = math.Abs(lnterm)*(math.Abs(gsl.Float64Eps*J1.val)+J1.err) + c.err/x
		return status
	} else if x < xmax {
		z := 32.0/(x*x) - 1.0
		ca, ct, cp := new(Result), new(Result), new(Result)
		stat_ca := bm1.Evaluate(z, ca)
		stat_ct := bth1.Evaluate(z, ct)
		stat_cp := bessel_cos_pi4_e(x, ct.val/x, cp)
		sqrtx := math.Sqrt(x)
		ampl := (0.75 + ca.val) / sqrtx
		result.val = -ampl * cp.val
		result.err = math.Abs(cp.val)*ca.err/sqrtx + math.Abs(ampl)*cp.err
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_ca, stat_ct, stat_cp)
	} else {
		return UnderflowError(result)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_Y1(x float64) float64 {
	result := new(Result)
	status := Bessel_Y1_e(x, result)
	return EvalResult(result, status)
}
