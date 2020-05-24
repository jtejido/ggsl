/* specfunc/shint.c
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

/* based on SLATEC shi.f, W. Fullerton

series for shi  on the interval  0.00000e+00 to  1.40625e-01
                                       with weighted error   4.67e-20
                                        log weighted error  19.33
                              significant figures required  17.07
                                   decimal places required  19.75
*/
var (
	shi_data = []float64{
		0.0078372685688900950695,
		0.0039227664934234563973,
		0.0000041346787887617267,
		0.0000000024707480372883,
		0.0000000000009379295591,
		0.0000000000000002451817,
		0.0000000000000000000467,
	}

	shi_cs = &chebyshevSeries{
		shi_data,
		6,
		-1, 1,
		6,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Shi_e(x float64, result *Result) err.GSLError {
	xsml := gsl.SqrtFloat64Eps /* sqrt (r1mach(3)) */
	ax := math.Abs(x)

	if ax < xsml {
		result.val = x
		result.err = 0.0
		return nil
	} else if ax <= 0.375 {
		result_c := new(Result)
		shi_cs.Evaluate(128.0*x*x/9.0-1.0, result_c)
		result.val = x * (1.0 + result_c.val)
		result.err = x * result_c.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result_Ei, result_E1 := new(Result), new(Result)
		status_Ei := Expint_Ei_e(x, result_Ei)
		status_E1 := Expint_E1_e(x, result_E1)
		var stat_E1, stat_Ei int
		if status_Ei != nil {
			stat_Ei = status_Ei.Status()
		}

		if status_E1 != nil {
			stat_E1 = status_E1.Status()
		}

		result.val = 0.5 * (result_Ei.val + result_E1.val)
		result.err = 0.5 * (result_Ei.err + result_E1.err)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

		if stat_Ei == err.EUNDRFLW && stat_E1 == err.EUNDRFLW {
			return err.ERROR("underflow", err.EUNDRFLW)
		} else if stat_Ei == err.EOVRFLW || stat_E1 == err.EOVRFLW {
			return err.ERROR("overflow", err.EOVRFLW)
		} else {
			return nil
		}
	}
}

func Chi_e(x float64, result *Result) err.GSLError {
	result_Ei, result_E1 := new(Result), new(Result)
	status_Ei := Expint_Ei_e(x, result_Ei)
	status_E1 := Expint_E1_e(x, result_E1)
	var stat_E1, stat_Ei int
	if status_Ei != nil {
		stat_Ei = status_Ei.Status()
	}

	if status_E1 != nil {
		stat_E1 = status_E1.Status()
	}

	if stat_Ei == err.EDOM || stat_E1 == err.EDOM {
		return DomainError(result)
	} else if stat_Ei == err.EUNDRFLW && stat_E1 == err.EUNDRFLW {
		return UnderflowError(result)
	} else if stat_Ei == err.EOVRFLW || stat_E1 == err.EOVRFLW {
		return OverflowError(result)
	} else {
		result.val = 0.5 * (result_Ei.val - result_E1.val)
		result.err = 0.5 * (result_Ei.err + result_E1.err)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Shi(x float64) float64 {
	result := new(Result)
	status := Shi_e(x, result)
	return EvalResult(result, status)
}

func Chi(x float64) float64 {
	result := new(Result)
	status := Chi_e(x, result)
	return EvalResult(result, status)
}
