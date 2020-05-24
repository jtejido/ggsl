/* specfunc/bessel_K0.c
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
	"github.com/lucky-se7en/ggsl/poly"
	"math"
)

/*
 Minimax rational approximation for [0,1), peak relative error = 2.04*GSL_DBL_EPSILON.
 Source: http://www.advanpix.com/?p=3812
*/
var (
	k0_poly = []float64{
		1.1593151565841244842077226e-01,
		2.7898287891460317300886539e-01,
		2.5248929932161220559969776e-02,
		8.4603509072136578707676406e-04,
		1.4914719243067801775856150e-05,
		1.6271068931224552553548933e-07,
		1.2082660336282566759313543e-09,
		6.6117104672254184399933971e-12,
	}

	i0_poly = []float64{
		1.0000000000000000044974165e+00,
		2.4999999999999822316775454e-01,
		2.7777777777892149148858521e-02,
		1.7361111083544590676709592e-03,
		6.9444476047072424198677755e-05,
		1.9288265756466775034067979e-06,
		3.9908220583262192851839992e-08,
	}

	/*
	   Chebyshev expansion for [1,8], peak relative error = 1.28*GSL_DBL_EPSILON.
	   Source: Pavel Holoborodko.
	*/
	ak0d = []float64{
		-3.28737867094650101e-02,
		-4.49369057710236880e-02,
		+2.98149992004308095e-03,
		-3.03693649396187920e-04,
		+3.91085569307646836e-05,
		-5.86872422399215952e-06,
		+9.82873709937322009e-07,
		-1.78978645055651171e-07,
		+3.48332306845240957e-08,
		-7.15909210462546599e-09,
		+1.54019930048919494e-09,
		-3.44555485579194210e-10,
		+7.97356101783753023e-11,
		-1.90090968913069735e-11,
		+4.65295609304114621e-12,
		-1.16614287433470780e-12,
		+2.98554375218596891e-13,
		-7.79276979512292169e-14,
		+2.07027467168948402e-14,
		-5.58987860393825313e-15,
		+1.53202965950646914e-15,
		-4.25737536712188186e-16,
		+1.19840238501357389e-16,
		-3.41407346762502397e-17,
	}

	ak0 = &chebyshevSeries{
		ak0d,
		23,
		-1, 1,
		10,
	}

	/*
	   Chebyshev expansion for [8,inf), peak relative error = 1.25*GSL_DBL_EPSILON.
	   Source: SLATEC/dbsk0e.f
	*/
	ak02d = []float64{
		-.1201869826307592240e-1,
		-.9174852691025695311e-2,
		+.1444550931775005821e-3,
		-.4013614175435709729e-5,
		+.1567831810852310673e-6,
		-.7770110438521737710e-8,
		+.4611182576179717883e-9,
		-.3158592997860565771e-10,
		+.2435018039365041128e-11,
		-.2074331387398347898e-12,
		+.1925787280589917085e-13,
		-.1927554805838956104e-14,
		+.2062198029197818278e-15,
		-.2341685117579242403e-16,
	}

	ak02 = &chebyshevSeries{
		ak02d,
		13,
		-1, 1,
		8,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_K0_scaled_e(x float64, result *Result) err.GSLError {

	if x <= 0.0 {
		return DomainError(result)
	} else if x < 1.0 {
		lx := math.Log(x)
		ex := math.Exp(x)
		x2 := x * x
		result.val = ex * (poly.PolyEval(k0_poly, 8, x2) - lx*(1.0+0.25*x2*poly.PolyEval(i0_poly, 7, 0.25*x2)))
		result.err = ex * (1.6 + math.Abs(lx)*0.6) * gsl.Float64Eps
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x <= 8.0 {
		sx := math.Sqrt(x)
		c := new(Result)
		ak0.Evaluate((16.0/x-9.0)/7.0, c)
		result.val = (1.203125 + c.val) / sx /* 1.203125 = 77/64 */
		result.err = c.err / sx
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}

	sx := math.Sqrt(x)
	c := new(Result)
	ak02.Evaluate(16.0/x-1.0, c)
	result.val = (1.25 + c.val) / sx
	result.err = (c.err + gsl.Float64Eps) / sx
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return nil
}

func Bessel_K0_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if x < 1.0 {
		lx := math.Log(x)
		x2 := x * x

		result.val = poly.PolyEval(k0_poly, 8, x2) - lx*(1.0+0.25*x2*poly.PolyEval(i0_poly, 7, 0.25*x2))
		result.err = (1.6 + math.Abs(lx)*0.6) * gsl.Float64Eps
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}

	K0_scaled := new(Result)
	stat_K0 := Bessel_K0_scaled_e(x, K0_scaled)
	stat_e := Exp_mult_err_e(-x, gsl.Float64Eps*math.Abs(x), K0_scaled.val, K0_scaled.err, result)
	return err.ErrorSelect(stat_e, stat_K0)
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_K0_scaled(x float64) float64 {
	result := new(Result)
	status := Bessel_K0_scaled_e(x, result)
	return EvalResult(result, status)
}

func Bessel_K0(x float64) float64 {
	result := new(Result)
	status := Bessel_K0_e(x, result)
	return EvalResult(result, status)
}
