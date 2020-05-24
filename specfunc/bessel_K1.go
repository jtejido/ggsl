/* specfunc/bessel_K1.c
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
	"github.com/jtejido/ggsl/poly"
	"math"
)

var (
	/*
	 Minimax rational approximation for [0,1), peak relative error = 1.83*GSL_DBL_EPSILON.
	 Source: http://www.advanpix.com/?p=3987
	*/
	k1_poly = []float64{
		-3.0796575782920622440538935e-01,
		-8.5370719728650778045782736e-02,
		-4.6421827664715603298154971e-03,
		-1.1253607036630425931072996e-04,
		-1.5592887702110907110292728e-06,
		-1.4030163679125934402498239e-08,
		-8.8718998640336832196558868e-11,
		-4.1614323580221539328960335e-13,
		-1.5261293392975541707230366e-15,
	}

	i1_poly = []float64{
		8.3333333333333325191635191e-02,
		6.9444444444467956461838830e-03,
		3.4722222211230452695165215e-04,
		1.1574075952009842696580084e-05,
		2.7555870002088181016676934e-07,
		4.9724386164128529514040614e-09,
	}
	/*
	   Chebyshev expansion for [1,8], peak relative error = 1.28*GSL_DBL_EPSILON.
	   Source: Pavel Holoborodko.
	*/
	ak1d = []float64{
		+2.07996868001418246e-01,
		+1.62581565017881476e-01,
		-5.87070423518863640e-03,
		+4.95021520115789501e-04,
		-5.78958347598556986e-05,
		+8.18614610209334726e-06,
		-1.31604832009487277e-06,
		+2.32546031520101213e-07,
		-4.42206518311557987e-08,
		+8.92163994883100361e-09,
		-1.89046270526983427e-09,
		+4.17568808108504702e-10,
		-9.55912361791375794e-11,
		+2.25769353153867758e-11,
		-5.48128000211158482e-12,
		+1.36386122546441926e-12,
		-3.46936690565986409e-13,
		+9.00354564415705942e-14,
		-2.37950577776254432e-14,
		+6.39447503964025336e-15,
		-1.74498363492322044e-15,
		+4.82994547989290473e-16,
		-1.35460927805445606e-16,
		+3.84604274446777234e-17,
		-1.10456856122581316e-17,
	}

	ak1 = &chebyshevSeries{
		ak1d,
		24,
		-1, 1,
		9,
	}
	/*
	   Chebyshev expansion for [8,inf), peak relative error = 1.25*GSL_DBL_EPSILON.
	   Source: SLATEC/dbsk1e.f
	*/
	ak12d = []float64{
		+.637930834373900104e-1,
		+.283288781304972094e-1,
		-.247537067390525035e-3,
		+.577197245160724882e-5,
		-.206893921953654830e-6,
		+.973998344138180418e-8,
		-.558533614038062498e-9,
		+.373299663404618524e-10,
		-.282505196102322545e-11,
		+.237201900248414417e-12,
		-.217667738799175398e-13,
		+.215791416161603245e-14,
		-.229019693071826928e-15,
		+.258288572982327496e-16,
	}

	ak12 = &chebyshevSeries{
		ak12d,
		13,
		-1, 1,
		7,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_K1_scaled_e(x float64, result *Result) err.GSLError {

	if x <= 0.0 {
		return DomainError(result)
	} else if x < 2.0*gsl.MinFloat64 {
		return OverflowError(result)
	} else if x < 1.0 {
		lx := math.Log(x)
		ex := math.Exp(x)
		x2 := x * x
		t := 0.25 * x2
		i1 := 0.5 * x * (1.0 + t*(0.5+t*poly.PolyEval(i1_poly, 6, t)))

		result.val = ex * (x2*poly.PolyEval(k1_poly, 9, x2) + x*lx*i1 + 1) / x
		result.err = ex * (1.6 + math.Abs(lx)*0.6) * gsl.Float64Eps
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x <= 8.0 {
		sx := math.Sqrt(x)
		c := new(Result)
		ak1.Evaluate((16.0/x-9.0)/7.0, c)

		result.val = (1.375 + c.val) / sx /* 1.375 = 11/8 */
		result.err = c.err / sx
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}

	sx := math.Sqrt(x)
	c := new(Result)
	ak12.Evaluate(16.0/x-1.0, c)

	result.val = (1.25 + c.val) / sx
	result.err = c.err / sx
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return nil

}

func Bessel_K1_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if x < 2.0*gsl.MinFloat64 {
		return OverflowError(result)
	} else if x < 1.0 {
		lx := math.Log(x)
		x2 := x * x
		t := 0.25 * x2
		i1 := 0.5 * x * (1.0 + t*(0.5+t*poly.PolyEval(i1_poly, 6, t)))

		result.val = (x2*poly.PolyEval(k1_poly, 9, x2) + x*lx*i1 + 1) / x
		result.err = (1.6 + math.Abs(lx)*0.6) * gsl.Float64Eps
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}

	K1_scaled := new(Result)
	stat_K1 := Bessel_K1_scaled_e(x, K1_scaled)
	stat_e := Exp_mult_err_e(-x, 0.0, K1_scaled.val, K1_scaled.err, result)
	result.err = math.Abs(result.val) * (gsl.Float64Eps*math.Abs(x) + K1_scaled.err/K1_scaled.val)
	return err.ErrorSelect(stat_e, stat_K1)

}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_K1_scaled(x float64) float64 {
	result := new(Result)
	status := Bessel_K1_scaled_e(x, result)
	return EvalResult(result, status)
}

func Bessel_K1(x float64) float64 {
	result := new(Result)
	status := Bessel_K1_e(x, result)
	return EvalResult(result, status)
}
