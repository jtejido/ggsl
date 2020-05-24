/* specfunc/debye.c
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

/* augmented to n=5 and 6 2005-11-08 by R. J. Mathar, http://www.strw.leidenuniv.nl/~mathar */
package specfunc

import (
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"math"
)

var (
	adeb1_data = []float64{
		2.4006597190381410194,
		0.1937213042189360089,
		-0.62329124554895770e-02,
		0.3511174770206480e-03,
		-0.228222466701231e-04,
		0.15805467875030e-05,
		-0.1135378197072e-06,
		0.83583361188e-08,
		-0.6264424787e-09,
		0.476033489e-10,
		-0.36574154e-11,
		0.2835431e-12,
		-0.221473e-13,
		0.17409e-14,
		-0.1376e-15,
		0.109e-16,
		-0.9e-18,
	}

	adeb1_cs = &chebyshevSeries{
		adeb1_data,
		16,
		-1.0, 1.0,
		9,
	}

	adeb2_data = []float64{
		2.5943810232570770282,
		0.2863357204530719834,
		-0.102062656158046713e-01,
		0.6049109775346844e-03,
		-0.405257658950210e-04,
		0.28633826328811e-05,
		-0.2086394303065e-06,
		0.155237875826e-07,
		-0.11731280087e-08,
		0.897358589e-10,
		-0.69317614e-11,
		0.5398057e-12,
		-0.423241e-13,
		0.33378e-14,
		-0.2645e-15,
		0.211e-16,
		-0.17e-17,
		0.1e-18,
	}

	adeb2_cs = &chebyshevSeries{
		adeb2_data,
		17,
		-1.0, 1.0,
		10,
	}

	adeb3_data = []float64{
		2.707737068327440945,
		0.340068135211091751,
		-0.12945150184440869e-01,
		0.7963755380173816e-03,
		-0.546360009590824e-04,
		0.39243019598805e-05,
		-0.2894032823539e-06,
		0.217317613962e-07,
		-0.16542099950e-08,
		0.1272796189e-09,
		-0.987963460e-11,
		0.7725074e-12,
		-0.607797e-13,
		0.48076e-14,
		-0.3820e-15,
		0.305e-16,
		-0.24e-17,
	}

	adeb3_cs = &chebyshevSeries{
		adeb3_data,
		16,
		-1.0, 1.0,
		10,
	}

	adeb4_data = []float64{
		2.781869415020523460,
		0.374976783526892863,
		-0.14940907399031583e-01,
		0.945679811437042e-03,
		-0.66132916138933e-04,
		0.4815632982144e-05,
		-0.3588083958759e-06,
		0.271601187416e-07,
		-0.20807099122e-08,
		0.1609383869e-09,
		-0.125470979e-10,
		0.9847265e-12,
		-0.777237e-13,
		0.61648e-14,
		-0.4911e-15,
		0.393e-16,
		-0.32e-17,
	}

	adeb4_cs = &chebyshevSeries{
		adeb4_data,
		16,
		-1.0, 1.0,
		10,
	}

	adeb5_data = []float64{
		2.8340269546834530149,
		0.3994098857106266445,
		-0.164566764773099646e-1,
		0.10652138340664541e-2,
		-0.756730374875418e-4,
		0.55745985240273e-5,
		-0.4190692330918e-6,
		0.319456143678e-7,
		-0.24613318171e-8,
		0.1912801633e-9,
		-0.149720049e-10,
		0.11790312e-11,
		-0.933329e-13,
		0.74218e-14,
		-0.5925e-15,
		0.475e-16,
		-0.39e-17,
	}

	adeb5_cs = &chebyshevSeries{
		adeb5_data,
		16,
		-1.0, 1.0,
		10,
	}

	adeb6_data = []float64{
		2.8726727134130122113,
		0.4174375352339027746,
		-0.176453849354067873e-1,
		0.11629852733494556e-2,
		-0.837118027357117e-4,
		0.62283611596189e-5,
		-0.4718644465636e-6,
		0.361950397806e-7,
		-0.28030368010e-8,
		0.2187681983e-9,
		-0.171857387e-10,
		0.13575809e-11,
		-0.1077580e-12,
		0.85893e-14,
		-0.6872e-15,
		0.552e-16,
		-0.44e-17,
	}

	adeb6_cs = &chebyshevSeries{
		adeb6_data,
		16,
		-1.0, 1.0,
		10,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Debye_1_e(x float64, result *Result) err.GSLError {
	val_infinity := 1.64493406684822644
	xcut := -gsl.LnMinFloat64

	if x < 0.0 {
		return DomainError(result)
	} else if x < 2.0*gsl.SqrtFloat64Eps {
		result.val = 1.0 - 0.25*x + x*x/36.0
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x <= 4.0 {
		t := x*x/8.0 - 1.0
		var c Result
		adeb1_cs.Evaluate(t, &c)
		result.val = c.val - 0.25*x
		result.err = c.err + 0.25*x*gsl.Float64Eps
		return nil
	} else if x < -(gsl.Ln2 + gsl.LnFloat64Eps) {
		nexp := math.Floor(xcut / x)
		ex := math.Exp(-x)
		sum := 0.0
		xk := nexp * x
		rk := nexp
		for i := nexp; i >= 1; i-- {
			sum *= ex
			sum += (1.0 + 1.0/xk) / rk
			rk -= 1.0
			xk -= x
		}
		result.val = val_infinity/x - sum*ex
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < xcut {
		result.val = (val_infinity - math.Exp(-x)*(x+1.0)) / x
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = val_infinity / x
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Debye_2_e(x float64, result *Result) err.GSLError {
	val_infinity := 4.80822761263837714
	xcut := -gsl.LnMinFloat64

	if x < 0.0 {
		return DomainError(result)
	} else if x < 2.0*gsl.Sqrt2*gsl.SqrtFloat64Eps {
		result.val = 1.0 - x/3.0 + x*x/24.0
		result.err = gsl.Float64Eps * result.val
		return nil
	} else if x <= 4.0 {
		t := x*x/8.0 - 1.0
		var c Result
		adeb2_cs.Evaluate(t, &c)
		result.val = c.val - x/3.0
		result.err = c.err + gsl.Float64Eps*x/3.0
		return nil
	} else if x < -(gsl.Ln2 + gsl.LnFloat64Eps) {
		nexp := math.Floor(xcut / x)
		ex := math.Exp(-x)
		xk := nexp * x
		rk := nexp
		sum := 0.0

		for i := nexp; i >= 1; i-- {
			sum *= ex
			sum += (1.0 + 2.0/xk + 2.0/(xk*xk)) / rk
			rk -= 1.0
			xk -= x
		}
		result.val = val_infinity/(x*x) - 2.0*sum*ex
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < xcut {
		x2 := x * x
		sum := 2.0 + 2.0*x + x2
		result.val = (val_infinity - 2.0*sum*math.Exp(-x)) / x2
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = (val_infinity / x) / x
		result.err = gsl.Float64Eps * result.val
		CheckUnderflow(result)
		return nil
	}
}

func Debye_3_e(x float64, result *Result) err.GSLError {
	val_infinity := 19.4818182068004875
	xcut := -gsl.LnMinFloat64

	if x < 0.0 {
		return DomainError(result)
	} else if x < 2.0*gsl.Sqrt2*gsl.SqrtFloat64Eps {
		result.val = 1.0 - 3.0*x/8.0 + x*x/20.0
		result.err = gsl.Float64Eps * result.val
		return nil
	} else if x <= 4.0 {
		t := x*x/8.0 - 1.0
		var c Result
		adeb3_cs.Evaluate(t, &c)
		result.val = c.val - 0.375*x
		result.err = c.err + gsl.Float64Eps*0.375*x
		return nil
	} else if x < -(gsl.Ln2 + gsl.LnFloat64Eps) {
		nexp := math.Floor(xcut / x)
		ex := math.Exp(-x)
		xk := nexp * x
		rk := nexp
		sum := 0.0

		for i := nexp; i >= 1; i-- {
			xk_inv := 1.0 / xk
			sum *= ex
			sum += (((6.0*xk_inv+6.0)*xk_inv+3.0)*xk_inv + 1.0) / rk
			rk -= 1.0
			xk -= x
		}
		result.val = val_infinity/(x*x*x) - 3.0*sum*ex
		result.err = gsl.Float64Eps * result.val
		return nil
	} else if x < xcut {
		x3 := x * x * x
		sum := 6.0 + 6.0*x + 3.0*x*x + x3
		result.val = (val_infinity - 3.0*sum*math.Exp(-x)) / x3
		result.err = gsl.Float64Eps * result.val
		return nil
	} else {
		result.val = ((val_infinity / x) / x) / x
		result.err = gsl.Float64Eps * result.val
		CheckUnderflow(result)
		return nil
	}
}

func Debye_4_e(x float64, result *Result) err.GSLError {
	val_infinity := 99.5450644937635129
	xcut := -gsl.LnMinFloat64

	if x < 0.0 {
		return DomainError(result)
	} else if x < 2.0*gsl.Sqrt2*gsl.SqrtFloat64Eps {
		result.val = 1.0 - 2.0*x/5.0 + x*x/18.0
		result.err = gsl.Float64Eps * result.val
		return nil
	} else if x <= 4.0 {
		t := x*x/8.0 - 1.0
		var c Result
		adeb4_cs.Evaluate(t, &c)
		result.val = c.val - 2.0*x/5.0
		result.err = c.err + gsl.Float64Eps*2.0*x/5.0
		return nil
	} else if x < -(gsl.Ln2 + gsl.LnFloat64Eps) {
		nexp := math.Floor(xcut / x)
		ex := math.Exp(-x)
		xk := nexp * x
		rk := nexp
		sum := 0.0

		for i := nexp; i >= 1; i-- {
			xk_inv := 1.0 / xk
			sum *= ex
			sum += ((((24.0*xk_inv+24.0)*xk_inv+12.0)*xk_inv+4.0)*xk_inv + 1.0) / rk
			rk -= 1.0
			xk -= x
		}
		result.val = val_infinity/(x*x*x*x) - 4.0*sum*ex
		result.err = gsl.Float64Eps * result.val
		return nil
	} else if x < xcut {
		x2 := x * x
		x4 := x2 * x2
		sum := 24.0 + 24.0*x + 12.0*x2 + 4.0*x2*x + x4
		result.val = (val_infinity - 4.0*sum*math.Exp(-x)) / x4
		result.err = gsl.Float64Eps * result.val
		return nil
	} else {
		result.val = (((val_infinity / x) / x) / x) / x
		result.err = gsl.Float64Eps * result.val
		CheckUnderflow(result)
		return nil
	}
}

func Debye_5_e(x float64, result *Result) err.GSLError {
	val_infinity := 610.405837190669483828710757875
	xcut := -gsl.LnMinFloat64

	if x < 0.0 {
		return DomainError(result)
	} else if x < 2.0*gsl.Sqrt2*gsl.SqrtFloat64Eps {
		result.val = 1.0 - 5.0*x/12.0 + 5.0*x*x/84.0
		result.err = gsl.Float64Eps * result.val
		return nil
	} else if x <= 4.0 {
		t := x*x/8.0 - 1.0
		var c Result
		adeb5_cs.Evaluate(t, &c)
		result.val = c.val - 5.0*x/12.0
		result.err = c.err + gsl.Float64Eps*5.0*x/12.0
		return nil
	} else if x < -(gsl.Ln2 + gsl.LnFloat64Eps) {
		nexp := math.Floor(xcut / x)
		ex := math.Exp(-x)
		xk := nexp * x
		rk := nexp
		sum := 0.0

		for i := nexp; i >= 1; i-- {
			xk_inv := 1.0 / xk
			sum *= ex
			sum += (((((120.0*xk_inv+120.0)*xk_inv+60.0)*xk_inv+20.0)*xk_inv+5.0)*xk_inv + 1.0) / rk
			rk -= 1.0
			xk -= x
		}
		result.val = val_infinity/(x*x*x*x*x) - 5.0*sum*ex
		result.err = gsl.Float64Eps * result.val
		return nil
	} else if x < xcut {
		x2 := x * x
		x4 := x2 * x2
		x5 := x4 * x
		sum := 120.0 + 120.0*x + 60.0*x2 + 20.0*x2*x + 5.0*x4 + x5
		result.val = (val_infinity - 5.0*sum*math.Exp(-x)) / x5
		result.err = gsl.Float64Eps * result.val
		return nil
	} else {
		result.val = ((((val_infinity / x) / x) / x) / x) / x
		result.err = gsl.Float64Eps * result.val
		CheckUnderflow(result)
		return nil
	}
}

func Debye_6_e(x float64, result *Result) err.GSLError {
	val_infinity := 4356.06887828990661194792541535
	xcut := -gsl.LnMinFloat64

	/* CHECK_POINTER(result) */

	if x < 0.0 {
		return DomainError(result)
	} else if x < 2.0*gsl.Sqrt2*gsl.SqrtFloat64Eps {
		result.val = 1.0 - 3.0*x/7.0 + x*x/16.0
		result.err = gsl.Float64Eps * result.val
		return nil
	} else if x <= 4.0 {
		t := x*x/8.0 - 1.0
		var c Result
		adeb6_cs.Evaluate(t, &c)
		result.val = c.val - 3.0*x/7.0
		result.err = c.err + gsl.Float64Eps*3.0*x/7.0
		return nil
	} else if x < -(gsl.Ln2 + gsl.LnFloat64Eps) {
		nexp := math.Floor(xcut / x)
		ex := math.Exp(-x)
		xk := nexp * x
		rk := nexp
		sum := 0.0

		for i := nexp; i >= 1; i-- {
			xk_inv := 1.0 / xk
			sum *= ex
			sum += ((((((720.0*xk_inv+720.0)*xk_inv+360.0)*xk_inv+120.0)*xk_inv+30.0)*xk_inv+6.0)*xk_inv + 1.0) / rk
			rk -= 1.0
			xk -= x
		}
		result.val = val_infinity/(x*x*x*x*x*x) - 6.0*sum*ex
		result.err = gsl.Float64Eps * result.val
		return nil
	} else if x < xcut {
		x2 := x * x
		x4 := x2 * x2
		x6 := x4 * x2
		sum := 720.0 + 720.0*x + 360.0*x2 + 120.0*x2*x + 30.0*x4 + 6.0*x4*x + x6
		result.val = (val_infinity - 6.0*sum*math.Exp(-x)) / x6
		result.err = gsl.Float64Eps * result.val
		return nil
	} else {
		result.val = (((((val_infinity / x) / x) / x) / x) / x) / x
		result.err = gsl.Float64Eps * result.val
		CheckUnderflow(result)
		return nil
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

func Debye_1(x float64) float64 {
	result := new(Result)
	status := Debye_1_e(x, result)
	return EvalResult(result, status)
}

func Debye_2(x float64) float64 {
	result := new(Result)
	status := Debye_2_e(x, result)
	return EvalResult(result, status)
}

func Debye_3(x float64) float64 {
	result := new(Result)
	status := Debye_3_e(x, result)
	return EvalResult(result, status)
}

func Debye_4(x float64) float64 {
	result := new(Result)
	status := Debye_4_e(x, result)
	return EvalResult(result, status)
}

func Debye_5(x float64) float64 {
	result := new(Result)
	status := Debye_5_e(x, result)
	return EvalResult(result, status)
}

func Debye_6(x float64) float64 {
	result := new(Result)
	status := Debye_6_e(x, result)
	return EvalResult(result, status)
}
