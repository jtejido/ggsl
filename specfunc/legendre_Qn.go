/* specfunc/legendre_Qn.c
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

/* Evaluate f_{ell+1}/f_ell
 * f_ell := Q^{b}_{a+ell}(x)
 * x > 1
 */
func legendreQ_CF1_xgt1(ell int, a, b, x float64, result *float64) err.GSLError {
	var (
		RECUR_BIG = gsl.SqrtMaxFloat64
		maxiter   = 5000
		n         = 1
		Anm2      = 1.0
		Bnm2      = 0.0
		Anm1      = 0.0
		Bnm1      = 1.0
		fell      = float64(ell)
		a1        = fell + 1.0 + a + b
		b1        = (2.0*(fell+1.0+a) + 1.0) * x
		An        = b1*Anm1 + a1*Anm2
		Bn        = b1*Bnm1 + a1*Bnm2
		an, bn    float64
		fn        = An / Bn
	)

	for n < maxiter {
		var old_fn, del, lna float64
		n++
		Anm2 = Anm1
		Bnm2 = Bnm1
		Anm1 = An
		Bnm1 = Bn
		lna = fell + float64(n) + a
		an = b*b - lna*lna
		bn = (2.0*lna + 1.0) * x
		An = bn*Anm1 + an*Anm2
		Bn = bn*Bnm1 + an*Bnm2

		if math.Abs(An) > RECUR_BIG || math.Abs(Bn) > RECUR_BIG {
			An /= RECUR_BIG
			Bn /= RECUR_BIG
			Anm1 /= RECUR_BIG
			Bnm1 /= RECUR_BIG
			Anm2 /= RECUR_BIG
			Bnm2 /= RECUR_BIG
		}

		old_fn = fn
		fn = An / Bn
		del = old_fn / fn

		if math.Abs(del-1.0) < 4.0*gsl.Float64Eps {
			break
		}
	}

	*result = fn

	if n == maxiter {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

/* Uniform asymptotic for Q_l(x).
 * Assumes x > -1.0 and x != 1.0.
 * Discards second order and higher terms.
 */
func legendre_Ql_asymp_unif(ell, x float64, result *Result) err.GSLError {
	if x < 1.0 {
		u := ell + 0.5
		th := math.Acos(x)
		var Y0, Y1 Result
		var stat_Y0, stat_Y1, stat_m err.GSLError
		var pre, B00, sum float64

		/* B00 = 1/8 (1 - th cot(th) / th^2
		 * pre = math.Sqrt(th/sin(th))
		 */
		if th < gsl.Root4Float64Eps {
			B00 = (1.0 + th*th/15.0) / 24.0
			pre = 1.0 + th*th/12.0
		} else {
			sin_th := math.Sqrt(1.0 - x*x)
			cot_th := x / sin_th
			B00 = 1.0 / 8.0 * (1.0 - th*cot_th) / (th * th)
			pre = math.Sqrt(th / sin_th)
		}

		stat_Y0 = Bessel_Y0_e(u*th, &Y0)
		stat_Y1 = Bessel_Y1_e(u*th, &Y1)

		sum = -0.5 * gsl.Pi * (Y0.val + th/u*Y1.val*B00)

		stat_m = Multiply_e(pre, sum, result)
		result.err += 0.5 * gsl.Pi * math.Abs(pre) * (Y0.err + math.Abs(th/u*B00)*Y1.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)

		return err.ErrorSelect(stat_m, stat_Y0, stat_Y1)
	} else {
		u := ell + 0.5
		xi := math.Acosh(x)
		var K0_scaled, K1_scaled Result
		var stat_K0, stat_K1, stat_e err.GSLError
		var pre, B00, sum float64

		/* B00 = -1/8 (1 - xi coth(xi) / xi^2
		 * pre = math.Sqrt(xi/sinh(xi))
		 */
		if xi < gsl.Root4Float64Eps {
			B00 = (1.0 - xi*xi/15.0) / 24.0
			pre = 1.0 - xi*xi/12.0
		} else {
			sinh_xi := math.Sqrt(x*x - 1.0)
			coth_xi := x / sinh_xi
			B00 = -1.0 / 8.0 * (1.0 - xi*coth_xi) / (xi * xi)
			pre = math.Sqrt(xi / sinh_xi)
		}

		stat_K0 = Bessel_K0_scaled_e(u*xi, &K0_scaled)
		stat_K1 = Bessel_K1_scaled_e(u*xi, &K1_scaled)

		sum = K0_scaled.val - xi/u*K1_scaled.val*B00

		stat_e = Exp_mult_e(-u*xi, pre*sum, result)
		result.err = gsl.Float64Eps * math.Abs(result.val) * math.Abs(u*xi)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

		return err.ErrorSelect(stat_e, stat_K0, stat_K1)
	}
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Legendre_Q0_e(x float64, result *Result) err.GSLError {
	if x <= -1.0 || x == 1.0 {
		return DomainError(result)
	} else if x*x < gsl.Root6Float64Eps { /* |x| <~ 0.05 */
		var (
			c3     = 1.0 / 3.0
			c5     = 1.0 / 5.0
			c7     = 1.0 / 7.0
			c9     = 1.0 / 9.0
			c11    = 1.0 / 11.0
			y      = x * x
			series = 1.0 + y*(c3+y*(c5+y*(c7+y*(c9+y*c11))))
		)
		result.val = x * series
		result.err = 2.0 * gsl.Float64Eps * math.Abs(x)
		return nil
	} else if x < 1.0 {
		result.val = 0.5 * math.Log((1.0+x)/(1.0-x))
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < 10.0 {
		result.val = 0.5 * math.Log((x+1.0)/(x-1.0))
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x*gsl.MinFloat64 < 2.0 {
		var (
			y  = 1.0 / (x * x)
			c1 = 1.0 / 3.0
			c2 = 1.0 / 5.0
			c3 = 1.0 / 7.0
			c4 = 1.0 / 9.0
			c5 = 1.0 / 11.0
			c6 = 1.0 / 13.0
			c7 = 1.0 / 15.0
		)
		result.val = (1.0 / x) * (1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7)))))))
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		return UnderflowError(result)
	}
}

func Legendre_Q1_e(x float64, result *Result) err.GSLError {

	if x <= -1.0 || x == 1.0 {
		return DomainError(result)
	} else if x*x < gsl.Root6Float64Eps { /* |x| <~ 0.05 */
		var (
			c3     = 1.0 / 3.0
			c5     = 1.0 / 5.0
			c7     = 1.0 / 7.0
			c9     = 1.0 / 9.0
			c11    = 1.0 / 11.0
			y      = x * x
			series = 1.0 + y*(c3+y*(c5+y*(c7+y*(c9+y*c11))))
		)
		result.val = x*x*series - 1.0
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < 1.0 {
		result.val = 0.5*x*(math.Log((1.0+x)/(1.0-x))) - 1.0
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < 6.0 {
		result.val = 0.5*x*math.Log((x+1.0)/(x-1.0)) - 1.0
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x*gsl.SqrtMinFloat64 < 0.99/gsl.Sqrt3 {
		var (
			y   = 1 / (x * x)
			c1  = 3.0 / 5.0
			c2  = 3.0 / 7.0
			c3  = 3.0 / 9.0
			c4  = 3.0 / 11.0
			c5  = 3.0 / 13.0
			c6  = 3.0 / 15.0
			c7  = 3.0 / 17.0
			c8  = 3.0 / 19.0
			sum = 1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*c8)))))))
		)
		result.val = sum / (3.0 * x * x)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		return UnderflowError(result)
	}
}

func Legendre_Ql_e(l int, x float64, result *Result) err.GSLError {

	if x <= -1.0 || x == 1.0 || l < 0 {
		return DomainError(result)
	} else if l == 0 {
		return Legendre_Q0_e(x, result)
	} else if l == 1 {
		return Legendre_Q1_e(x, result)
	} else if l > 100000 {
		return legendre_Ql_asymp_unif(float64(l), x, result)
	} else if x < 1.0 {
		/* Forward recurrence.
		 */
		var Q0, Q1 Result
		stat_Q0 := Legendre_Q0_e(x, &Q0)
		stat_Q1 := Legendre_Q1_e(x, &Q1)
		Qellm1 := Q0.val
		Qell := Q1.val
		var Qellp1 float64
		var ell int
		for ell = 1; ell < l; ell++ {
			fell := float64(ell)
			Qellp1 = (x*(2.0*fell+1.0)*Qell - fell*Qellm1) / (fell + 1.0)
			Qellm1 = Qell
			Qell = Qellp1
		}
		result.val = Qell
		result.err = gsl.Float64Eps * float64(l) * math.Abs(result.val)
		return err.ErrorSelect(stat_Q0, stat_Q1)
	} else {
		/* x > 1.0 */
		var rat float64
		stat_CF1 := legendreQ_CF1_xgt1(l, 0.0, 0.0, x, &rat)
		var stat_Q err.GSLError
		Qellp1 := rat * gsl.SqrtMinFloat64
		Qell := gsl.SqrtMinFloat64
		var Qellm1 float64
		var ell int
		for ell = l; ell > 0; ell-- {
			fell := float64(ell)
			Qellm1 = (x*(2.0*fell+1.0)*Qell - (fell+1.0)*Qellp1) / fell
			Qellp1 = Qell
			Qell = Qellm1
		}

		if math.Abs(Qell) > math.Abs(Qellp1) {
			var Q0 Result
			stat_Q = Legendre_Q0_e(x, &Q0)
			result.val = gsl.SqrtMinFloat64 * Q0.val / Qell
			result.err = float64(l) * gsl.Float64Eps * math.Abs(result.val)
		} else {
			var Q1 Result
			stat_Q = Legendre_Q1_e(x, &Q1)
			result.val = gsl.SqrtMinFloat64 * Q1.val / Qellp1
			result.err = float64(l) * gsl.Float64Eps * math.Abs(result.val)
		}

		return err.ErrorSelect(stat_Q, stat_CF1)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Legendre_Q0(x float64) float64 {
	result := new(Result)
	status := Legendre_Q0_e(x, result)
	return EvalResult(result, status)
}

func Legendre_Q1(x float64) float64 {
	result := new(Result)
	status := Legendre_Q1_e(x, result)
	return EvalResult(result, status)
}

func Legendre_Ql(l int, x float64) float64 {
	result := new(Result)
	status := Legendre_Ql_e(l, x, result)
	return EvalResult(result, status)
}
