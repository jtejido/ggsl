/* specfunc/hyperg_0F1.c
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

var (
	locEPS = 1000.0 * gsl.Float64Eps
)

/* Evaluate bessel_I(nu, x), allowing nu < 0.
 * This is fine here because we do not not allow
 * nu to be a negative integer.
 * x > 0.
 */
func hyperg_0F1_bessel_I(nu, x float64, result *Result) err.GSLError {
	if x > gsl.LnMaxFloat64 {
		return OverflowError(result)
	}

	if nu < 0.0 {
		anu := -nu
		s := 2.0 / gsl.Pi * math.Sin(anu*gsl.Pi)
		ex := math.Exp(x)
		I, K := new(Result), new(Result)
		stat_I := Bessel_Inu_scaled_e(anu, x, I)
		stat_K := Bessel_Knu_scaled_e(anu, x, K)
		result.val = ex*I.val + s*(K.val/ex)
		result.err = ex*I.err + math.Abs(s*K.err/ex)
		result.err += math.Abs(s*(K.val/ex)) * gsl.Float64Eps * anu * gsl.Pi
		return err.ErrorSelect(stat_K, stat_I)
	} else {
		ex := math.Exp(x)
		I := new(Result)
		stat_I := Bessel_Inu_scaled_e(nu, x, I)
		result.val = ex * I.val
		result.err = ex*I.err + gsl.Float64Eps*math.Abs(result.val)
		return stat_I
	}
}

/* Evaluate bessel_J(nu, x), allowing nu < 0.
 * This is fine here because we do not not allow
 * nu to be a negative integer.
 * x > 0.
 */
func hyperg_0F1_bessel_J(nu, x float64, result *Result) err.GSLError {
	if nu < 0.0 {
		anu := -nu
		s := math.Sin(anu * gsl.Pi)
		c := math.Cos(anu * gsl.Pi)
		J, Y := new(Result), new(Result)
		stat_J := Bessel_Jnu_e(anu, x, J)
		stat_Y := Bessel_Ynu_e(anu, x, Y)
		result.val = c*J.val - s*Y.val
		result.err = math.Abs(c*J.err) + math.Abs(s*Y.err)
		result.err += math.Abs(anu*gsl.Pi) * gsl.Float64Eps * math.Abs(J.val+Y.val)
		return err.ErrorSelect(stat_Y, stat_J)
	}
	return Bessel_Jnu_e(nu, x, result)

}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Hyperg_0F1_e(c, x float64, result *Result) err.GSLError {
	rintc := math.Floor(c + 0.5)
	c_neg_integer := (c < 0.0 && math.Abs(c-rintc) < locEPS)

	if c == 0.0 || c_neg_integer {
		return DomainError(result)
	} else if x < 0.0 {
		Jcm1, lg_c := new(Result), new(Result)
		var sgn float64
		stat_g := Lngamma_sgn_e(c, lg_c, &sgn)
		stat_J := hyperg_0F1_bessel_J(c-1.0, 2.0*math.Sqrt(-x), Jcm1)
		if stat_g != nil {
			result.val = 0.0
			result.err = 0.0
			return stat_g
		} else if Jcm1.val == 0.0 {
			result.val = 0.0
			result.err = 0.0
			return stat_J
		} else {
			tl := math.Log(-x) * 0.5 * (1.0 - c)
			ln_pre_val := lg_c.val + tl
			ln_pre_err := lg_c.err + 2.0*gsl.Float64Eps*math.Abs(tl)
			return Exp_mult_err_e(ln_pre_val, ln_pre_err, sgn*Jcm1.val, Jcm1.err, result)
		}
	} else if x == 0.0 {
		result.val = 1.0
		result.err = 1.0
		return nil
	} else {
		Icm1, lg_c := new(Result), new(Result)
		var sgn float64

		stat_g := Lngamma_sgn_e(c, lg_c, &sgn)
		stat_I := hyperg_0F1_bessel_I(c-1.0, 2.0*math.Sqrt(x), Icm1)

		if stat_g != nil {
			result.val = 0.0
			result.err = 0.0
			return stat_g
		} else if Icm1.val == 0.0 {
			result.val = 0.0
			result.err = 0.0
			return stat_I
		} else {
			tl := math.Log(x) * 0.5 * (1.0 - c)
			ln_pre_val := lg_c.val + tl
			ln_pre_err := lg_c.err + 2.0*gsl.Float64Eps*math.Abs(tl)
			return Exp_mult_err_e(ln_pre_val, ln_pre_err, sgn*Icm1.val, Icm1.err, result)
		}
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Hyperg_0F1(c, x float64) float64 {
	result := new(Result)
	status := Hyperg_0F1_e(c, x, result)
	return EvalResult(result, status)
}
