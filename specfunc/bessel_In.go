/* specfunc/bessel_In.c
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

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_In_scaled_e(n int, x float64, result *Result) err.GSLError {
	ax := math.Abs(x)
	n = gsl.AbsInt(n)

	if n == 0 {
		return Bessel_I0_scaled_e(x, result)
	} else if n == 1 {
		return Bessel_I1_scaled_e(x, result)
	} else if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if x*x < 10.0*(float64(n)+1.0)/math.E {
		t := new(Result)
		ex := math.Exp(-ax)
		stat_In := Bessel_IJ_taylor_e(float64(n), ax, 1, 50, gsl.Float64Eps, t)
		result.val = t.val * ex
		result.err = t.err * ex
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		if x < 0.0 && gsl.IsOdd(n) {
			result.val = -result.val
		}

		return stat_In
	} else if n < 150 && ax < 1e7 {
		I0_scaled := new(Result)
		stat_I0 := Bessel_I0_scaled_e(ax, I0_scaled)
		var rat float64
		stat_CF1 := bessel_I_CF1_ser(float64(n), ax, &rat)
		Ikp1 := rat * gsl.SqrtMinFloat64
		Ik := gsl.SqrtMinFloat64

		for k := n; k >= 1; k-- {
			Ikm1 := Ikp1 + 2.0*float64(k)/ax*Ik
			Ikp1 = Ik
			Ik = Ikm1
		}

		result.val = I0_scaled.val * (gsl.SqrtMinFloat64 / Ik)
		result.err = I0_scaled.err * (gsl.SqrtMinFloat64 / Ik)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		if x < 0.0 && gsl.IsOdd(n) {
			result.val = -result.val
		}

		return err.ErrorSelect(stat_I0, stat_CF1)
	} else if math.Min(0.29/(float64(n)*float64(n)), 0.5/(float64(n)*float64(n)+x*x)) < 0.5*gsl.Root3Float64Eps {
		stat_as := bessel_Inu_scaled_asymp_unif_e(float64(n), ax, result)

		if x < 0.0 && gsl.IsOdd(n) {
			result.val = -result.val
		}

		return stat_as
	}

	nhi := 2 + (1.2 / gsl.Root6Float64Eps)
	r_Ikp1, r_Ik := new(Result), new(Result)
	stat_a1 := bessel_Inu_scaled_asymp_unif_e(nhi+1.0, ax, r_Ikp1)
	stat_a2 := bessel_Inu_scaled_asymp_unif_e(nhi, ax, r_Ik)

	Ikp1 := r_Ikp1.val
	Ik := r_Ik.val

	for k := int(nhi); k > n; k-- {
		Ikm1 := Ikp1 + 2.0*float64(k)/ax*Ik
		Ikp1 = Ik
		Ik = Ikm1
	}

	result.val = Ik
	result.err = Ik * (r_Ikp1.err/r_Ikp1.val + r_Ik.err/r_Ik.val)

	if x < 0.0 && gsl.IsOdd(n) {
		result.val = -result.val
	}
	return err.ErrorSelect(stat_a1, stat_a2)
}

func Bessel_In_e(n_in int, x float64, result *Result) err.GSLError {
	ax := math.Abs(x)
	n := gsl.AbsInt(n_in)
	In_scaled := new(Result)
	stat_In_scaled := Bessel_In_scaled_e(n, ax, In_scaled)

	if ax > gsl.LnMaxFloat64-1.0 {
		return OverflowError(result)
	}

	ex := math.Exp(ax)
	result.val = ex * In_scaled.val
	result.err = ex * In_scaled.err
	result.err += ax * gsl.Float64Eps * math.Abs(result.val)
	if x < 0.0 && gsl.IsOdd(n) {
		result.val = -result.val
	}

	return stat_In_scaled

}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_In_scaled(n int, x float64) float64 {
	result := new(Result)
	status := Bessel_In_scaled_e(n, x, result)
	return EvalResult(result, status)
}

func Bessel_In(n int, x float64) float64 {
	result := new(Result)
	status := Bessel_In_e(n, x, result)
	return EvalResult(result, status)
}
