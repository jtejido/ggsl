/* specfunc/bessel_Kn.c
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

func bessel_Kn_scaled_small_x(n int, x float64, result *Result) err.GSLError {
	var term2 float64
	y := 0.25 * x * x
	ln_x_2 := math.Log(0.5 * x)
	ex := math.Exp(x)
	ln_nm1_fact := new(Result)
	Lnfact_e(uint(n-1), ln_nm1_fact)

	ln_pre1 := -float64(n)*ln_x_2 + ln_nm1_fact.val
	if ln_pre1 > gsl.LnMaxFloat64-3.0 {
		return err.ERROR("error", err.EOVRFLW)
	}

	sum1 := 1.0
	k_term := 1.0
	for k := 1; k <= n-1; k++ {
		k_term *= -y / (float64(k) * (float64(n) - float64(k)))
		sum1 += k_term
	}
	term1 := 0.5 * math.Exp(ln_pre1) * sum1

	pre2 := 0.5 * math.Exp(float64(n)*ln_x_2)
	if pre2 > 0.0 {
		KMAX := 20
		psi_n, npk_fact := new(Result), new(Result)
		yk := 1.0
		k_fact := 1.0
		psi_kp1 := -gsl.Euler
		Psi_int_e(n, psi_n)
		Fact_e(uint(n), npk_fact)
		psi_npkp1 := psi_n.val + 1.0/float64(n)
		sum2 := (psi_kp1 + psi_npkp1 - 2.0*ln_x_2) / npk_fact.val
		for k := 1; k < KMAX; k++ {
			psi_kp1 += 1.0 / float64(k)
			psi_npkp1 += 1.0 / (float64(n) + float64(k))
			k_fact *= float64(k)
			npk_fact.val *= float64(n) + float64(k)
			yk *= y
			k_term = yk * (psi_kp1 + psi_npkp1 - 2.0*ln_x_2) / (k_fact * npk_fact.val)
			sum2 += k_term
		}
		term2 = 1.0
		if gsl.IsOdd(n) {
			term2 = -1.0
		}

		term2 *= pre2 * sum2
	} else {
		term2 = 0.0
	}

	result.val = ex * (term1 + term2)
	result.err = ex * gsl.Float64Eps * (math.Abs(ln_pre1)*math.Abs(term1) + math.Abs(term2))
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	return nil
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_Kn_scaled_e(n int, x float64, result *Result) err.GSLError {
	n = gsl.AbsInt(n)

	if x <= 0.0 {
		return DomainError(result)
	} else if n == 0 {
		return Bessel_K0_scaled_e(x, result)
	} else if n == 1 {
		return Bessel_K1_scaled_e(x, result)
	} else if x <= 5.0 {
		return bessel_Kn_scaled_small_x(n, x, result)
	} else if gsl.Root3Float64Eps*x > 0.25*(float64(n)*float64(n)+1) {
		return bessel_Knu_scaled_asympx_e(float64(n), x, result)
	} else if math.Min(0.29/(float64(n)*float64(n)), 0.5/(float64(n)*float64(n)+x*x)) < gsl.Root3Float64Eps {
		return bessel_Knu_scaled_asymp_unif_e(float64(n), x, result)
	}

	/* Upward recurrence. [Gradshteyn + Ryzhik, 8.471.1] */
	two_over_x := 2.0 / x
	r_b_jm1, r_b_j := new(Result), new(Result)
	stat_0 := Bessel_K0_scaled_e(x, r_b_jm1)
	stat_1 := Bessel_K1_scaled_e(x, r_b_j)
	b_jm1 := r_b_jm1.val
	b_j := r_b_j.val

	for j := 1; j < n; j++ {
		b_jp1 := b_jm1 + float64(j)*two_over_x*b_j
		b_jm1 = b_j
		b_j = b_jp1
	}

	result.val = b_j
	result.err = float64(n) * (math.Abs(b_j) * (math.Abs(r_b_jm1.err/r_b_jm1.val) + math.Abs(r_b_j.err/r_b_j.val)))
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	return err.ErrorSelect(stat_0, stat_1)

}

func Bessel_Kn_e(n int, x float64, result *Result) err.GSLError {
	status := Bessel_Kn_scaled_e(n, x, result)
	ex := math.Exp(-x)
	result.val *= ex
	result.err *= ex
	result.err += x * gsl.Float64Eps * math.Abs(result.val)
	return status
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_Kn_scaled(n int, x float64) float64 {
	result := new(Result)
	status := Bessel_Kn_scaled_e(n, x, result)
	return EvalResult(result, status)
}

func Bessel_Kn(n int, x float64) float64 {
	result := new(Result)
	status := Bessel_Kn_e(n, x, result)
	return EvalResult(result, status)
}
