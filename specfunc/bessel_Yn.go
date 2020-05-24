/* specfunc/bessel_Yn.c
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

/* assumes n >= 1 */
func bessel_Yn_small_x(n int, x float64, result *Result) err.GSLError {
	y := 0.25 * x * x
	ln_x_2 := math.Log(0.5 * x)
	ln_nm1_fact := new(Result)
	var k_term float64
	var term1, sum1, ln_pre1 float64
	var term2, sum2, pre2 float64

	Lnfact_e(uint(n-1), ln_nm1_fact)

	ln_pre1 = -float64(n)*ln_x_2 + ln_nm1_fact.val
	if ln_pre1 > gsl.LnMaxFloat64-3.0 {
		return err.ERROR("error", err.EOVRFLW)
	}

	sum1 = 1.0
	k_term = 1.0
	for k := 1; k <= n-1; k++ {
		k_term *= y / float64(k*(n-k))
		sum1 += k_term
	}
	term1 = -math.Exp(ln_pre1) * sum1 / gsl.Pi

	pre2 = -math.Exp(float64(n)*ln_x_2) / gsl.Pi
	if math.Abs(pre2) > 0.0 {
		KMAX := 20
		psi_n, npk_fact := new(Result), new(Result)
		yk := 1.0
		k_fact := 1.0
		psi_kp1 := -gsl.Euler
		Psi_int_e(n, psi_n)
		Fact_e(uint(n), npk_fact)
		psi_npkp1 := psi_n.val + 1.0/float64(n)
		sum2 = (psi_kp1 + psi_npkp1 - 2.0*ln_x_2) / npk_fact.val
		for k := 1; k < KMAX; k++ {
			psi_kp1 += 1. / float64(k)
			psi_npkp1 += 1. / float64(n+k)
			k_fact *= float64(k)
			npk_fact.val *= float64(n + k)
			yk *= -y
			k_term = yk * (psi_kp1 + psi_npkp1 - 2.0*ln_x_2) / (k_fact * npk_fact.val)
			sum2 += k_term
		}
		term2 = pre2 * sum2
	} else {
		term2 = 0.0
	}

	result.val = term1 + term2
	result.err = gsl.Float64Eps * (math.Abs(ln_pre1)*math.Abs(term1) + math.Abs(term2))
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	return nil
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_Yn_e(n int, x float64, result *Result) err.GSLError {
	sign := 1.

	if n < 0 {
		/* reduce to case n >= 0 */
		n = -n
		if gsl.IsOdd(n) {
			sign = -1.
		}
	}

	if n == 0 {
		status := Bessel_Y0_e(x, result)
		result.val *= sign
		return status
	} else if n == 1 {
		status := Bessel_Y1_e(x, result)
		result.val *= sign
		return status
	} else {
		if x <= 0.0 {
			return DomainError(result)
		}
		if x < 5.0 {
			status := bessel_Yn_small_x(n, x, result)
			result.val *= sign
			return status
		} else if gsl.Root3Float64Eps*x > (float64(n)*float64(n) + 1.0) {
			status := bessel_Ynu_asympx_e(float64(n), x, result)
			result.val *= sign
			return status
		} else if n > 50 {
			status := bessel_Ynu_asymp_Olver_e(float64(n), x, result)
			result.val *= sign
			return status
		} else {
			two_over_x := 2.0 / x
			r_by, r_bym := new(Result), new(Result)
			stat_1 := Bessel_Y1_e(x, r_by)
			stat_0 := Bessel_Y0_e(x, r_bym)
			bym := r_bym.val
			by := r_by.val

			for j := 1; j < n; j++ {
				byp := float64(j)*two_over_x*by - bym
				bym = by
				by = byp
			}
			result.val = sign * by
			result.err = math.Abs(result.val) * (math.Abs(r_by.err/r_by.val) + math.Abs(r_bym.err/r_bym.val))
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

			return err.ErrorSelect(stat_1, stat_0)
		}
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_Yn(n int, x float64) float64 {
	result := new(Result)
	status := Bessel_Yn_e(n, x, result)
	return EvalResult(result, status)
}
