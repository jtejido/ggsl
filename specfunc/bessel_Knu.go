/* specfunc/bessel_Knu.c
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
func Bessel_Knu_scaled_e(nu, x float64, result *Result) err.GSLError {
	if x <= 0.0 || nu < 0.0 {
		return DomainError(result)
	}

	result_e10 := new(Result_e10)
	status := Bessel_Knu_scaled_e10_e(nu, x, result_e10)
	status2 := Result_smash_e(result_e10, result)
	return err.ErrorSelect(status, status2)
}

func Bessel_Knu_scaled_e10_e(nu, x float64, result *Result_e10) err.GSLError {
	if x <= 0.0 || nu < 0.0 {
		return DomainError_e10(result)
	}

	N := int(nu + 0.5)
	mu := nu - float64(N) /* -1/2 <= mu <= 1/2 */
	var K_mu, K_mup1, Kp_mu float64
	var K_nu, K_nup1, K_num1 float64
	var e10 int

	if x < 2.0 {

		bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu)
	} else {
		bessel_K_scaled_steed_temme_CF2(mu, x, &K_mu, &K_mup1, &Kp_mu)
	}

	/* recurse forward to obtain K_num1, val */
	K_nu = K_mu
	K_nup1 = K_mup1

	for n := 0; n < N; n++ {
		K_num1 = K_nu
		K_nu = K_nup1
		/* rescale the recurrence to avoid overflow */
		if math.Abs(K_nu) > gsl.SqrtMaxFloat64 {
			p := math.Floor(math.Log(math.Abs(K_nu)) / math.Ln10)
			factor := math.Pow(10.0, p)
			K_num1 /= factor
			K_nu /= factor
			e10 += int(p)
		}
		K_nup1 = 2.0*(mu+float64(n)+1)/x*K_nu + K_num1
	}

	result.val = K_nu
	result.err = 2.0 * gsl.Float64Eps * (float64(N) + 4.0) * math.Abs(result.val)
	result.e10 = e10
	return nil
}

func Bessel_Knu_e(nu, x float64, result *Result) err.GSLError {
	b := new(Result)
	stat_K := Bessel_Knu_scaled_e(nu, x, b)
	stat_e := Exp_mult_err_e(-x, 0.0, b.val, b.err, result)
	return err.ErrorSelect(stat_e, stat_K)
}

func Bessel_lnKnu_e(nu, x float64, result *Result) err.GSLError {
	if x <= 0.0 || nu < 0.0 {
		return DomainError(result)
	} else if nu == 0.0 {
		K_scaled := new(Result)
		/* This cannot underflow, and
		 * it will not throw GSL_EDOM
		 * since that is already checked.
		 */
		Bessel_K0_scaled_e(x, K_scaled)
		result.val = -x + math.Log(math.Abs(K_scaled.val))
		result.err = gsl.Float64Eps*math.Abs(x) + math.Abs(K_scaled.err/K_scaled.val)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < 2.0 && nu > 1.0 {
		/* Make use of the inequality
		 * Knu(x) <= 1/2 (2/x)^nu Gamma(nu),
		 * which follows from the integral representation
		 * [Abramowitz+Stegun, 9.6.23 (2)]. With this
		 * we decide whether or not there is an overflow
		 * problem because x is small.
		 */

		lg_nu := new(Result)
		Lngamma_e(nu, lg_nu)
		ln_bound := -math.Ln2 - nu*math.Log(0.5*x) + lg_nu.val
		if ln_bound > gsl.LnMaxFloat64-20.0 {
			/* x must be very small or nu very large (or both).
			 */
			xi := 0.25 * x * x
			sum := 1.0 - xi/(nu-1.0)
			if nu > 2.0 {
				sum += (xi / (nu - 1.0)) * (xi / (nu - 2.0))
			}
			result.val = ln_bound + math.Log(sum)
			result.err = lg_nu.err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return nil
		}
		/* can drop-through here */
	}

	/* We passed the above tests, so no problem.
	 * Evaluate as usual. Note the possible drop-through
	 * in the above code!
	 */
	K_scaled := new(Result_e10)
	status := Bessel_Knu_scaled_e10_e(nu, x, K_scaled)
	result.val = -x + math.Log(math.Abs(K_scaled.val)) + float64(K_scaled.e10)*math.Ln10
	result.err = gsl.Float64Eps*math.Abs(x) + math.Abs(K_scaled.err/K_scaled.val)
	result.err += gsl.Float64Eps * math.Abs(result.val)
	return status

}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_Knu_scaled(nu, x float64) float64 {
	result := new(Result)
	status := Bessel_Knu_scaled_e(nu, x, result)
	return EvalResult(result, status)
}

func Bessel_Knu(nu, x float64) float64 {
	result := new(Result)
	status := Bessel_Knu_e(nu, x, result)
	return EvalResult(result, status)
}

func Bessel_lnKnu(nu, x float64) float64 {
	result := new(Result)
	status := Bessel_lnKnu_e(nu, x, result)
	return EvalResult(result, status)
}
