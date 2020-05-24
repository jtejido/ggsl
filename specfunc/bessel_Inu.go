/* specfunc/bessel_Inu.c
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

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_Inu_scaled_e(nu, x float64, result *Result) err.GSLError {

	if x < 0.0 || nu < 0.0 {
		return DomainError(result)
	} else if x*x < 10.0*(nu+1.0) {
		b := new(Result)
		ex := math.Exp(-x)
		stat := Bessel_IJ_taylor_e(nu, x, 1, 100, gsl.Float64Eps, b)

		result.val = b.val * ex
		result.err = b.err * ex
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat
	} else if 0.5/(nu*nu+x*x) < gsl.Root3Float64Eps {
		return bessel_Inu_scaled_asymp_unif_e(nu, x, result)
	}

	N := int(nu + 0.5)
	mu := nu - float64(N)
	var K_mu, K_mup1, Kp_mu, K_nu, K_nup1, K_num1, I_nu_ratio float64
	var stat_Irat, stat_Kmu err.GSLError

	if x < 2.0 {
		stat_Kmu = bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu)
	} else {
		stat_Kmu = bessel_K_scaled_steed_temme_CF2(mu, x, &K_mu, &K_mup1, &Kp_mu)
	}

	K_nu = K_mu
	K_nup1 = K_mup1

	for n := 0; n < N; n++ {
		K_num1 = K_nu
		K_nu = K_nup1
		K_nup1 = 2.0*(mu+float64(n)+1)/x*K_nu + K_num1
	}

	stat_Irat = bessel_I_CF1_ser(nu, x, &I_nu_ratio)

	/* solve for I_nu */
	result.val = 1.0 / (x * (K_nup1 + I_nu_ratio*K_nu))
	result.err = gsl.Float64Eps * (0.5*float64(N) + 2.0) * math.Abs(result.val)

	return err.ErrorSelect(stat_Kmu, stat_Irat)
}

func Bessel_Inu_e(nu, x float64, result *Result) err.GSLError {
	b := new(Result)
	stat_I := Bessel_Inu_scaled_e(nu, x, b)
	stat_e := Exp_mult_err_e(x, math.Abs(x*gsl.Float64Eps), b.val, b.err, result)
	return err.ErrorSelect(stat_e, stat_I)
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_Inu_scaled(nu, x float64) float64 {
	result := new(Result)
	status := Bessel_Inu_scaled_e(nu, x, result)
	return EvalResult(result, status)
}

func Bessel_Inu(nu, x float64) float64 {
	result := new(Result)
	status := Bessel_Inu_e(nu, x, result)
	return EvalResult(result, status)
}
