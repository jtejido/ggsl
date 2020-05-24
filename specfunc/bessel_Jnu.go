/* specfunc/bessel_Jnu.c
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
func bessel_Jnupos_e(nu, x float64, result *Result) err.GSLError {

	if x == 0.0 {
		if nu == 0.0 {
			result.val = 1.0
			result.err = 0.0
		} else {
			result.val = 0.0
			result.err = 0.0
		}

		return nil
	} else if x*x < 10.0*(nu+1.0) {
		return Bessel_IJ_taylor_e(nu, x, -1, 100, gsl.Float64Eps, result)
	} else if nu > 50.0 {
		return bessel_Jnu_asymp_Olver_e(nu, x, result)
	} else if x > 1000.0 {
		return bessel_Jnu_asympx_e(nu, x, result)
	}

	N := int(nu + 0.5)
	mu := nu - float64(N)
	var Jnup1_Jnu, sgn_Jnu float64
	stat_CF1 := bessel_J_CF1(nu, x, &Jnup1_Jnu, &sgn_Jnu)

	if x < 2.0 {
		Y_mu, Y_mup1 := new(Result), new(Result)
		stat_mu := bessel_Y_temme(mu, x, Y_mu, Y_mup1)
		Ynm1 := Y_mu.val
		Yn := Y_mup1.val
		var Ynp1 float64

		for n := 1; n < N; n++ {
			Ynp1 = 2.0*(mu+float64(n))/x*Yn - Ynm1
			Ynm1 = Yn
			Yn = Ynp1
		}

		result.val = 2.0 / (gsl.Pi * x) / (Jnup1_Jnu*Yn - Ynp1)
		result.err = gsl.Float64Eps * (float64(N) + 2.0) * math.Abs(result.val)
		return err.ErrorSelect(stat_mu, stat_CF1)
	}
	var P, Q float64
	stat_CF2 := bessel_JY_steed_CF2(mu, x, &P, &Q)
	Jnp1 := sgn_Jnu * gsl.SqrtMinFloat64 * Jnup1_Jnu
	Jn := sgn_Jnu * gsl.SqrtMinFloat64

	for n := N; n > 0; n-- {
		Jnm1 := 2.0*(mu+float64(n))/x*Jn - Jnp1
		Jnp1 = Jn
		Jn = Jnm1
	}

	Jmup1_Jmu := Jnp1 / Jn
	sgn_Jmu := gsl.Sign(Jn)
	Jmuprime_Jmu := mu/x - Jmup1_Jmu

	gamma := (P - Jmuprime_Jmu) / Q
	Jmu := sgn_Jmu * math.Sqrt(2.0/(gsl.Pi*x)/(Q+gamma*(P-Jmuprime_Jmu)))

	result.val = Jmu * (sgn_Jnu * gsl.SqrtMinFloat64) / Jn
	result.err = 2.0 * gsl.Float64Eps * (float64(N) + 2.0) * math.Abs(result.val)

	return err.ErrorSelect(stat_CF2, stat_CF1)

}

func Bessel_Jnu_e(nu, x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if nu < 0.0 {
		Jstatus := bessel_Jnupos_e(-nu, x, result)
		Jval := result.val
		Jerr := result.err
		Ystatus := bessel_Ynupos_e(-nu, x, result)
		Yval := result.val
		Yerr := result.err
		sinstatus := Sin_pi_e(nu, result)
		s := result.val
		serr := result.err
		cosstatus := Cos_pi_e(nu, result)
		c := result.val
		cerr := result.err
		result.val = s*Yval + c*Jval
		result.err = math.Abs(c*Yerr) + math.Abs(s*Jerr) + math.Abs(cerr*Yval) + math.Abs(serr*Jval)
		return err.ErrorSelect(Jstatus, Ystatus, sinstatus, cosstatus)
	}

	return bessel_Jnupos_e(nu, x, result)
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

func Bessel_Jnu(nu, x float64) float64 {
	result := new(Result)
	status := Bessel_Jnu_e(nu, x, result)
	return EvalResult(result, status)
}
