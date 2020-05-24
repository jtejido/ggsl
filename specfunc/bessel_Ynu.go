/* specfunc/bessel_Ynu.c
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
func bessel_Ynupos_e(nu, x float64, result *Result) err.GSLError {
	if nu > 50.0 {
		return bessel_Ynu_asymp_Olver_e(nu, x, result)
	} else {

		/* -1/2 <= mu <= 1/2 */
		N := int(nu + 0.5)
		mu := nu - float64(N)

		Y_mu, Y_mup1 := new(Result), new(Result)
		var stat_mu err.GSLError

		if x < 2.0 {
			/* Determine Ymu, Ymup1 directly. This is really
			 * an optimization since this case could as well
			 * be handled by a call to gsl_sf_bessel_JY_mu_restricted(),
			 * as below.
			 */
			stat_mu = bessel_Y_temme(mu, x, Y_mu, Y_mup1)
		} else {
			/* Determine Ymu, Ymup1 and Jmu, Jmup1.
			 */
			J_mu, J_mup1 := new(Result), new(Result)
			stat_mu = bessel_JY_mu_restricted(mu, x, J_mu, J_mup1, Y_mu, Y_mup1)
		}

		/* Forward recursion to get Ynu, Ynup1.
		 */
		Ynm1 := Y_mu.val
		Yn := Y_mup1.val
		for n := 1; n <= N; n++ {
			Ynp1 := 2.0*(mu+float64(n))/x*Yn - Ynm1
			Ynm1 = Yn
			Yn = Ynp1
		}

		result.val = Ynm1 /* Y_nu */
		result.err = (float64(N) + 1.0) * math.Abs(Ynm1) * (math.Abs(Y_mu.err/Y_mu.val) + math.Abs(Y_mup1.err/Y_mup1.val))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(Ynm1)

		return stat_mu
	}
}

func Bessel_Ynu_e(nu, x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if nu < 0.0 {
		Jstatus := bessel_Jnupos_e(-nu, x, result)
		Jval := result.val
		Jerr := result.err
		Ystatus := bessel_Ynupos_e(-nu, x, result)
		Yval := result.val
		Yerr := result.err
		/* double s = sin(M_PI*nu), c = cos(M_PI*nu); */
		sinstatus := Sin_pi_e(nu, result)
		s := result.val
		serr := result.err
		cosstatus := Cos_pi_e(nu, result)
		c := result.val
		cerr := result.err
		result.val = c*Yval - s*Jval
		result.err = math.Abs(c*Yerr) + math.Abs(s*Jerr) + math.Abs(cerr*Yval) + math.Abs(serr*Jval)
		return err.ErrorSelect(Jstatus, Ystatus, sinstatus, cosstatus)
	}

	return bessel_Ynupos_e(nu, x, result)
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

func Bessel_Ynu(nu, x float64) float64 {
	result := new(Result)
	status := Bessel_Ynu_e(nu, x, result)
	return EvalResult(result, status)

}
