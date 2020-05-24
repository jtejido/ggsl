/* specfunc/bessel_temme.c
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

/* Calculate series for Y_nu and K_nu for small x and nu.
 * This is applicable for x < 2 and |nu|<=1/2.
 * These functions assume x > 0.
 */

var (
	g1_dat = []float64{
		-1.14516408366268311786898152867,
		0.00636085311347084238122955495,
		0.00186245193007206848934643657,
		0.000152833085873453507081227824,
		0.000017017464011802038795324732,
		-6.4597502923347254354668326451e-07,
		-5.1819848432519380894104312968e-08,
		4.5189092894858183051123180797e-10,
		3.2433227371020873043666259180e-11,
		6.8309434024947522875432400828e-13,
		2.8353502755172101513119628130e-14,
		-7.9883905769323592875638087541e-16,
		-3.3726677300771949833341213457e-17,
		-3.6586334809210520744054437104e-20,
	}

	g1b_cs = &chebyshevSeries{
		g1_dat,
		13,
		-1, 1,
		7,
	}

	g2_dat = []float64{
		1.882645524949671835019616975350,
		-0.077490658396167518329547945212,
		-0.018256714847324929419579340950,
		0.0006338030209074895795923971731,
		0.0000762290543508729021194461175,
		-9.5501647561720443519853993526e-07,
		-8.8927268107886351912431512955e-08,
		-1.9521334772319613740511880132e-09,
		-9.4003052735885162111769579771e-11,
		4.6875133849532393179290879101e-12,
		2.2658535746925759582447545145e-13,
		-1.1725509698488015111878735251e-15,
		-7.0441338200245222530843155877e-17,
		-2.4377878310107693650659740228e-18,
		-7.5225243218253901727164675011e-20,
	}

	g2b_cs = &chebyshevSeries{
		g2_dat,
		14,
		-1, 1,
		8,
	}
)

func temme_gamma(nu float64, g_1pnu, g_1mnu, g1, g2 *float64) err.GSLError {
	anu := math.Abs(nu)
	x := 4.0*anu - 1.0
	r_g1, r_g2 := new(Result), new(Result)
	g1b_cs.Evaluate(x, r_g1)
	g2b_cs.Evaluate(x, r_g2)
	*g1 = r_g1.val
	*g2 = r_g2.val
	*g_1mnu = 1.0 / (r_g2.val + nu*r_g1.val)
	*g_1pnu = 1.0 / (r_g2.val - nu*r_g1.val)
	return nil
}

func bessel_Y_temme(nu, x float64, Ynu, Ynup1 *Result) err.GSLError {
	max_iter := 15000
	half_x := 0.5 * x
	ln_half_x := math.Log(half_x)
	half_x_nu := math.Exp(nu * ln_half_x)
	pi_nu := gsl.Pi * nu
	alpha := pi_nu / 2.0
	sigma := -nu * ln_half_x
	sinrat := 1.0
	if math.Abs(pi_nu) >= gsl.Float64Eps {
		sinrat = pi_nu / math.Sin(pi_nu)
	}
	sinhrat := 1.0
	if math.Abs(sigma) >= gsl.Float64Eps {
		sinhrat = math.Sinh(sigma) / sigma
	}
	sinhalf := 1.0
	if math.Abs(alpha) >= gsl.Float64Eps {
		sinhalf = math.Sin(alpha) / alpha
	}

	sin_sqr := nu * gsl.Pi * gsl.Pi * 0.5 * sinhalf * sinhalf

	var kk int
	var stat_iter err.GSLError
	var g_1pnu, g_1mnu, g1, g2 float64

	stat_g := temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2)

	fk := 2.0 / gsl.Pi * sinrat * (math.Cosh(sigma)*g1 - sinhrat*ln_half_x*g2)
	pk := 1.0 / gsl.Pi / half_x_nu * g_1pnu
	qk := 1.0 / gsl.Pi * half_x_nu * g_1mnu
	hk := pk
	ck := 1.0

	sum0 := fk + sin_sqr*qk
	sum1 := pk

	for kk < max_iter {
		kk++
		k := float64(kk)
		fk = (k*fk + pk + qk) / (k*k - nu*nu)
		ck *= -half_x * half_x / k
		pk /= (k - nu)
		qk /= (k + nu)
		gk := fk + sin_sqr*qk
		hk = -k*gk + pk
		del0 := ck * gk
		del1 := ck * hk
		sum0 += del0
		sum1 += del1
		if math.Abs(del0) < 0.5*(1.0+math.Abs(sum0))*gsl.Float64Eps {
			break
		}
	}

	Ynu.val = -sum0
	Ynu.err = (2.0 + 0.5*float64(kk)) * gsl.Float64Eps * math.Abs(Ynu.val)
	Ynup1.val = -sum1 * 2.0 / x
	Ynup1.err = (2.0 + 0.5*float64(kk)) * gsl.Float64Eps * math.Abs(Ynup1.val)

	if kk >= max_iter {
		stat_iter = err.MaxIteration()
	}
	return err.ErrorSelect(stat_iter, stat_g)
}

func bessel_K_scaled_temme(nu, x float64, K_nu, K_nup1, Kp_nu *float64) err.GSLError {
	max_iter := 15000
	half_x := 0.5 * x
	ln_half_x := math.Log(half_x)
	half_x_nu := math.Exp(nu * ln_half_x)
	pi_nu := gsl.Pi * nu
	sigma := -nu * ln_half_x
	var sinrat, sinhrat float64

	if math.Abs(pi_nu) < gsl.Float64Eps {
		sinrat = 1.
	} else {
		sinrat = pi_nu / math.Sin(pi_nu)
	}

	if math.Abs(sigma) < gsl.Float64Eps {
		sinhrat = 1.
	} else {
		sinhrat = math.Sinh(sigma) / sigma
	}

	ex := math.Exp(x)

	var sum0, sum1 float64
	var fk, pk, qk, hk, ck float64
	var kk int
	var stat_iter err.GSLError
	var g_1pnu, g_1mnu, g1, g2 float64
	stat_g := temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2)

	fk = sinrat * (math.Cosh(sigma)*g1 - sinhrat*ln_half_x*g2)
	pk = 0.5 / half_x_nu * g_1pnu
	qk = 0.5 * half_x_nu * g_1mnu
	hk = pk
	ck = 1.0
	sum0 = fk
	sum1 = hk
	for kk < max_iter {
		kk++
		k := float64(kk)
		fk = (k*fk + pk + qk) / (k*k - nu*nu)
		ck *= half_x * half_x / k
		pk /= (k - nu)
		qk /= (k + nu)
		hk = -k*fk + pk
		del0 := ck * fk
		del1 := ck * hk
		sum0 += del0
		sum1 += del1
		if math.Abs(del0) < 0.5*math.Abs(sum0)*gsl.Float64Eps {
			break
		}
	}

	*K_nu = sum0 * ex
	*K_nup1 = sum1 * 2.0 / x * ex
	*Kp_nu = -*K_nup1 + nu/x**K_nu

	if kk == max_iter {
		stat_iter = err.MaxIteration()
	}

	return err.ErrorSelect(stat_iter, stat_g)
}
