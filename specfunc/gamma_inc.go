/* specfunc/gamma_inc.c
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

func gamma_inc_D(a, x float64, result *Result) err.GSLError {
	if a < 10.0 {
		lg := new(Result)
		Lngamma_e(a+1.0, lg)
		lnr := a*math.Log(x) - x - lg.val
		result.val = math.Exp(lnr)
		result.err = 2.0 * gsl.Float64Eps * (math.Abs(lnr) + 1.0) * math.Abs(result.val)
		return nil
	}

	gstar, ln_term := new(Result), new(Result)

	if x < 0.5*a {
		u := x / a
		ln_u := math.Log(u)
		ln_term.val = ln_u - u + 1.0
		ln_term.err = (math.Abs(ln_u) + math.Abs(u) + 1.0) * gsl.Float64Eps

	} else {
		mu := (x - a) / a
		Log_1plusx_mx_e(mu, ln_term)
		ln_term.err += gsl.Float64Eps * math.Abs(mu)
	}

	Gammastar_e(a, gstar)
	term1 := math.Exp(a*ln_term.val) / math.Sqrt(2.0*gsl.Pi*a)
	result.val = term1 / gstar.val
	result.err = 2.0 * gsl.Float64Eps * (math.Abs(a*ln_term.val) + 1.0) * math.Abs(result.val)
	/* Include propagated error from log term */
	result.err += math.Abs(a) * ln_term.err * math.Abs(result.val)
	result.err += gstar.err / math.Abs(gstar.val) * math.Abs(result.val)
	return nil

}

func gamma_inc_P_series(a, x float64, result *Result) err.GSLError {
	nmax := 10000
	D := new(Result)
	stat_D := gamma_inc_D(a, x, D)
	/* Approximating the terms of the series using Stirling's
	   approximation gives t_n = (x/a)^n * exp(-n(n+1)/(2a)), so the
	   convergence condition is n^2 / (2a) + (1-(x/a) + (1/2a)) n >>
	   -log(GSL_DBL_EPS) if we want t_n < O(1e-16) t_0. The condition
	   below detects cases where the minimum value of n is > 5000 */

	if x > 0.995*a && a > 1e5 { /* Difficult case: try continued fraction */
		cf_res := new(Result)
		status := Exprel_n_CF_e(a, x, cf_res)
		result.val = D.val * cf_res.val
		result.err = math.Abs(D.val*cf_res.err) + math.Abs(D.err*cf_res.val)
		return status
	}

	/* Series would require excessive number of terms */

	if x > (a + float64(nmax)) {
		return err.ERROR("gamma_inc_P_series x>>a exceeds range", err.EMAXITER)
	}

	/* Normal case: sum the series */

	sum := 1.0
	term := 1.0
	var remainder float64
	var n int

	/* Handle lower part of the series where t_n is increasing, |x| > a+n */

	var nlow int

	if x > a {
		nlow = int(x) - int(a)
	}

	for n = 1; n < nlow; n++ {
		term *= x / (a + float64(n))
		sum += term
	}

	/* Handle upper part of the series where t_n is decreasing, |x| < a+n */

	for ; n < nmax; n++ {
		term *= x / (a + float64(n))
		sum += term
		if math.Abs(term/sum) < gsl.Float64Eps {
			break
		}
	}

	/*  Estimate remainder of series ~ t_(n+1)/(1-x/(a+n+1)) */

	tnp1 := (x / (a + float64(n))) * term
	remainder = tnp1 / (1.0 - x/(a+float64(n)+1.0))

	result.val = D.val * sum
	result.err = D.err*math.Abs(sum) + math.Abs(D.val*remainder)
	result.err += (1.0 + float64(n)) * gsl.Float64Eps * math.Abs(result.val)

	if n == nmax && math.Abs(remainder/sum) > gsl.SqrtFloat64Eps {
		return err.ERROR("gamma_inc_P_series failed to converge", err.EMAXITER)
	}

	return stat_D

}

/* Q large x asymptotic
 */
func gamma_inc_Q_large_x(a, x float64, result *Result) err.GSLError {
	nmax := 5000
	D := new(Result)
	stat_D := gamma_inc_D(a, x, D)

	sum := 1.0
	term := 1.0
	last := 1.0
	var n int
	for n = 1; n < nmax; n++ {
		term *= (a - float64(n)) / x
		if math.Abs(term/last) > 1.0 {
			break
		}

		if math.Abs(term/sum) < gsl.Float64Eps {
			break
		}
		sum += term
		last = term
	}

	result.val = D.val * (a / x) * sum
	result.err = D.err * math.Abs((a/x)*sum)
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	if n == nmax {
		return err.ERROR("error in large x asymptotic", err.EMAXITER)
	}

	return stat_D
}

/* Uniform asymptotic for x near a, a and x large.
 * See [Temme, p. 285]
 */
func gamma_inc_Q_asymp_unif(a, x float64, result *Result) err.GSLError {
	rta := math.Sqrt(a)
	eps := (x - a) / a
	ln_term := new(Result)
	stat_ln := Log_1plusx_mx_e(eps, ln_term)
	eta := gsl.Sign(eps) * math.Sqrt(-2.0*ln_term.val)

	erfc := new(Result)
	var c0, c1 float64

	/* This used to say erfc(eta*M_SQRT2*rta), which is wrong.
	 * The sqrt(2) is in the denominator. Oops.
	 * Fixed: [GJ] Mon Nov 15 13:25:32 MST 2004
	 */
	Erfc_e(eta*rta/math.Sqrt2, erfc)

	if math.Abs(eps) < gsl.Root5Float64Eps {
		c0 = -1.0/3.0 + eps*(1.0/12.0-eps*(23.0/540.0-eps*(353.0/12960.0-eps*589.0/30240.0)))
		c1 = -1.0/540.0 - eps/288.0
	} else {
		rt_term := math.Sqrt(-2.0 * ln_term.val / (eps * eps))
		lam := x / a
		c0 = (1.0 - 1.0/rt_term) / eps
		c1 = -(eta*eta*eta*(lam*lam+10.0*lam+1.0) - 12.0*eps*eps*eps) / (12.0 * eta * eta * eta * eps * eps * eps)
	}

	R := math.Exp(-0.5*a*eta*eta) / (math.Sqrt2 * math.SqrtPi * rta) * (c0 + c1/a)

	result.val = 0.5*erfc.val + R
	result.err = gsl.Float64Eps*math.Abs(R*0.5*a*eta*eta) + 0.5*erfc.err
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	return stat_ln
}

/* Continued fraction which occurs in evaluation
 * of Q(a,x) or Gamma(a,x).
 *
 *              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
 *   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
 *             1 +   1 +     1 +   1 +      1 +   1 +
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
 *
 * Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
 * See gamma_inc_Q_CF() below.
 *
 */
func gamma_inc_F_CF(a, x float64, result *Result) err.GSLError {
	nmax := 5000
	small := gsl.Float64Eps * gsl.Float64Eps * gsl.Float64Eps

	hn := 1.0
	Cn := 1.0 / small
	Dn := 1.0
	var n int

	for n = 2; n < nmax; n++ {
		var an, delta float64

		if gsl.IsOdd(n) {
			an = 0.5 * (float64(n) - 1) / x
		} else {
			an = (0.5*float64(n) - a) / x
		}

		Dn = 1.0 + an*Dn
		if math.Abs(Dn) < small {
			Dn = small
		}
		Cn = 1.0 + an/Cn
		if math.Abs(Cn) < small {
			Cn = small
		}
		Dn = 1.0 / Dn
		delta = Cn * Dn
		hn *= delta
		if math.Abs(delta-1.0) < gsl.Float64Eps {
			break
		}
	}

	result.val = hn
	result.err = 2.0 * gsl.Float64Eps * math.Abs(hn)
	result.err += gsl.Float64Eps * (2.0 + 0.5*float64(n)) * math.Abs(result.val)

	if n == nmax {
		err.ERROR("error in CF for F(a,x)", err.EMAXITER)
	}

	return nil
}

/* Continued fraction for Q.
 *
 * Q(a,x) = D(a,x) a/x F(a,x)
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no):
 *
 * Since the Gautschi equivalent series method for CF evaluation may lead
 * to singularities, I have replaced it with the modified Lentz algorithm
 * given in
 *
 * I J Thompson and A R Barnett
 * Coulomb and Bessel Functions of Complex Arguments and Order
 * J Computational Physics 64:490-509 (1986)
 *
 * In consequence, gamma_inc_Q_CF_protected() is now obsolete and has been
 * removed.
 *
 * Identification of terms between the above equation for F(a, x) and
 * the first equation in the appendix of Thompson&Barnett is as follows:
 *
 *    b_0 = 0, b_n = 1 for all n > 0
 *
 *    a_1 = 1
 *    a_n = (n/2-a)/x    for n even
 *    a_n = (n-1)/(2x)   for n odd
 *
 */
func gamma_inc_Q_CF(a, x float64, result *Result) err.GSLError {
	D, F := new(Result), new(Result)
	stat_D := gamma_inc_D(a, x, D)
	stat_F := gamma_inc_F_CF(a, x, F)

	result.val = D.val * (a / x) * F.val
	result.err = D.err*math.Abs((a/x)*F.val) + math.Abs(D.val*a/x*F.err)

	return err.ErrorSelect(stat_F, stat_D)
}

func gamma_inc_Q_series(a, x float64, result *Result) err.GSLError {
	var term1, sum, term2 float64
	var stat_sum err.GSLError

	pg21 := -2.404113806319188570799476
	lnx := math.Log(x)
	el := gsl.Euler + lnx
	c1 := -el
	c2 := gsl.Pi*gsl.Pi/12.0 - 0.5*el*el
	c3 := el*(gsl.Pi*gsl.Pi/12.0-el*el/6.0) + pg21/6.0
	c4 := -0.04166666666666666667 * (-1.758243446661483480 + lnx) * (-0.764428657272716373 + lnx) * (0.723980571623507657 + lnx) * (4.107554191916823640 + lnx)
	c5 := -0.0083333333333333333 * (-2.06563396085715900 + lnx) * (-1.28459889470864700 + lnx) * (-0.27583535756454143 + lnx) * (1.33677371336239618 + lnx) * (5.17537282427561550 + lnx)
	c6 := -0.0013888888888888889 * (-2.30814336454783200 + lnx) * (-1.65846557706987300 + lnx) * (-0.88768082560020400 + lnx) * (0.17043847751371778 + lnx) * (1.92135970115863890 + lnx) * (6.22578557795474900 + lnx)
	c7 := -0.00019841269841269841 * (-2.5078657901291800 + lnx) * (-1.9478900888958200 + lnx) * (-1.3194837322612730 + lnx) * (-0.5281322700249279 + lnx) * (0.5913834939078759 + lnx) * (2.4876819633378140 + lnx) * (7.2648160783762400 + lnx)
	c8 := -0.00002480158730158730 * (-2.677341544966400 + lnx) * (-2.182810448271700 + lnx) * (-1.649350342277400 + lnx) * (-1.014099048290790 + lnx) * (-0.191366955370652 + lnx) * (0.995403817918724 + lnx) * (3.041323283529310 + lnx) * (8.295966556941250 + lnx)
	c9 := -2.75573192239859e-6 * (-2.8243487670469080 + lnx) * (-2.3798494322701120 + lnx) * (-1.9143674728689960 + lnx) * (-1.3814529102920370 + lnx) * (-0.7294312810261694 + lnx) * (0.1299079285269565 + lnx) * (1.3873333251885240 + lnx) * (3.5857258865210760 + lnx) * (9.3214237073814600 + lnx)
	c10 := -2.75573192239859e-7 * (-2.9540329644556910 + lnx) * (-2.5491366926991850 + lnx) * (-2.1348279229279880 + lnx) * (-1.6741881076349450 + lnx) * (-1.1325949616098420 + lnx) * (-0.4590034650618494 + lnx) * (0.4399352987435699 + lnx) * (1.7702236517651670 + lnx) * (4.1231539047474080 + lnx) * (10.342627908148680 + lnx)

	term1 = a * (c1 + a*(c2+a*(c3+a*(c4+a*(c5+a*(c6+a*(c7+a*(c8+a*(c9+a*c10)))))))))

	nmax := 5000
	t := 1.0
	var n int
	sum = 1.0

	for n = 1; n < nmax; n++ {
		t *= -x / (float64(n) + 1.0)
		sum += (a + 1.0) / (a + float64(n) + 1.0) * t
		if math.Abs(t/sum) < gsl.Float64Eps {
			break
		}
	}

	if n == nmax {
		stat_sum = err.MaxIteration()
	}

	term2 = (1.0 - term1) * a / (a + 1.0) * x * sum
	result.val = term1 + term2
	result.err = gsl.Float64Eps * (math.Abs(term1) + 2.0*math.Abs(term2))
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return stat_sum
}

/* series for small a and x, but not defined for a == 0 */
func gamma_inc_series(a, x float64, result *Result) err.GSLError {
	Q, G := new(Result), new(Result)

	stat_Q := gamma_inc_Q_series(a, x, Q)

	stat_G := Gamma_e(a, G)

	result.val = Q.val * G.val
	result.err = math.Abs(Q.val*G.err) + math.Abs(Q.err*G.val)
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	return err.ErrorSelect(stat_Q, stat_G)
}

func gamma_inc_a_gt_0(a, x float64, result *Result) err.GSLError {
	/* x > 0 and a > 0; use result for Q */
	Q, G := new(Result), new(Result)
	stat_Q := Gamma_inc_Q_e(a, x, Q)
	stat_G := Gamma_e(a, G)

	result.val = G.val * Q.val
	result.err = math.Abs(G.val*Q.err) + math.Abs(G.err*Q.val)
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	return err.ErrorSelect(stat_G, stat_Q)
}

func gamma_inc_CF(a, x float64, result *Result) err.GSLError {
	F, pre := new(Result), new(Result)
	am1lgx := (a - 1.0) * math.Log(x)
	stat_F := gamma_inc_F_CF(a, x, F)
	stat_E := Exp_err_e(am1lgx-x, gsl.Float64Eps*math.Abs(am1lgx), pre)

	result.val = F.val * pre.val
	result.err = math.Abs(F.err*pre.val) + math.Abs(F.val*pre.err)
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	return err.ErrorSelect(stat_F, stat_E)
}

/* evaluate Gamma(0,x), x > 0 */
func gamma_Inc_A0(x float64, result *Result) err.GSLError {
	return Expint_E1_e(x, result)
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

func Gamma_inc_Q_e(a, x float64, result *Result) err.GSLError {
	if a < 0.0 || x < 0.0 {
		return DomainError(result)
	} else if x == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if a == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if x <= 0.5*a {
		/* If the series is quick, do that. It is
		 * robust and simple.
		 */
		P := new(Result)
		stat_P := gamma_inc_P_series(a, x, P)
		result.val = 1.0 - P.val
		result.err = P.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_P
	} else if a >= 1.0e+06 && (x-a)*(x-a) < a {
		/* Then try the difficult asymptotic regime.
		 * This is the only way to do this region.
		 */
		return gamma_inc_Q_asymp_unif(a, x, result)
	} else if a < 0.2 && x < 5.0 {
		/* Cancellations at small a must be handled
		 * analytically; x should not be too big
		 * either since the series terms grow
		 * with x and log(x).
		 */
		return gamma_inc_Q_series(a, x, result)
	} else if a <= x {
		if x <= 1.0e+06 {
			/* Continued fraction is excellent for x >~ a.
			 * We do not let x be too large when x > a since
			 * it is somewhat pointless to try this there;
			 * the function is rapidly decreasing for
			 * x large and x > a, and it will just
			 * underflow in that region anyway. We
			 * catch that case in the standard
			 * large-x method.
			 */
			return gamma_inc_Q_CF(a, x, result)
		}

		return gamma_inc_Q_large_x(a, x, result)
	}

	if x > a-math.Sqrt(a) {
		/* Continued fraction again. The convergence
		 * is a little slower here, but that is fine.
		 * We have to trade that off against the slow
		 * convergence of the series, which is the
		 * only other option.
		 */
		return gamma_inc_Q_CF(a, x, result)
	}

	P := new(Result)
	stat_P := gamma_inc_P_series(a, x, P)
	result.val = 1.0 - P.val
	result.err = P.err
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return stat_P

}

func Gamma_inc_P_e(a, x float64, result *Result) err.GSLError {
	if a <= 0.0 || x < 0.0 {
		return DomainError(result)
	} else if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if x < 20.0 || x < 0.5*a {
		return gamma_inc_P_series(a, x, result)
	} else if a > 1.0e+06 && (x-a)*(x-a) < a {
		/* Crossover region. Note that Q and P are
		 * roughly the same order of magnitude here,
		 * so the subtraction is stable.
		 */
		Q := new(Result)
		stat_Q := gamma_inc_Q_asymp_unif(a, x, Q)
		result.val = 1.0 - Q.val
		result.err = Q.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_Q
	} else if a <= x {
		/* Q <~ P in this area, so the
		 * subtractions are stable.
		 */
		Q := new(Result)
		var stat_Q err.GSLError
		if a > 0.2*x {
			stat_Q = gamma_inc_Q_CF(a, x, Q)
		} else {
			stat_Q = gamma_inc_Q_large_x(a, x, Q)
		}
		result.val = 1.0 - Q.val
		result.err = Q.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_Q
	}

	if (x-a)*(x-a) < a {
		/* This condition is meant to insure
		 * that Q is not very close to 1,
		 * so the subtraction is stable.
		 */
		Q := new(Result)
		stat_Q := gamma_inc_Q_CF(a, x, Q)
		result.val = 1.0 - Q.val
		result.err = Q.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_Q
	}

	return gamma_inc_P_series(a, x, result)

}

func Gamma_inc_e(a, x float64, result *Result) err.GSLError {

	if x < 0.0 {
		return DomainError(result)
	} else if x == 0.0 {
		return Gamma_e(a, result)
	} else if a == 0.0 {
		return gamma_Inc_A0(x, result)
	} else if a > 0.0 {
		return gamma_inc_a_gt_0(a, x, result)
	} else if x > 0.25 {
		/* continued fraction seems to fail for x too small; otherwise
		   it is ok, independent of the value of |x/a|, because of the
		   non-oscillation in the expansion, i.e. the CF is
		   un-conditionally convergent for a < 0 and x > 0
		*/
		return gamma_inc_CF(a, x, result)
	} else if math.Abs(a) < 0.5 {
		return gamma_inc_series(a, x, result)
	}

	fa := math.Floor(a)
	da := a - fa
	g_da := new(Result)
	var stat_g_da err.GSLError

	if da > 0.0 {
		stat_g_da = gamma_inc_a_gt_0(da, x, g_da)
	} else {
		stat_g_da = gamma_Inc_A0(x, g_da)
	}

	alpha := da
	gax := g_da.val

	for ok := true; ok; ok = alpha > a {
		shift := math.Exp(-x + (alpha-1.0)*math.Log(x))
		gax = (gax - shift) / (alpha - 1.0)
		alpha -= 1.0
	}

	result.val = gax
	result.err = 2.0 * (1.0 + math.Abs(a)) * gsl.Float64Eps * math.Abs(gax)
	return stat_g_da

}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Gamma_inc_P(a, x float64) float64 {
	result := new(Result)
	status := Gamma_inc_P_e(a, x, result)
	return EvalResult(result, status)
}

func Gamma_inc_Q(a, x float64) float64 {
	result := new(Result)
	status := Gamma_inc_Q_e(a, x, result)
	return EvalResult(result, status)
}

func Gamma_inc(a, x float64) float64 {
	result := new(Result)
	status := Gamma_inc_e(a, x, result)
	return EvalResult(result, status)
}
