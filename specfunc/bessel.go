/* specfunc/bessel.c
 *
 * Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003 Gerard Jungman
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

const (
	CubeRoot2_ = 1.25992104989487316476721060728
)

func debye_u1(tpow [16]float64) float64 {
	return (3.0*tpow[1] - 5.0*tpow[3]) / 24.0
}

func debye_u2(tpow [16]float64) float64 {
	return (81.0*tpow[2] - 462.0*tpow[4] + 385.0*tpow[6]) / 1152.0
}

func debye_u3(tpow [16]float64) float64 {
	return (30375.0*tpow[3] - 369603.0*tpow[5] + 765765.0*tpow[7] - 425425.0*tpow[9]) / 414720.0
}

func debye_u4(tpow [16]float64) float64 {
	return (4465125.0*tpow[4] - 94121676.0*tpow[6] + 349922430.0*tpow[8] -
		446185740.0*tpow[10] + 185910725.0*tpow[12]) / 39813120.0
}

func debye_u5(tpow [16]float64) float64 {
	return (1519035525.0*tpow[5] - 49286948607.0*tpow[7] +
		284499769554.0*tpow[9] - 614135872350.0*tpow[11] +
		566098157625.0*tpow[13] - 188699385875.0*tpow[15]) / 6688604160.0
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_IJ_taylor_e(v, x, sign float64, kmax int, threshold float64, result *Result) err.GSLError {
	if v < 0.0 || x < 0.0 {
		return DomainError(result)
	} else if x == 0.0 {
		if v == 0.0 {
			result.val = 1.0
			result.err = 0.0
		} else {
			result.val = 0.0
			result.err = 0.0
		}

		return nil
	}

	prefactor, sum := new(Result), new(Result)
	var stat_pre, stat_sum, stat_mul err.GSLError
	if v == 0.0 {
		prefactor.val = 1.0
		prefactor.err = 0.0
		stat_pre = nil
	} else if v < math.MaxInt32-1 {
		poch_factor, tc_factor := new(Result), new(Result)
		N := int(math.Floor(v + 0.5))
		f := v - float64(N)
		stat_poch := Poch_e(float64(N)+1.0, f, poch_factor)
		stat_tc := Taylorcoeff_e(N, 0.5*x, tc_factor)
		p := math.Pow(0.5*x, f)
		prefactor.val = tc_factor.val * p / poch_factor.val
		prefactor.err = tc_factor.err * p / poch_factor.val
		prefactor.err += math.Abs(prefactor.val) / poch_factor.val * poch_factor.err
		prefactor.err += 2.0 * gsl.Float64Eps * math.Abs(prefactor.val)
		stat_pre = err.ErrorSelect(stat_tc, stat_poch)
	} else {
		lg := new(Result)
		stat_lg := Lngamma_e(v+1.0, lg)
		term1 := v * math.Log(0.5*x)
		term2 := lg.val
		ln_pre := term1 - term2
		ln_pre_err := gsl.Float64Eps*(math.Abs(term1)+math.Abs(term2)) + lg.err
		stat_ex := Exp_err_e(ln_pre, ln_pre_err, prefactor)
		stat_pre = err.ErrorSelect(stat_ex, stat_lg)
	}

	y := sign * 0.25 * x * x
	sumk := 1.0
	term := 1.0
	var k int

	for k = 1; k <= kmax; k++ {
		term *= y / ((v + float64(k)) * float64(k))
		sumk += term
		if math.Abs(term/sumk) < threshold {
			break
		}
	}
	sum.val = sumk
	sum.err = threshold * math.Abs(sumk)
	if k >= kmax {
		stat_sum = err.MaxIteration()
	}

	stat_mul = Multiply_err_e(prefactor.val, prefactor.err, sum.val, sum.err, result)

	return err.ErrorSelect(stat_mul, stat_pre, stat_sum)

}

/* Hankel's Asymptotic Expansion - A&S 9.2.5
 *
 * x >> nu*nu+1
 * error ~ O( ((nu*nu+1)/x)^4 )
 *
 * empirical error analysis:
 *   choose  GSL_ROOT4_MACH_EPS * x > (nu*nu + 1)
 *
 * This is not especially useful. When the argument gets
 * large enough for this to apply, the cos() and sin()
 * start loosing digits. However, this seems inevitable
 * for this particular method.
 *
 * Wed Jun 25 14:39:38 MDT 2003 [GJ]
 * This function was inconsistent since the Q term did not
 * go to relative order eps^2. That's why the error estimate
 * originally given was screwy (it didn't make sense that the
 * "empirical" error was coming out O(eps^3)).
 * With Q to proper order, the error is O(eps^4).
 *
 * Sat Mar 15 05:16:18 GMT 2008 [BG]
 * Extended to use additional terms in the series to gain
 * higher accuracy.
 *
 */
func bessel_Jnu_asympx_e(v, x float64, result *Result) err.GSLError {
	mu := 4.0 * v * v
	chi := x - (0.5*v+0.25)*gsl.Pi

	P := 0.
	Q := 0.

	k := 0.
	t := 1.

	for k < 1000 {
		if k != 0 {
			t *= -(mu - (2*k-1)*(2*k-1)) / (k * (8 * x))
		}

		convP := math.Abs(t) < gsl.Float64Eps*math.Abs(P)
		P += t

		k++

		t *= (mu - (2*k-1)*(2*k-1)) / (k * (8 * x))
		convQ := math.Abs(t) < gsl.Float64Eps*math.Abs(Q)
		Q += t

		/* To preserve the consistency of the series we need to exit
		   when P and Q have the same number of terms */

		if convP && convQ && k > (v/2) {
			break
		}

		k++
	}

	pre := math.Sqrt(2.0 / (gsl.Pi * x))
	c := math.Cos(chi)
	s := math.Sin(chi)

	result.val = pre * (c*P - s*Q)
	result.err = pre * gsl.Float64Eps * (math.Abs(c*P) + math.Abs(s*Q) + math.Abs(t)) * (1 + math.Abs(x))
	return nil
}

func bessel_Ynu_asympx_e(v, x float64, result *Result) err.GSLError {
	var ampl, theta float64
	alpha := x
	beta := -0.5 * v * gsl.Pi
	stat_a := bessel_asymp_Mnu_e(v, x, &ampl)
	stat_t := bessel_asymp_thetanu_corr_e(v, x, &theta)
	sin_alpha := math.Sin(alpha)
	cos_alpha := math.Cos(alpha)
	sin_chi := math.Sin(beta + theta)
	cos_chi := math.Cos(beta + theta)
	sin_term := sin_alpha*cos_chi + sin_chi*cos_alpha
	sin_term_mag := math.Abs(sin_alpha*cos_chi) + math.Abs(sin_chi*cos_alpha)
	result.val = ampl * sin_term
	result.err = math.Abs(ampl) * gsl.Float64Eps * sin_term_mag
	result.err += math.Abs(result.val) * 2.0 * gsl.Float64Eps

	if math.Abs(alpha) > 1.0/gsl.Float64Eps {
		result.err *= 0.5 * math.Abs(alpha)
	} else if math.Abs(alpha) > 1.0/gsl.SqrtFloat64Eps {
		result.err *= 256.0 * math.Abs(alpha) * gsl.SqrtFloat64Eps
	}

	return err.ErrorSelect(stat_t, stat_a)
}

func bessel_Inu_scaled_asympx_e(v, x float64, result *Result) err.GSLError {
	mu := 4.0 * v * v
	mum1 := mu - 1.0
	mum9 := mu - 9.0
	pre := 1.0 / math.Sqrt(2.0*gsl.Pi*x)
	r := mu / x
	result.val = pre * (1.0 - mum1/(8.0*x) + mum1*mum9/(128.0*x*x))
	result.err = 2.0*gsl.Float64Eps*math.Abs(result.val) + pre*math.Abs(0.1*r*r*r)
	return nil
}

func bessel_Knu_scaled_asympx_e(v, x float64, result *Result) err.GSLError {
	mu := 4.0 * v * v
	mum1 := mu - 1.0
	mum9 := mu - 9.0
	pre := math.Sqrt(gsl.Pi / (2.0 * x))
	r := v / x
	result.val = pre * (1.0 + mum1/(8.0*x) + mum1*mum9/(128.0*x*x))
	result.err = 2.0*gsl.Float64Eps*math.Abs(result.val) + pre*math.Abs(0.1*r*r*r)
	return nil
}

/* nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.7]
 *
 * error:
 *   The error has the form u_N(t)/nu^N  where  0 <= t <= 1.
 *   It is not hard to show that |u_N(t)| is small for such t.
 *   We have N=6 here, and |u_6(t)| < 0.025, so the error is clearly
 *   bounded by 0.025/nu^6. This gives the asymptotic bound on nu
 *   seen below as nu ~ 100. For general MACH_EPS it will be
 *                     nu > 0.5 / MACH_EPS^(1/6)
 *   When t is small, the bound is even better because |u_N(t)| vanishes
 *   as t->0. In fact u_N(t) ~ C t^N as t->0, with C ~= 0.1.
 *   We write
 *                     err_N <= min(0.025, C(1/(1+(x/nu)^2))^3) / nu^6
 *   therefore
 *                     min(0.29/nu^2, 0.5/(nu^2+x^2)) < MACH_EPS^{1/3}
 *   and this is the general form.
 *
 * empirical error analysis, assuming 14 digit requirement:
 *   choose   x > 50.000 nu   ==>  nu >   3
 *   choose   x > 10.000 nu   ==>  nu >  15
 *   choose   x >  2.000 nu   ==>  nu >  50
 *   choose   x >  1.000 nu   ==>  nu >  75
 *   choose   x >  0.500 nu   ==>  nu >  80
 *   choose   x >  0.100 nu   ==>  nu >  83
 *
 * This makes sense. For x << nu, the error will be of the form u_N(1)/nu^N,
 * since the polynomial term will be evaluated near t=1, so the bound
 * on nu will become constant for small x. Furthermore, increasing x with
 * nu fixed will decrease the error.
 */
func bessel_Inu_scaled_asymp_unif_e(v, x float64, result *Result) err.GSLError {
	z := x / v
	root_term := math.Hypot(1.0, z)
	pre := 1.0 / math.Sqrt(2.0*gsl.Pi*v*root_term)
	eta := root_term + math.Log(z/(1.0+root_term))

	var ex_arg float64

	if z < 1.0/gsl.Root3Float64Eps {
		ex_arg = v * (-z + eta)
	} else {
		ex_arg = -0.5 * v / z * (1.0 - 1.0/(12.0*z*z))
	}

	ex_result := new(Result)
	stat_ex := Exp_e(ex_arg, ex_result)

	if stat_ex == nil {
		t := 1.0 / root_term
		var tpow [16]float64
		tpow[0] = 1.0
		for i := 1; i < 16; i++ {
			tpow[i] = t * tpow[i-1]
		}

		sum := 1.0 + debye_u1(tpow)/v + debye_u2(tpow)/(v*v) + debye_u3(tpow)/(v*v*v) + debye_u4(tpow)/(v*v*v*v) + debye_u5(tpow)/(v*v*v*v*v)

		result.val = pre * ex_result.val * sum
		result.err = pre * ex_result.val / (v * v * v * v * v * v)
		result.err += pre * ex_result.err * math.Abs(sum)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = 0.0
		result.err = 0.0
		return stat_ex
	}

}

/* nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.8]
 *
 * error:
 *   identical to that above for Inu_scaled
 */
func bessel_Knu_scaled_asymp_unif_e(v, x float64, result *Result) err.GSLError {
	var ex_arg float64
	z := x / v
	root_term := math.Hypot(1.0, z)
	pre := math.Sqrt(gsl.Pi / (2.0 * v * root_term))
	eta := root_term + math.Log(z/(1.0+root_term))
	if z < 1.0/gsl.Root3Float64Eps {
		ex_arg = v * (z - eta)
	} else {
		ex_arg = 0.5 * v / z * (1.0 + 1.0/(12.0*z*z))
	}

	ex_result := new(Result)
	stat_ex := Exp_e(ex_arg, ex_result)

	if stat_ex == nil {
		t := 1.0 / root_term
		var tpow [16]float64
		tpow[0] = 1.0
		for i := 1; i < 16; i++ {
			tpow[i] = t * tpow[i-1]
		}

		sum := 1.0 - debye_u1(tpow)/v + debye_u2(tpow)/(v*v) - debye_u3(tpow)/(v*v*v) + debye_u4(tpow)/(v*v*v*v) - debye_u5(tpow)/(v*v*v*v*v)

		result.val = pre * ex_result.val * sum
		result.err = pre * ex_result.err * math.Abs(sum)
		result.err += pre * ex_result.val / (v * v * v * v * v * v)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = 0.0
		result.err = 0.0
		return stat_ex
	}

}

func bessel_JY_mu_restricted(mu, x float64, Jmu, Jmup1, Ymu, Ymup1 *Result) err.GSLError {

	if x <= 0.0 || math.Abs(mu) > 0.5 {
		Jmu.val = 0.0
		Jmu.err = 0.0
		Jmup1.val = 0.0
		Jmup1.err = 0.0
		Ymu.val = 0.0
		Ymu.err = 0.0
		Ymup1.val = 0.0
		Ymup1.err = 0.0
		return err.ERROR("error", err.EDOM)
	} else if x == 0.0 {
		if mu == 0.0 {
			Jmu.val = 1.0
			Jmu.err = 0.0
		} else {
			Jmu.val = 0.0
			Jmu.err = 0.0
		}
		Jmup1.val = 0.0
		Jmup1.err = 0.0
		Ymu.val = 0.0
		Ymu.err = 0.0
		Ymup1.val = 0.0
		Ymup1.err = 0.0
		return err.ERROR("error", err.EDOM)
	} else {
		var stat_Y, stat_J err.GSLError
		if x < 2.0 {
			Jmup2 := new(Result)
			stat_J1 := Bessel_IJ_taylor_e(mu+1.0, x, -1, 100, gsl.Float64Eps, Jmup1)
			stat_J2 := Bessel_IJ_taylor_e(mu+2.0, x, -1, 100, gsl.Float64Eps, Jmup2)
			c := 2.0 * (mu + 1.0) / x
			Jmu.val = c*Jmup1.val - Jmup2.val
			Jmu.err = c*Jmup1.err + Jmup2.err
			Jmu.err += 2.0 * gsl.Float64Eps * math.Abs(Jmu.val)
			stat_J = err.ErrorSelect(stat_J1, stat_J2)
			stat_Y = bessel_Y_temme(mu, x, Ymu, Ymup1)
			return err.ErrorSelect(stat_J, stat_Y)
		} else if x < 1000.0 {
			var P, Q, J_ratio, J_sgn float64
			stat_CF1 := bessel_J_CF1(mu, x, &J_ratio, &J_sgn)
			stat_CF2 := bessel_JY_steed_CF2(mu, x, &P, &Q)
			Jprime_J_ratio := mu/x - J_ratio
			gamma := (P - Jprime_J_ratio) / Q
			Jmu.val = J_sgn * math.Sqrt(2.0/(gsl.Pi*x)/(Q+gamma*(P-Jprime_J_ratio)))
			Jmu.err = 4.0 * gsl.Float64Eps * math.Abs(Jmu.val)
			Jmup1.val = J_ratio * Jmu.val
			Jmup1.err = math.Abs(J_ratio) * Jmu.err
			Ymu.val = gamma * Jmu.val
			Ymu.err = math.Abs(gamma) * Jmu.err
			Ymup1.val = Ymu.val * (mu/x - P - Q/gamma)
			Ymup1.err = Ymu.err*math.Abs(mu/x-P-Q/gamma) + 4.0*gsl.Float64Eps*math.Abs(Ymup1.val)
			return err.ErrorSelect(stat_CF1, stat_CF2)
		}

		/* Use asymptotics for large argument.
		 */
		stat_J0 := bessel_Jnu_asympx_e(mu, x, Jmu)
		stat_J1 := bessel_Jnu_asympx_e(mu+1.0, x, Jmup1)
		stat_Y0 := bessel_Ynu_asympx_e(mu, x, Ymu)
		stat_Y1 := bessel_Ynu_asympx_e(mu+1.0, x, Ymup1)
		stat_J = err.ErrorSelect(stat_J0, stat_J1)
		stat_Y = err.ErrorSelect(stat_Y0, stat_Y1)
		return err.ErrorSelect(stat_J, stat_Y)
	}

}

func bessel_J_CF1(v, x float64, ratio, sgn *float64) err.GSLError {
	RECUR_BIG := gsl.SqrtMaxFloat64
	RECUR_SMALL := gsl.SqrtMinFloat64
	maxiter := 10000
	n := 1
	Anm2 := 1.0
	Bnm2 := 0.0
	Anm1 := 0.0
	Bnm1 := 1.0
	a1 := x / (2.0 * (v + 1.0))
	An := Anm1 + a1*Anm2
	Bn := Bnm1 + a1*Bnm2

	fn := An / Bn
	dn := a1
	s := 1.0

	for n < maxiter {
		n++
		Anm2 = Anm1
		Bnm2 = Bnm1
		Anm1 = An
		Bnm1 = Bn
		an := -x * x / (4.0 * (v + float64(n) - 1.0) * (v + float64(n)))
		An = Anm1 + an*Anm2
		Bn = Bnm1 + an*Bnm2

		if math.Abs(An) > RECUR_BIG || math.Abs(Bn) > RECUR_BIG {
			An /= RECUR_BIG
			Bn /= RECUR_BIG
			Anm1 /= RECUR_BIG
			Bnm1 /= RECUR_BIG
			Anm2 /= RECUR_BIG
		} else if math.Abs(An) < RECUR_SMALL || math.Abs(Bn) < RECUR_SMALL {
			An /= RECUR_SMALL
			Bn /= RECUR_SMALL
			Anm1 /= RECUR_SMALL
			Bnm1 /= RECUR_SMALL
			Anm2 /= RECUR_SMALL
			Bnm2 /= RECUR_SMALL
		}

		old_fn := fn
		fn = An / Bn
		del := old_fn / fn

		dn = 1.0 / (2.0*(v+float64(n))/x - dn)
		if dn < 0.0 {
			s = -s
		}

		if math.Abs(del-1.0) < 2.0*gsl.Float64Eps {
			break
		}

	}

	/* FIXME: we should return an error term here as well, because the
	   error from this recurrence affects the overall error estimate. */

	*ratio = fn
	*sgn = s

	if n >= maxiter {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

func bessel_I_CF1_ser(v, x float64, ratio *float64) err.GSLError {
	maxk := 20000
	tk := 1.0
	sum := 1.0
	rhok := 0.0
	var k int
	for k = 1; k < maxk; k++ {
		ak := 0.25 * (x / (v + float64(k))) * x / (v + float64(k) + 1.0)
		rhok = -ak * (1.0 + rhok) / (1.0 + ak*(1.0+rhok))
		tk *= rhok
		sum += tk

		if math.Abs(tk/sum) < gsl.Float64Eps {
			break
		}

	}

	*ratio = x / (2.0 * (v + 1.0)) * sum

	if k == maxk {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

func bessel_JY_steed_CF2(v, x float64, P, Q *float64) err.GSLError {
	var i int
	max_iter := 10000
	SMALL := 1.0e-100
	x_inv := 1.0 / x
	a := 0.25 - v*v
	p := -0.5 * x_inv
	q := 1.0
	br := 2.0 * x
	bi := 2.0
	fact := a * x_inv / (p*p + q*q)
	cr := br + q*fact
	ci := bi + p*fact
	den := br*br + bi*bi
	dr := br / den
	di := -bi / den
	dlr := cr*dr - ci*di
	dli := cr*di + ci*dr
	temp := p*dlr - q*dli
	q = p*dli + q*dlr
	p = temp
	for i = 2; i <= max_iter; i++ {

		a += 2 * (float64(i) - 1)
		bi += 2.0
		dr = a*dr + br
		di = a*di + bi
		if math.Abs(dr)+math.Abs(di) < SMALL {
			dr = SMALL
		}
		fact = a / (cr*cr + ci*ci)
		cr = br + cr*fact
		ci = bi - ci*fact
		if math.Abs(cr)+math.Abs(ci) < SMALL {
			cr = SMALL
		}
		den = dr*dr + di*di
		dr /= den
		di /= -den
		dlr = cr*dr - ci*di
		dli = cr*di + ci*dr
		temp = p*dlr - q*dli
		q = p*dli + q*dlr
		p = temp
		if math.Abs(dlr-1.0)+math.Abs(dli) < gsl.Float64Eps {
			break
		}
	}

	*P = p
	*Q = q

	if i == max_iter {
		err.ERROR("error", err.EMAXITER)
	}

	return nil
}

/* Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
 * to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
 *
 * This is unstable for small x; x > 2 is a good cutoff.
 * Also requires |nu| < 1/2.
 */
func bessel_K_scaled_steed_temme_CF2(v, x float64, K_nu, K_nup1, Kp_nu *float64) err.GSLError {
	maxiter := 10000
	bi := 2.0 * (1.0 + x)
	di := 1.0 / bi
	delhi := di
	hi := di

	qi := 0.0
	qip1 := 1.0

	ai := -(0.25 - v*v)
	a1 := ai
	ci := -ai
	Qi := -ai
	i := 1

	s := 1.0 + Qi*delhi

	for i = 2; i <= maxiter; i++ {
		var dels, tmp float64
		ai -= 2.0 * (float64(i) - 1)
		ci = -ai * ci / float64(i)
		tmp = (qi - bi*qip1) / ai
		qi = qip1
		qip1 = tmp
		Qi += ci * qip1
		bi += 2.0
		di = 1.0 / (bi + ai*di)
		delhi = (bi*di - 1.0) * delhi
		hi += delhi
		dels = Qi * delhi
		s += dels
		if math.Abs(dels/s) < gsl.Float64Eps {
			break
		}
	}

	hi *= -a1

	*K_nu = math.Sqrt(gsl.Pi/(2.0*x)) / s
	*K_nup1 = *K_nu * (v + x + 0.5 - hi) / x
	*Kp_nu = -*K_nup1 + v/x**K_nu
	if i == maxiter {
		err.ERROR("error", err.EMAXITER)
	}

	return nil
}

func bessel_cos_pi4_e(y, eps float64, result *Result) err.GSLError {
	sy := math.Sin(y)
	cy := math.Cos(y)
	s := sy + cy
	d := sy - cy
	abs_sum := math.Abs(cy) + math.Abs(sy)
	var seps, ceps float64
	if math.Abs(eps) < gsl.Root5Float64Eps {
		e2 := eps * eps
		seps = eps * (1.0 - e2/6.0*(1.0-e2/20.0))
		ceps = 1.0 - e2/2.0*(1.0-e2/12.0)
	} else {
		seps = math.Sin(eps)
		ceps = math.Cos(eps)
	}
	result.val = (ceps*s - seps*d) / gsl.Sqrt2
	result.err = 2.0 * gsl.Float64Eps * (math.Abs(ceps) + math.Abs(seps)) * abs_sum / gsl.Sqrt2

	/* Try to account for error in evaluation of sin(y), cos(y).
	 * This is a little sticky because we don't really know
	 * how the library routines are doing their argument reduction.
	 * However, we will make a reasonable guess.
	 * FIXME ?
	 */
	if y > 1.0/gsl.Float64Eps {
		result.err *= 0.5 * y
	} else if y > 1.0/gsl.SqrtFloat64Eps {
		result.err *= 256.0 * y * gsl.SqrtFloat64Eps
	}

	return nil
}

func bessel_sin_pi4_e(y, eps float64, result *Result) err.GSLError {
	sy := math.Sin(y)
	cy := math.Cos(y)
	s := sy + cy
	d := sy - cy
	abs_sum := math.Abs(cy) + math.Abs(sy)
	var seps, ceps float64
	if math.Abs(eps) < gsl.Root5Float64Eps {
		e2 := eps * eps
		seps = eps * (1.0 - e2/6.0*(1.0-e2/20.0))
		ceps = 1.0 - e2/2.0*(1.0-e2/12.0)
	} else {
		seps = math.Sin(eps)
		ceps = math.Cos(eps)
	}
	result.val = (ceps*d + seps*s) / gsl.Sqrt2
	result.err = 2.0 * gsl.Float64Eps * (math.Abs(ceps) + math.Abs(seps)) * abs_sum / gsl.Sqrt2

	/* Try to account for error in evaluation of sin(y), cos(y).
	 * See above.
	 * FIXME ?
	 */
	if y > 1.0/gsl.Float64Eps {
		result.err *= 0.5 * y
	} else if y > 1.0/gsl.SqrtFloat64Eps {
		result.err *= 256.0 * y * gsl.SqrtFloat64Eps
	}

	return nil
}
