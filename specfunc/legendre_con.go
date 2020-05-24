/* specfunc/legendre_con.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2010 Brian Gough
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
	root_2OverPi_ = 0.797884560802865355879892
	recurse_large = (1.0e-5 * gsl.MaxFloat64)
	recurse_small = (1.0e+5 * gsl.MinFloat64)
)

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* Continued fraction for f_{ell+1}/f_ell
 * f_ell := P^{-mu-ell}_{-1/2 + I tau}(x),  x < 1.0
 *
 * Uses standard CF method from Temme's book.
 */
func conicalP_negmu_xlt1_CF1(mu float64, ell int, tau, x float64, result *Result) err.GSLError {
	var (
		RECUR_BIG = gsl.SqrtMaxFloat64
		maxiter   = 5000
		n         = 1
		xi        = x / (math.Sqrt(1.0-x) * math.Sqrt(1.0+x))
		Anm2      = 1.0
		Bnm2      = 0.0
		Anm1      = 0.0
		Bnm1      = 1.0
		a1        = 1.0
		b1        = 2.0 * (mu + float64(ell) + 1.0) * xi
		An        = b1*Anm1 + a1*Anm2
		Bn        = b1*Bnm1 + a1*Bnm2
		an, bn    float64
		fn        = An / Bn
	)

	for n < maxiter {
		var old_fn, del float64
		n++
		Anm2 = Anm1
		Bnm2 = Bnm1
		Anm1 = An
		Bnm1 = Bn
		an = tau*tau + (mu-0.5+float64(ell)+float64(n))*(mu-0.5+float64(ell)+float64(n))
		bn = 2.0 * (float64(ell) + mu + float64(n)) * xi
		An = bn*Anm1 + an*Anm2
		Bn = bn*Bnm1 + an*Bnm2

		if math.Abs(An) > RECUR_BIG || math.Abs(Bn) > RECUR_BIG {
			An /= RECUR_BIG
			Bn /= RECUR_BIG
			Anm1 /= RECUR_BIG
			Bnm1 /= RECUR_BIG
			Anm2 /= RECUR_BIG
			Bnm2 /= RECUR_BIG
		}

		old_fn = fn
		fn = An / Bn
		del = old_fn / fn

		if math.Abs(del-1.0) < 2.0*gsl.Float64Eps {
			break
		}
	}

	result.val = fn
	result.err = 4.0 * gsl.Float64Eps * (math.Sqrt(float64(n)) + 1.0) * math.Abs(fn)

	if n >= maxiter {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

/* Continued fraction for f_{ell+1}/f_ell
 * f_ell := P^{-mu-ell}_{-1/2 + I tau}(x),  x >= 1.0
 *
 * Uses Gautschi (Euler) equivalent series.
 */
func conicalP_negmu_xgt1_CF1(mu float64, ell int, tau float64, x float64, result *Result) err.GSLError {
	var (
		maxk  = 20000
		gamma = 1.0 - 1.0/(x*x)
		pre   = math.Sqrt(x-1.0) * math.Sqrt(x+1.0) / (x * (2.0 * (float64(ell) + mu + 1.0)))
		tk    = 1.0
		sum   = 1.0
		rhok  = 0.0
		k     int
	)

	for k = 1; k < maxk; k++ {
		fk := float64(k)
		tlk := 2.0 * (float64(ell) + mu + fk)
		l1k := (float64(ell) + mu - 0.5 + 1.0 + fk)
		ak := -(tau*tau + l1k*l1k) / (tlk * (tlk + 2.0)) * gamma
		rhok = -ak * (1.0 + rhok) / (1.0 + ak*(1.0+rhok))
		tk *= rhok
		sum += tk
		if math.Abs(tk/sum) < gsl.Float64Eps {
			break
		}
	}

	result.val = pre * sum
	result.err = math.Abs(pre * tk)
	result.err += 2.0 * gsl.Float64Eps * (math.Sqrt(float64(k)) + 1.0) * math.Abs(pre*sum)

	if k >= maxk {
		return err.ERROR("error", err.EMAXITER)
	}
	return nil
}

/* Implementation of large negative mu asymptotic
 * [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
 */

func olver_U1(beta2, p float64) float64 {
	return (p - 1.0) / (24.0 * (1.0 + beta2)) * (3.0 + beta2*(2.0+5.0*p*(1.0+p)))
}

func olver_U2(beta2, p float64) float64 {
	var (
		beta4 = beta2 * beta2
		p2    = p * p
		poly1 = 4.0*beta4 + 84.0*beta2 - 63.0
		poly2 = 16.0*beta4 + 90.0*beta2 - 81.0
		poly3 = beta2 * p2 * (97.0*beta2 - 432.0 + 77.0*p*(beta2-6.0) - 385.0*beta2*p2*(1.0+p))
	)
	return (1.0 - p) / (1152.0 * (1.0 + beta2)) * (poly1 + poly2 + poly3)
}

var (
	u3c1 = []float64{-1307.0, -1647.0, 3375.0, 3675.0}
	u3c2 = []float64{29366.0, 35835.0, -252360.0, -272630.0,
		276810.0, 290499.0}
	u3c3 = []float64{-29748.0, -8840.0, 1725295.0, 1767025.0,
		-7313470.0, -754778.0, 6309875.0, 6480045.0}
	u3c4 = []float64{2696.0, -16740.0, -524250.0, -183975.0,
		14670540.0, 14172939.0, -48206730.0, -48461985.0,
		36756720.0, 37182145.0}
	u3c5 = []float64{9136.0, 22480.0, 12760.0,
		-252480.0, -4662165.0, -1705341.0,
		92370135.0, 86244015.0, -263678415.0,
		-260275015.0, 185910725.0, 185910725.0}
)

/* Large negative mu asymptotic
 * P^{-mu}_{-1/2 + I tau}, mu . Inf
 * |x| < 1
 *
 * [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
 */
func ConicalP_xlt1_large_neg_mu_e(mu, tau, x float64, result *Result, ln_multiplier *float64) err.GSLError {
	var (
		beta     = tau / mu
		beta2    = beta * beta
		S        = beta * math.Acos((1.0-beta2)/(1.0+beta2))
		p        = x / math.Sqrt(beta2*(1.0-x*x)+1.0)
		lg_mup1  Result
		lg_stat  = Lngamma_e(mu+1.0, &lg_mup1)
		ln_pre_1 = 0.5*mu*(S-math.Log(1.0+beta2)+math.Log((1.0-p)/(1.0+p))) - lg_mup1.val
		ln_pre_2 = -0.25 * math.Log(1.0+beta2*(1.0-x))
		ln_pre_3 = -tau * math.Atan(p*beta)
		ln_pre   = ln_pre_1 + ln_pre_2 + ln_pre_3
		sum      = 1.0 - olver_U1(beta2, p)/mu + olver_U2(beta2, p)/(mu*mu)
	)

	if sum == 0.0 {
		result.val = 0.0
		result.err = 0.0
		*ln_multiplier = 0.0
		return nil
	} else {
		stat_e := Exp_mult_e(ln_pre, sum, result)
		if stat_e != nil {
			result.val = sum
			result.err = 2.0 * gsl.Float64Eps * math.Abs(sum)
			*ln_multiplier = ln_pre
		} else {
			*ln_multiplier = 0.0
		}
		return lg_stat
	}
}

/* Implementation of large tau asymptotic
 *
 * A_n^{-mu}, B_n^{-mu}  [Olver, p.465, 469]
 */
func olver_B0_xi(mu, xi float64) float64 {
	return (1.0 - 4.0*mu*mu) / (8.0 * xi) * (1.0/math.Tanh(xi) - 1.0/xi)
}

func olver_A1_xi(mu, xi, x float64) float64 {
	B := olver_B0_xi(mu, xi)
	var psi float64
	if math.Abs(x-1.0) < gsl.Root4Float64Eps {
		y := x - 1.0
		s := -1.0/3.0 + y*(2.0/15.0-y*(61.0/945.0-452.0/14175.0*y))
		psi = (4.0*mu*mu - 1.0) / 16.0 * s
	} else {
		psi = (4.0*mu*mu - 1.0) / 16.0 * (1.0/(x*x-1.0) - 1.0/(xi*xi))
	}
	return 0.5*xi*xi*B*B + (mu+0.5)*B - psi + mu/6.0*(0.25-mu*mu)
}

func olver_B0_th(mu, theta float64) float64 {
	return -(1.0 - 4.0*mu*mu) / (8.0 * theta) * (1.0/math.Tan(theta) - 1.0/theta)
}

func olver_A1_th(mu, theta, x float64) float64 {
	B := olver_B0_th(mu, theta)
	var psi float64
	if math.Abs(x-1.0) < gsl.Root4Float64Eps {
		y := 1.0 - x
		s := -1.0/3.0 + y*(2.0/15.0-y*(61.0/945.0-452.0/14175.0*y))
		psi = (4.0*mu*mu - 1.0) / 16.0 * s
	} else {
		psi = (4.0*mu*mu - 1.0) / 16.0 * (1.0/(x*x-1.0) + 1.0/(theta*theta))
	}
	return -0.5*theta*theta*B*B + (mu+0.5)*B - psi + mu/6.0*(0.25-mu*mu)
}

/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}
 * 1 < x
 * tau . Inf
 * [Olver, p. 469]
 */
func ConicalP_xgt1_neg_mu_largetau_e(mu, tau, x, acosh_x float64, result *Result, ln_multiplier *float64) err.GSLError {
	var (
		xi              = acosh_x
		ln_xi_pre       float64
		ln_pre          float64
		sumA, sumB, sum float64
		arg             float64
		J_mup1, J_mu    Result
		J_mum1          float64
	)

	if xi < gsl.Root4Float64Eps {
		ln_xi_pre = -xi * xi / 6.0 /* math.Log(1.0 - xi*xi/6.0) */
	} else {
		var lnshxi Result
		Lnsinh_e(xi, &lnshxi)
		ln_xi_pre = math.Log(xi) - lnshxi.val /* math.Log(xi/sinh(xi) */
	}

	ln_pre = 0.5*ln_xi_pre - mu*math.Log(tau)

	arg = tau * xi

	Bessel_Jnu_e(mu+1.0, arg, &J_mup1)
	Bessel_Jnu_e(mu, arg, &J_mu)
	J_mum1 = -J_mup1.val + 2.0*mu/arg*J_mu.val /* careful of mu < 1 */

	sumA = 1.0 - olver_A1_xi(-mu, xi, x)/(tau*tau)
	sumB = olver_B0_xi(-mu, xi)
	sum = J_mu.val*sumA - xi/tau*J_mum1*sumB

	if sum == 0.0 {
		result.val = 0.0
		result.err = 0.0
		*ln_multiplier = 0.0
		return nil
	} else {
		stat_e := Exp_mult_e(ln_pre, sum, result)
		if stat_e != nil {
			result.val = sum
			result.err = 2.0 * gsl.Float64Eps * math.Abs(sum)
			*ln_multiplier = ln_pre
		} else {
			*ln_multiplier = 0.0
		}
		return nil
	}
}

/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}
 * -1 < x < 1
 * tau . Inf
 * [Olver, p. 473]
 */
func ConicalP_xlt1_neg_mu_largetau_e(mu, tau, x, acos_x float64, result *Result, ln_multiplier *float64) err.GSLError {
	var (
		theta                   = acos_x
		ln_th_pre               float64
		ln_pre                  float64
		sumA, sumB, sum, sumerr float64
		arg                     float64
		I_mup1, I_mu            Result
		I_mum1                  float64
	)

	if theta < gsl.Root4Float64Eps {
		ln_th_pre = theta * theta / 6.0 /* math.Log(1.0 + theta*theta/6.0) */
	} else {
		ln_th_pre = math.Log(theta / math.Sin(theta))
	}

	ln_pre = 0.5*ln_th_pre - mu*math.Log(tau)

	arg = tau * theta
	Bessel_Inu_e(mu+1.0, arg, &I_mup1)
	Bessel_Inu_e(mu, arg, &I_mu)
	I_mum1 = I_mup1.val + 2.0*mu/arg*I_mu.val /* careful of mu < 1 */

	sumA = 1.0 - olver_A1_th(-mu, theta, x)/(tau*tau)
	sumB = olver_B0_th(-mu, theta)
	sum = I_mu.val*sumA - theta/tau*I_mum1*sumB
	sumerr = math.Abs(I_mu.err * sumA)
	sumerr += math.Abs(I_mup1.err * theta / tau * sumB)
	sumerr += math.Abs(I_mu.err * theta / tau * sumB * 2.0 * mu / arg)

	if sum == 0.0 {
		result.val = 0.0
		result.err = 0.0
		*ln_multiplier = 0.0
		return nil
	} else {
		stat_e := Exp_mult_e(ln_pre, sum, result)
		if stat_e != nil {
			result.val = sum
			result.err = sumerr
			result.err += gsl.Float64Eps * math.Abs(sum)
			*ln_multiplier = ln_pre
		} else {
			*ln_multiplier = 0.0
		}
		return nil
	}
}

/* Hypergeometric function which appears in the
 * large x expansion below:
 *
 *   2F1(1/4 - mu/2 - I tau/2, 3/4 - mu/2 - I tau/2, 1 - I tau, y)
 *
 * Note that for the usage below y = 1/x^2;
 */
func conicalP_hyperg_large_x(mu, tau, y float64, reF, imF *float64) err.GSLError {
	var (
		kmax = 1000
		re_a = 0.25 - 0.5*mu
		re_b = 0.75 - 0.5*mu
		re_c = 1.0
		im_a = -0.5 * tau
		im_b = -0.5 * tau
		im_c = -tau

		re_sum  = 1.0
		im_sum  = 0.0
		re_term = 1.0
		im_term = 0.0
		k       int
	)

	for k = 1; k <= kmax; k++ {
		var (
			fk            = float64(k)
			re_ak         = re_a + fk - 1.0
			re_bk         = re_b + fk - 1.0
			re_ck         = re_c + fk - 1.0
			im_ak         = im_a
			im_bk         = im_b
			im_ck         = im_c
			den           = re_ck*re_ck + im_ck*im_ck
			re_multiplier = ((re_ak*re_bk-im_ak*im_bk)*re_ck + im_ck*(im_ak*re_bk+re_ak*im_bk)) / den
			im_multiplier = ((im_ak*re_bk+re_ak*im_bk)*re_ck - im_ck*(re_ak*re_bk-im_ak*im_bk)) / den
			re_tmp        = re_multiplier*re_term - im_multiplier*im_term
			im_tmp        = im_multiplier*re_term + re_multiplier*im_term
			asum          = math.Abs(re_sum) + math.Abs(im_sum)
		)
		re_term = y / fk * re_tmp
		im_term = y / fk * im_tmp
		if math.Abs(re_term/asum) < gsl.Float64Eps && math.Abs(im_term/asum) < gsl.Float64Eps {
			break
		}
		re_sum += re_term
		im_sum += im_term
	}

	*reF = re_sum
	*imF = im_sum

	if k == kmax {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

/* P^{mu}_{-1/2 + I tau}
 * x.Inf
 */
func ConicalP_large_x_e(mu, tau, x float64, result *Result, ln_multiplier *float64) err.GSLError {
	/* 2F1 term
	 */

	var y float64
	if x < 0.5*gsl.SqrtMaxFloat64 {
		y = 1 / (x * x)
	}
	var (
		reF, imF float64
		stat_F   = conicalP_hyperg_large_x(mu, tau, y, &reF, &imF)

		/* f = Gamma(+i tau)/Gamma(1/2 - mu + i tau)
		 * FIXME: shift so it's better for tau. 0
		 */
		lgr_num, lgth_num Result
		lgr_den, lgth_den Result
		stat_gn           = Lngamma_complex_e(0.0, tau, &lgr_num, &lgth_num)
		stat_gd           = Lngamma_complex_e(0.5-mu, tau, &lgr_den, &lgth_den)

		angle = lgth_num.val - lgth_den.val + math.Atan2(imF, reF)

		lnx         = math.Log(x)
		lnxp1       = math.Log(x + 1.0)
		lnxm1       = math.Log(x - 1.0)
		lnpre_const = 0.5*gsl.Ln2 - 0.5*gsl.LnPi
		lnpre_comm  = (mu-0.5)*lnx - 0.5*mu*(lnxp1+lnxm1)
		lnpre_err   = gsl.Float64Eps*(0.5*gsl.Ln2+0.5*gsl.LnPi) + gsl.Float64Eps*math.Abs((mu-0.5)*lnx) + gsl.Float64Eps*math.Abs(0.5*mu)*(math.Abs(lnxp1)+math.Abs(lnxm1))

		/*  result = pre*|F|*|f| * cos(angle - tau * (math.Log(x)+gsl.Ln2))
		 */
		cos_result Result
		stat_cos   = Cos_e(angle+tau*(math.Log(x)+gsl.Ln2), &cos_result)
		status     = err.ErrorSelect(stat_cos, stat_gd, stat_gn, stat_F)
	)

	if cos_result.val == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return status
	} else {
		lnFf_val := 0.5*math.Log(reF*reF+imF*imF) + lgr_num.val - lgr_den.val
		lnFf_err := lgr_num.err + lgr_den.err + gsl.Float64Eps*math.Abs(lnFf_val)
		lnnoc_val := lnpre_const + lnpre_comm + lnFf_val
		lnnoc_err := lnpre_err + lnFf_err + gsl.Float64Eps*math.Abs(lnnoc_val)
		stat_e := Exp_mult_err_e(lnnoc_val, lnnoc_err, cos_result.val, cos_result.err, result)
		if stat_e == nil {
			*ln_multiplier = 0.0
		} else {
			result.val = cos_result.val
			result.err = cos_result.err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			*ln_multiplier = lnnoc_val
		}

		return status
	}
}

/* P^{mu}_{-1/2 + I tau}  first hypergeometric representation
 * -1 < x < 1
 * This is more effective for |x| small, however it will work w/o
 * reservation for any x < 0 because everything is positive
 * definite in that case.
 *
 * [Kolbig,   (3)] (note typo in args of gamma functions)
 * [Bateman, (22)] (correct form)
 */
func conicalP_xlt1_hyperg_A(mu, tau, x float64, result *Result) err.GSLError {
	var (
		x2                           = x * x
		err_amp                      = 1.0 + 1.0/(gsl.Float64Eps+math.Abs(1.0-math.Abs(x)))
		pre_val                      = gsl.SqrtPi / math.Pow(0.5*math.Sqrt(1-x2), mu)
		pre_err                      = err_amp * gsl.Float64Eps * (math.Abs(mu) + 1.0) * math.Abs(pre_val)
		ln_g1, ln_g2, arg_g1, arg_g2 Result
		F1, F2                       Result
		pre1, pre2                   Result
		t1_val, t1_err               float64
		t2_val, t2_err               float64
		stat_F1                      = Hyperg_2F1_conj_e(0.25-0.5*mu, 0.5*tau, 0.5, x2, &F1)
		stat_F2                      = Hyperg_2F1_conj_e(0.75-0.5*mu, 0.5*tau, 1.5, x2, &F2)
		status                       = err.ErrorSelect(stat_F1, stat_F2)
	)

	Lngamma_complex_e(0.75-0.5*mu, -0.5*tau, &ln_g1, &arg_g1)
	Lngamma_complex_e(0.25-0.5*mu, -0.5*tau, &ln_g2, &arg_g2)

	Exp_err_e(-2.0*ln_g1.val, 2.0*ln_g1.err, &pre1)
	Exp_err_e(-2.0*ln_g2.val, 2.0*ln_g2.err, &pre2)
	pre2.val *= -2.0 * x
	pre2.err *= 2.0 * math.Abs(x)
	pre2.err += gsl.Float64Eps * math.Abs(pre2.val)

	t1_val = pre1.val * F1.val
	t1_err = math.Abs(pre1.val)*F1.err + pre1.err*math.Abs(F1.val)
	t2_val = pre2.val * F2.val
	t2_err = math.Abs(pre2.val)*F2.err + pre2.err*math.Abs(F2.val)

	result.val = pre_val * (t1_val + t2_val)
	result.err = pre_val * (t1_err + t2_err)
	result.err += pre_err * math.Abs(t1_val+t2_val)
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	return status
}

/* V0, V1 from Kolbig, m = 0
 */
func conicalP_0_V(t, f, tau, sgn float64, V0, V1 *float64) err.GSLError {
	var (
		C [8]float64
		T [8]float64
		H [8]float64
		V [12]float64
		i int
	)
	T[0] = 1.0
	H[0] = 1.0
	V[0] = 1.0
	for i = 1; i <= 7; i++ {
		T[i] = T[i-1] * t
		H[i] = H[i-1] * (t * f)
	}

	for i = 1; i <= 11; i++ {
		V[i] = V[i-1] * tau
	}

	C[0] = 1.0
	C[1] = (H[1] - 1.0) / (8.0 * T[1])
	C[2] = (9.0*H[2] + 6.0*H[1] - 15.0 - sgn*8.0*T[2]) / (128.0 * T[2])
	C[3] = 5.0 * (15.0*H[3] + 27.0*H[2] + 21.0*H[1] - 63.0 - sgn*T[2]*(16.0*H[1]+24.0)) / (1024.0 * T[3])
	C[4] = 7.0 * (525.0*H[4] + 1500.0*H[3] + 2430.0*H[2] + 1980.0*H[1] - 6435.0 + 192.0*T[4] - sgn*T[2]*(720.0*H[2]+1600.0*H[1]+2160.0)) / (32768.0 * T[4])
	C[5] = 21.0 * (2835.0*H[5] + 11025.0*H[4] + 24750.0*H[3] + 38610.0*H[2] + 32175.0*H[1] - 109395.0 + T[4]*(1984.0*H[1]+4032.0) - sgn*T[2]*(4800.0*H[3]+15120.0*H[2]+26400.0*H[1]+34320.0)) / (262144.0 * T[5])
	C[6] = 11.0 * (218295.0*H[6] + 1071630.0*H[5] + 3009825.0*H[4] + 6142500.0*H[3] + 9398025.0*H[2] + 7936110.0*H[1] - 27776385.0 + T[4]*(254016.0*H[2]+749952.0*H[1]+1100736.0) - sgn*T[2]*(441000.0*H[4]+1814400.0*H[3]+4127760.0*H[2]+6552000.0*H[1]+8353800.0+31232.0*T[4])) / (4194304.0 * T[6])

	*V0 = C[0] + (-4.0*C[3]/T[1]+C[4])/V[4] + (-192.0*C[5]/T[3]+144.0*C[6]/T[2])/V[8] + sgn*(-C[2]/V[2]+(-24.0*C[4]/T[2]+12.0*C[5]/T[1]-C[6])/V[6]+(-1920.0*C[6]/T[4])/V[10])
	*V1 = C[1]/V[1] + (8.0*(C[3]/T[2]-C[4]/T[1])+C[5])/V[5] + (384.0*C[5]/T[4]-768.0*C[6]/T[3])/V[9] + sgn*((2.0*C[2]/T[1]-C[3])/V[3]+(48.0*C[4]/T[3]-72.0*C[5]/T[2]+18.0*C[6]/T[1])/V[7]+(3840.0*C[6]/T[5])/V[11])

	return nil
}

/* V0, V1 from Kolbig, m = 1
 */
func conicalP_1_V(t, f, tau, sgn float64, V0, V1 *float64) err.GSLError {
	var (
		Cm1 float64
		C   [8]float64
		T   [8]float64
		H   [8]float64
		V   [12]float64
		i   int
	)
	T[0] = 1.0
	H[0] = 1.0
	V[0] = 1.0
	for i = 1; i <= 7; i++ {
		T[i] = T[i-1] * t
		H[i] = H[i-1] * (t * f)
	}

	for i = 1; i <= 11; i++ {
		V[i] = V[i-1] * tau
	}

	Cm1 = -1.0
	C[0] = 3.0 * (1.0 - H[1]) / (8.0 * T[1])
	C[1] = (-15.0*H[2] + 6.0*H[1] + 9.0 + sgn*8.0*T[2]) / (128.0 * T[2])
	C[2] = 3.0 * (-35.0*H[3] - 15.0*H[2] + 15.0*H[1] + 35.0 + sgn*T[2]*(32.0*H[1]+8.0)) / (1024.0 * T[3])
	C[3] = (-4725.0*H[4] - 6300.0*H[3] - 3150.0*H[2] + 3780.0*H[1] + 10395.0 - 1216.0*T[4] + sgn*T[2]*(6000.0*H[2]+5760.0*H[1]+1680.0)) / (32768.0 * T[4])
	C[4] = 7.0 * (-10395.0*H[5] - 23625.0*H[4] - 28350.0*H[3] - 14850.0*H[2] + 19305.0*H[1] + 57915.0 - T[4]*(6336.0*H[1]+6080.0) + sgn*T[2]*(16800.0*H[3]+30000.0*H[2]+25920.0*H[1]+7920.0)) / (262144.0 * T[5])
	C[5] = (-2837835.0*H[6] - 9168390.0*H[5] - 16372125.0*H[4] - 18918900*H[3] - 10135125.0*H[2] + 13783770.0*H[1] + 43648605.0 - T[4]*(3044160.0*H[2]+5588352.0*H[1]+4213440.0) + sgn*T[2]*(5556600.0*H[4]+14817600.0*H[3]+20790000.0*H[2]+17297280.0*H[1]+5405400.0+323072.0*T[4])) / (4194304.0 * T[6])
	C[6] = 0.0

	*V0 = C[0] + (-4.0*C[3]/T[1]+C[4])/V[4] + (-192.0*C[5]/T[3]+144.0*C[6]/T[2])/V[8] + sgn*(-C[2]/V[2]+(-24.0*C[4]/T[2]+12.0*C[5]/T[1]-C[6])/V[6])
	*V1 = C[1]/V[1] + (8.0*(C[3]/T[2]-C[4]/T[1])+C[5])/V[5] + (384.0*C[5]/T[4]-768.0*C[6]/T[3])/V[9] + sgn*(Cm1*V[1]+(2.0*C[2]/T[1]-C[3])/V[3]+(48.0*C[4]/T[3]-72.0*C[5]/T[2]+18.0*C[6]/T[1])/V[7])

	return nil
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/* P^0_{-1/2 + I lambda}
 */
func ConicalP_0_e(lambda, x float64, result *Result) err.GSLError {

	if x <= -1.0 {
		return DomainError(result)
	} else if x == 1.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if lambda == 0.0 {
		var K Result
		var stat_K err.GSLError
		if x < 1.0 {
			th := math.Acos(x)
			s := math.Sin(0.5 * th)
			stat_K = Ellint_Kcomp_e(s, gsl.MODE_DEFAULT, &K)
			result.val = 2.0 / gsl.Pi * K.val
			result.err = 2.0 / gsl.Pi * K.err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return stat_K
		} else {
			xi := math.Acosh(x)
			c := math.Cosh(0.5 * xi)
			t := math.Tanh(0.5 * xi)
			stat_K = Ellint_Kcomp_e(t, gsl.MODE_DEFAULT, &K)
			result.val = 2.0 / gsl.Pi / c * K.val
			result.err = 2.0 / gsl.Pi / c * K.err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return stat_K
		}
	} else if (x <= 0.0 && lambda < 1000.0) || (x < 0.1 && lambda < 17.0) || (x < 0.2 && lambda < 5.0) {
		return conicalP_xlt1_hyperg_A(0.0, lambda, x, result)
	} else if (x <= 0.2 && lambda < 17.0) || (x <= 1.5 && lambda < 20.0) {
		return Hyperg_2F1_conj_e(0.5, lambda, 1.0, (1.0-x)/2, result)
	} else if 1.5 < x && lambda < gsl.Max(x, 20.0) {
		var P Result
		var lm float64
		stat_P := ConicalP_large_x_e(0.0, lambda, x, &P, &lm)
		stat_e := Exp_mult_err_e(lm, 2.0*gsl.Float64Eps*math.Abs(lm), P.val, P.err, result)
		return err.ErrorSelect(stat_e, stat_P)
	} else {
		var V0, V1 float64
		if x < 1.0 {
			th := math.Acos(x)
			sth := math.Sqrt(1.0 - x*x) /* sin(th) */
			var I0, I1 Result
			stat_I0 := Bessel_I0_scaled_e(th*lambda, &I0)
			stat_I1 := Bessel_I1_scaled_e(th*lambda, &I1)
			stat_I := err.ErrorSelect(stat_I0, stat_I1)
			stat_V := conicalP_0_V(th, x/sth, lambda, -1.0, &V0, &V1)
			bessterm := V0*I0.val + V1*I1.val
			besserr := math.Abs(V0)*I0.err + math.Abs(V1)*I1.err
			arg1 := th * lambda
			sqts := math.Sqrt(th / sth)
			stat_e := Exp_mult_err_e(arg1, 4.0*gsl.Float64Eps*math.Abs(arg1), sqts*bessterm, sqts*besserr, result)
			return err.ErrorSelect(stat_e, stat_V, stat_I)
		} else {
			sh := math.Sqrt(x-1.0) * math.Sqrt(x+1.0) /* sinh(xi)      */
			xi := math.Log(x + sh)                    /* xi = acosh(x) */
			var J0, J1 Result
			stat_J0 := Bessel_J0_e(xi*lambda, &J0)
			stat_J1 := Bessel_J1_e(xi*lambda, &J1)
			stat_J := err.ErrorSelect(stat_J0, stat_J1)
			stat_V := conicalP_0_V(xi, x/sh, lambda, 1.0, &V0, &V1)
			bessterm := V0*J0.val + V1*J1.val
			besserr := math.Abs(V0)*J0.err + math.Abs(V1)*J1.err
			pre_val := math.Sqrt(xi / sh)
			pre_err := 2.0 * math.Abs(pre_val)
			result.val = pre_val * bessterm
			result.err = pre_val * besserr
			result.err += pre_err * math.Abs(bessterm)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return err.ErrorSelect(stat_V, stat_J)
		}
	}
}

/* P^1_{-1/2 + I lambda}
 */
func ConicalP_1_e(lambda, x float64, result *Result) err.GSLError {
	if x <= -1.0 {
		return DomainError(result)
	} else if lambda == 0.0 {
		var K, E Result
		var stat_K, stat_E err.GSLError
		if x == 1.0 {
			result.val = 0.0
			result.err = 0.0
			return nil
		} else if x < 1.0 {
			if 1.0-x < gsl.SqrtFloat64Eps {
				err_amp := gsl.Max(1.0, 1.0/(gsl.Float64Eps+math.Abs(1.0-x)))
				result.val = 0.25 / gsl.Sqrt2 * math.Sqrt(1.0-x) * (1.0 + 5.0/16.0*(1.0-x))
				result.err = err_amp * 3.0 * gsl.Float64Eps * math.Abs(result.val)
				return nil
			} else {
				var (
					th  = math.Acos(x)
					s   = math.Sin(0.5 * th)
					c2  = 1.0 - s*s
					sth = math.Sin(th)
					pre = 2.0 / (gsl.Pi * sth)
				)
				stat_K = Ellint_Kcomp_e(s, gsl.MODE_DEFAULT, &K)
				stat_E = Ellint_Ecomp_e(s, gsl.MODE_DEFAULT, &E)
				result.val = pre * (E.val - c2*K.val)
				result.err = pre * (E.err + math.Abs(c2)*K.err)
				result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
				return err.ErrorSelect(stat_K, stat_E)
			}
		} else {
			if x-1.0 < gsl.SqrtFloat64Eps {
				err_amp := gsl.Max(1.0, 1.0/(gsl.Float64Eps+math.Abs(1.0-x)))
				result.val = -0.25 / gsl.Sqrt2 * math.Sqrt(x-1.0) * (1.0 - 5.0/16.0*(x-1.0))
				result.err = err_amp * 3.0 * gsl.Float64Eps * math.Abs(result.val)
				return nil
			} else {
				var (
					xi  = math.Acosh(x)
					c   = math.Cosh(0.5 * xi)
					t   = math.Tanh(0.5 * xi)
					sxi = math.Sinh(xi)
					pre = 2.0 / (gsl.Pi * sxi) * c
				)
				stat_K = Ellint_Kcomp_e(t, gsl.MODE_DEFAULT, &K)
				stat_E = Ellint_Ecomp_e(t, gsl.MODE_DEFAULT, &E)
				result.val = pre * (E.val - K.val)
				result.err = pre * (E.err + K.err)
				result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
				return err.ErrorSelect(stat_K, stat_E)
			}
		}
	} else if (x <= 0.0 && lambda < 1000.0) || (x < 0.1 && lambda < 17.0) || (x < 0.2 && lambda < 5.0) {
		return conicalP_xlt1_hyperg_A(1.0, lambda, x, result)
	} else if (x <= 0.2 && lambda < 17.0) || (x < 1.5 && lambda < 20.0) {
		var (
			arg    = math.Abs(x*x - 1.0)
			sgn    = gsl.Sign(1.0 - x)
			pre    = 0.5 * (lambda*lambda + 0.25) * sgn * math.Sqrt(arg)
			F      Result
			stat_F = Hyperg_2F1_conj_e(1.5, lambda, 2.0, (1.0-x)/2, &F)
		)

		result.val = pre * F.val
		result.err = math.Abs(pre) * F.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_F
	} else if 1.5 <= x && lambda < gsl.Max(x, 20.0) {
		var P Result
		var lm float64
		stat_P := ConicalP_large_x_e(1.0, lambda, x, &P, &lm)
		stat_e := Exp_mult_err_e(lm, 2.0*gsl.Float64Eps*math.Abs(lm), P.val, P.err, result)
		return err.ErrorSelect(stat_e, stat_P)
	} else {
		var V0, V1 float64
		if x < 1.0 {
			var (
				sqrt_1mx = math.Sqrt(1.0 - x)
				sqrt_1px = math.Sqrt(1.0 + x)
				th       = math.Acos(x)
				sth      = sqrt_1mx * sqrt_1px /* sin(th) */
				I0, I1   Result
			)
			stat_I0 := Bessel_I0_scaled_e(th*lambda, &I0)
			stat_I1 := Bessel_I1_scaled_e(th*lambda, &I1)
			stat_I := err.ErrorSelect(stat_I0, stat_I1)
			stat_V := conicalP_1_V(th, x/sth, lambda, -1.0, &V0, &V1)
			bessterm := V0*I0.val + V1*I1.val
			besserr := math.Abs(V0)*I0.err + math.Abs(V1)*I1.err + 2.0*gsl.Float64Eps*math.Abs(V0*I0.val) + 2.0*gsl.Float64Eps*math.Abs(V1*I1.val)
			arg1 := th * lambda
			sqts := math.Sqrt(th / sth)
			stat_e := Exp_mult_err_e(arg1, 2.0*gsl.Float64Eps*math.Abs(arg1), sqts*bessterm, sqts*besserr, result)
			result.err *= 1.0 / sqrt_1mx
			return err.ErrorSelect(stat_e, stat_V, stat_I)
		} else {
			var (
				sqrt_xm1 = math.Sqrt(x - 1.0)
				sqrt_xp1 = math.Sqrt(x + 1.0)
				sh       = sqrt_xm1 * sqrt_xp1 /* sinh(xi)      */
				xi       = math.Log(x + sh)    /* xi = acosh(x) */
				xi_lam   = xi * lambda
				J0, J1   Result
				stat_J0  = Bessel_J0_e(xi_lam, &J0)
				stat_J1  = Bessel_J1_e(xi_lam, &J1)
				stat_J   = err.ErrorSelect(stat_J0, stat_J1)
				stat_V   = conicalP_1_V(xi, x/sh, lambda, 1.0, &V0, &V1)
				bessterm = V0*J0.val + V1*J1.val
				besserr  = math.Abs(V0)*J0.err + math.Abs(V1)*J1.err + 512.0*2.0*gsl.Float64Eps*math.Abs(V0*J0.val) + 512.0*2.0*gsl.Float64Eps*math.Abs(V1*J1.val) + gsl.Float64Eps*math.Abs(xi_lam*V0*J1.val) + gsl.Float64Eps*math.Abs(xi_lam*V1*J0.val)
				pre      = math.Sqrt(xi / sh)
			)
			result.val = pre * bessterm
			result.err = pre * besserr * sqrt_xp1 / sqrt_xm1
			result.err += 4.0 * gsl.Float64Eps * math.Abs(result.val)
			return err.ErrorSelect(stat_V, stat_J)
		}
	}
}

/* P^{1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.8, 8.6.12]
 * checked OK [GJ] Fri May  8 12:24:36 MDT 1998
 */
func ConicalP_half_e(lambda, x float64, result *Result) err.GSLError {

	if x <= -1.0 {
		return DomainError(result)
	} else if x < 1.0 {
		err_amp := 1.0 + 1.0/(gsl.Float64Eps+math.Abs(1.0-math.Abs(x)))
		ac := math.Acos(x)
		den := math.Sqrt(math.Sqrt(1.0-x) * math.Sqrt(1.0+x))
		result.val = root_2OverPi_ / den * math.Cosh(ac*lambda)
		result.err = err_amp * 3.0 * gsl.Float64Eps * math.Abs(result.val)
		result.err *= math.Abs(ac*lambda) + 1.0
		return nil
	} else if x == 1.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		/* x > 1 */
		err_amp := 1.0 + 1.0/(gsl.Float64Eps+math.Abs(1.0-math.Abs(x)))
		sq_term := math.Sqrt(x-1.0) * math.Sqrt(x+1.0)
		ln_term := math.Log(x + sq_term)
		den := math.Sqrt(sq_term)
		carg_val := lambda * ln_term
		carg_err := 2.0 * gsl.Float64Eps * math.Abs(carg_val)
		var cos_result Result
		stat_cos := Cos_err_e(carg_val, carg_err, &cos_result)
		result.val = root_2OverPi_ / den * cos_result.val
		result.err = err_amp * root_2OverPi_ / den * cos_result.err
		result.err += 4.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_cos
	}
}

/* P^{-1/2}_{-1/2 + I lambda} (x)
 * [Abramowitz+Stegun 8.6.9, 8.6.14]
 * checked OK [GJ] Fri May  8 12:24:43 MDT 1998
 */
func ConicalP_mhalf_e(lambda, x float64, result *Result) err.GSLError {

	if x <= -1.0 {
		return DomainError(result)
	} else if x < 1.0 {
		ac := math.Acos(x)
		den := math.Sqrt(math.Sqrt(1.0-x) * math.Sqrt(1.0+x))
		arg := ac * lambda
		err_amp := 1.0 + 1.0/(gsl.Float64Eps+math.Abs(1.0-math.Abs(x)))
		if math.Abs(arg) < gsl.SqrtFloat64Eps {
			result.val = root_2OverPi_ / den * ac
			result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
			result.err *= err_amp
		} else {
			result.val = root_2OverPi_ / (den * lambda) * math.Sinh(arg)
			result.err = gsl.Float64Eps * (math.Abs(arg) + 1.0) * math.Abs(result.val)
			result.err *= err_amp
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		}
		return nil
	} else if x == 1.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		/* x > 1 */
		sq_term := math.Sqrt(x-1.0) * math.Sqrt(x+1.0)
		ln_term := math.Log(x + sq_term)
		den := math.Sqrt(sq_term)
		arg_val := lambda * ln_term
		arg_err := 2.0 * gsl.Float64Eps * math.Abs(arg_val)
		if arg_val < gsl.SqrtFloat64Eps {
			result.val = root_2OverPi_ / den * ln_term
			result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return nil
		} else {
			var sin_result Result
			stat_sin := Sin_err_e(arg_val, arg_err, &sin_result)
			result.val = root_2OverPi_ / (den * lambda) * sin_result.val
			result.err = root_2OverPi_ / math.Abs(den*lambda) * sin_result.err
			result.err += 3.0 * gsl.Float64Eps * math.Abs(result.val)
			return stat_sin
		}
	}
}

func ConicalP_sph_reg_e(l int, lambda, x float64, result *Result) err.GSLError {
	if x <= -1.0 || l < -1 {
		return DomainError(result)
	} else if l == -1 {
		return ConicalP_half_e(lambda, x, result)
	} else if l == 0 {
		return ConicalP_mhalf_e(lambda, x, result)
	} else if x == 1.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if x < 0.0 {
		c := 1.0 / math.Sqrt(1.0-x*x)
		var r_Pellm1, r_Pell Result
		stat_0 := ConicalP_half_e(lambda, x, &r_Pellm1) /* P^( 1/2) */
		stat_1 := ConicalP_mhalf_e(lambda, x, &r_Pell)  /* P^(-1/2) */
		stat_P := err.ErrorSelect(stat_0, stat_1)
		Pellm1 := r_Pellm1.val
		Pell := r_Pell.val
		var Pellp1 float64
		var ell int

		for ell = 0; ell < l; ell++ {
			fell := float64(ell)
			d := (fell+1.0)*(fell+1.0) + lambda*lambda
			Pellp1 = (Pellm1 - (2.0*fell+1.0)*c*x*Pell) / d
			Pellm1 = Pell
			Pell = Pellp1
		}

		result.val = Pell
		result.err = (0.5*float64(l) + 1.0) * gsl.Float64Eps * math.Abs(Pell)
		result.err += gsl.Float64Eps * float64(l) * math.Abs(result.val)
		return stat_P
	} else if x < 1.0 {
		xi := x / (math.Sqrt(1.0-x) * math.Sqrt(1.0+x))
		var rat, Phf Result
		stat_CF1 := conicalP_negmu_xlt1_CF1(0.5, l, lambda, x, &rat)
		stat_Phf := ConicalP_half_e(lambda, x, &Phf)
		Pellp1 := rat.val * gsl.SqrtMinFloat64
		Pell := gsl.SqrtMinFloat64
		var Pellm1 float64
		var ell int

		for ell = l; ell >= 0; ell-- {
			fell := float64(ell)
			d := (fell+1.0)*(fell+1.0) + lambda*lambda
			Pellm1 = (2.0*fell+1.0)*xi*Pell + d*Pellp1
			Pellp1 = Pell
			Pell = Pellm1
		}

		result.val = gsl.SqrtMinFloat64 * Phf.val / Pell
		result.err = gsl.SqrtMinFloat64 * Phf.err / math.Abs(Pell)
		result.err += math.Abs(rat.err/rat.val) * (float64(l) + 1.0) * math.Abs(result.val)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

		return err.ErrorSelect(stat_Phf, stat_CF1)
	} else if x == 1.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		/* x > 1.0 */
		xi := x / math.Sqrt((x-1.0)*(x+1.0))
		var rat Result
		stat_CF1 := conicalP_negmu_xgt1_CF1(0.5, l, lambda, x, &rat)
		var stat_P err.GSLError
		Pellp1 := rat.val * gsl.SqrtMinFloat64
		Pell := gsl.SqrtMinFloat64
		var Pellm1 float64
		var ell int

		for ell = l; ell >= 0; ell-- {
			fell := float64(ell)
			d := (fell+1.0)*(fell+1.0) + lambda*lambda
			Pellm1 = (2.0*fell+1.0)*xi*Pell - d*Pellp1
			Pellp1 = Pell
			Pell = Pellm1
		}

		if math.Abs(Pell) > math.Abs(Pellp1) {
			var Phf Result
			stat_P = ConicalP_half_e(lambda, x, &Phf)
			result.val = gsl.SqrtMinFloat64 * Phf.val / Pell
			result.err = 2.0 * gsl.SqrtMinFloat64 * Phf.err / math.Abs(Pell)
			result.err += 2.0 * math.Abs(rat.err/rat.val) * (float64(l) + 1.0) * math.Abs(result.val)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		} else {
			var Pmhf Result
			stat_P = ConicalP_mhalf_e(lambda, x, &Pmhf)
			result.val = gsl.SqrtMinFloat64 * Pmhf.val / Pellp1
			result.err = 2.0 * gsl.SqrtMinFloat64 * Pmhf.err / math.Abs(Pellp1)
			result.err += 2.0 * math.Abs(rat.err/rat.val) * (float64(l) + 1.0) * math.Abs(result.val)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		}

		return err.ErrorSelect(stat_P, stat_CF1)
	}
}

func ConicalP_cyl_reg_e(m int, lambda, x float64, result *Result) err.GSLError {

	if x <= -1.0 || m < -1 {
		return DomainError(result)
	} else if m == -1 {
		return ConicalP_1_e(lambda, x, result)
	} else if m == 0 {
		return ConicalP_0_e(lambda, x, result)
	} else if x == 1.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if x < 0.0 {
		c := 1.0 / math.Sqrt(1.0-x*x)
		var r_Pkm1, r_Pk Result
		stat_0 := ConicalP_1_e(lambda, x, &r_Pkm1) /* P^1 */
		stat_1 := ConicalP_0_e(lambda, x, &r_Pk)   /* P^0 */
		stat_P := err.ErrorSelect(stat_0, stat_1)
		Pkm1 := r_Pkm1.val
		Pk := r_Pk.val
		var Pkp1 float64
		var k int

		for k = 0; k < m; k++ {
			fk := float64(k)
			d := (fk+0.5)*(fk+0.5) + lambda*lambda
			Pkp1 = (Pkm1 - 2.0*fk*c*x*Pk) / d
			Pkm1 = Pk
			Pk = Pkp1
		}

		result.val = Pk
		result.err = (float64(m) + 2.0) * gsl.Float64Eps * math.Abs(Pk)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

		return stat_P
	} else if x < 1.0 {
		xi := x / (math.Sqrt(1.0-x) * math.Sqrt(1.0+x))
		var rat, P0 Result
		stat_CF1 := conicalP_negmu_xlt1_CF1(0.0, m, lambda, x, &rat)
		stat_P0 := ConicalP_0_e(lambda, x, &P0)
		Pkp1 := rat.val * gsl.SqrtMinFloat64
		Pk := gsl.SqrtMinFloat64
		var Pkm1 float64
		var k int

		for k = m; k > 0; k-- {
			fk := float64(k)
			d := (fk+0.5)*(fk+0.5) + lambda*lambda
			Pkm1 = 2.0*fk*xi*Pk + d*Pkp1
			Pkp1 = Pk
			Pk = Pkm1
		}

		result.val = gsl.SqrtMinFloat64 * P0.val / Pk
		result.err = 2.0 * gsl.SqrtMinFloat64 * P0.err / math.Abs(Pk)
		result.err += 2.0 * math.Abs(rat.err/rat.val) * (float64(m) + 1.0) * math.Abs(result.val)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

		return err.ErrorSelect(stat_P0, stat_CF1)
	} else if x == 1.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		/* x > 1.0 */

		xi := x / math.Sqrt((x-1.0)*(x+1.0))
		var rat Result
		stat_CF1 := conicalP_negmu_xgt1_CF1(0.0, m, lambda, x, &rat)
		var stat_P err.GSLError
		Pkp1 := rat.val * gsl.SqrtMinFloat64
		Pk := gsl.SqrtMinFloat64
		var Pkm1 float64
		var k int

		for k = m; k > -1; k-- {
			fk := float64(k)
			d := (fk+0.5)*(fk+0.5) + lambda*lambda
			Pkm1 = 2.0*fk*xi*Pk - d*Pkp1
			Pkp1 = Pk
			Pk = Pkm1
		}

		if math.Abs(Pk) > math.Abs(Pkp1) {
			var P1 Result
			stat_P = ConicalP_1_e(lambda, x, &P1)
			result.val = gsl.SqrtMinFloat64 * P1.val / Pk
			result.err = 2.0 * gsl.SqrtMinFloat64 * P1.err / math.Abs(Pk)
			result.err += 2.0 * math.Abs(rat.err/rat.val) * (float64(m) + 2.0) * math.Abs(result.val)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		} else {
			var P0 Result
			stat_P = ConicalP_0_e(lambda, x, &P0)
			result.val = gsl.SqrtMinFloat64 * P0.val / Pkp1
			result.err = 2.0 * gsl.SqrtMinFloat64 * P0.err / math.Abs(Pkp1)
			result.err += 2.0 * math.Abs(rat.err/rat.val) * (float64(m) + 2.0) * math.Abs(result.val)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		}

		return err.ErrorSelect(stat_P, stat_CF1)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func ConicalP_0(lambda, x float64) float64 {
	result := new(Result)
	status := ConicalP_0_e(lambda, x, result)
	return EvalResult(result, status)
}

func ConicalP_1(lambda, x float64) float64 {
	result := new(Result)
	status := ConicalP_1_e(lambda, x, result)
	return EvalResult(result, status)
}

func ConicalP_half(lambda, x float64) float64 {
	result := new(Result)
	status := ConicalP_half_e(lambda, x, result)
	return EvalResult(result, status)
}

func ConicalP_mhalf(lambda, x float64) float64 {
	result := new(Result)
	status := ConicalP_mhalf_e(lambda, x, result)
	return EvalResult(result, status)
}

func ConicalP_sph_reg(l int, lambda, x float64) float64 {
	result := new(Result)
	status := ConicalP_sph_reg_e(l, lambda, x, result)
	return EvalResult(result, status)
}
