/* specfunc/legendre_H3d.c
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

/* See [Abbott+Schaefer, Ap.J. 308, 546 (1986)] for
 * enough details to follow what is happening here.
 */

/* Logarithm of normalization factor, Log[N(ell,lambda)].
 * N(ell,lambda) = Product[ lambda^2 + n^2, {n,0,ell} ]
 *               = |Gamma(ell + 1 + I lambda)|^2  lambda sinh(Pi lambda) / Pi
 * Assumes ell >= 0.
 */
func legendre_H3d_lnnorm(ell int, lambda float64, result *float64) err.GSLError {
	abs_lam := math.Abs(lambda)
	fell := float64(ell)
	if abs_lam == 0.0 {
		*result = 0.0
		return err.ERROR("error", err.EDOM)
	} else if lambda > (fell+1.0)/gsl.Root3Float64Eps {
		/*
		 * There is a cancellation between the sinh(Pi lambda)
		 * term and the math.Log(gamma(ell + 1 + i lambda) in the
		 * result below, so we show some care and save some digits.
		 * Note that the above guarantees that lambda is large,
		 * since ell >= 0. We use Stirling and a simple expansion
		 * of sinh.
		 */
		rat := (fell + 1.0) / lambda
		ln_lam2ell2 := 2.0*math.Log(lambda) + math.Log(1.0+rat*rat)
		lg_corrected := -2.0*(fell+1.0) + gsl.LnPi + (fell+0.5)*ln_lam2ell2 + 1.0/(288.0*lambda*lambda)
		angle_terms := lambda * 2.0 * rat * (1.0 - rat*rat/3.0)
		*result = math.Log(abs_lam) + lg_corrected + angle_terms - gsl.LnPi
		return nil
	} else {
		var lg_r, lg_theta, ln_sinh Result
		Lngamma_complex_e(fell+1.0, lambda, &lg_r, &lg_theta)
		Lnsinh_e(gsl.Pi*abs_lam, &ln_sinh)
		*result = math.Log(abs_lam) + ln_sinh.val + 2.0*lg_r.val - gsl.LnPi
		return nil
	}
}

/* Calculate series for small eta*lambda.
 * Assumes eta > 0, lambda != 0.
 *
 * This is just the defining hypergeometric for the Legendre function.
 *
 * P^{mu}_{-1/2 + I lam}(z) = 1/Gamma(l+3/2) ((z+1)/(z-1)^(mu/2)
 *                            2F1(1/2 - I lam, 1/2 + I lam; l+3/2; (1-z)/2)
 * We use
 *       z = cosh(eta)
 * (z-1)/2 = sinh^2(eta/2)
 *
 * And recall
 * H3d = sqrt(Pi Norm /(2 lam^2 sinh(eta))) P^{-l-1/2}_{-1/2 + I lam}(cosh(eta))
 */
func legendre_H3d_series(ell int, lambda, eta float64, result *Result) err.GSLError {
	var (
		nmax                           = 5000
		shheta                         = math.Sinh(0.5 * eta)
		ln_zp1                         = gsl.Ln2 + math.Log(1.0+shheta*shheta)
		ln_zm1                         = gsl.Ln2 + 2.0*math.Log(shheta)
		zeta                           = -shheta * shheta
		lg_lp32                        Result
		term                           = 1.0
		sum                            = 1.0
		sum_err                        = 0.0
		lnsheta                        Result
		lnN                            float64
		lnpre_val, lnpre_err, lnprepow float64
		stat_e                         err.GSLError
		n                              int
		fell                           = float64(ell)
	)

	Lngamma_e(fell+3.0/2.0, &lg_lp32)
	Lnsinh_e(eta, &lnsheta)
	legendre_H3d_lnnorm(ell, lambda, &lnN)
	lnprepow = 0.5 * (fell + 0.5) * (ln_zm1 - ln_zp1)
	lnpre_val = lnprepow + 0.5*(lnN+gsl.LnPi-gsl.Ln2-lnsheta.val) - lg_lp32.val - math.Log(math.Abs(lambda))
	lnpre_err = lnsheta.err + lg_lp32.err + gsl.Float64Eps*math.Abs(lnpre_val)
	lnpre_err += 2.0 * gsl.Float64Eps * (math.Abs(lnN) + gsl.LnPi + gsl.Ln2)
	lnpre_err += 2.0 * gsl.Float64Eps * (0.5 * (fell + 0.5) * (math.Abs(ln_zm1) + math.Abs(ln_zp1)))
	for n = 1; n < nmax; n++ {
		fn := float64(n)
		aR := fn - 0.5
		term *= (aR*aR + lambda*lambda) * zeta / (fell + fn + 0.5) / fn
		sum += term
		sum_err += 2.0 * gsl.Float64Eps * math.Abs(term)
		if math.Abs(term/sum) < 2.0*gsl.Float64Eps {
			break
		}
	}

	stat_e = Exp_mult_err_e(lnpre_val, lnpre_err, sum, math.Abs(term)+sum_err, result)
	var s err.GSLError

	if n == nmax {
		s = err.MaxIteration()
	}

	return err.ErrorSelect(stat_e, s)
}

/* Evaluate legendre_H3d(ell+1)/legendre_H3d(ell)
 * by continued fraction. Use the Gautschi (Euler)
 * equivalent series.
 */
/* FIXME: Maybe we have to worry about this. The a_k are
 * not positive and there can be a blow-up. It happened
 * for J_nu once or twice. Then we should probably use
 * the method above.
 */
func legendre_H3d_CF1_ser(ell int, lambda, coth_eta float64, result *Result) err.GSLError {
	var (
		fell    = float64(ell)
		pre     = math.Hypot(lambda, fell+1.0) / ((2.0*fell + 3) * coth_eta)
		maxk    = 20000
		tk      = 1.0
		sum     = 1.0
		rhok    = 0.0
		sum_err = 0.0
		k       int
	)

	for k = 1; k < maxk; k++ {
		fk := float64(k)
		tlk := (2.0*fell + 1.0 + 2.0*fk)
		l1k := (fell + 1.0 + fk)
		ak := -(lambda*lambda + l1k*l1k) / (tlk * (tlk + 2.0) * coth_eta * coth_eta)
		rhok = -ak * (1.0 + rhok) / (1.0 + ak*(1.0+rhok))
		tk *= rhok
		sum += tk
		sum_err += 2.0 * gsl.Float64Eps * fk * math.Abs(tk)
		if math.Abs(tk/sum) < gsl.Float64Eps {
			break
		}
	}

	result.val = pre * sum
	result.err = math.Abs(pre * tk)
	result.err += math.Abs(pre * sum_err)
	result.err += 4.0 * gsl.Float64Eps * math.Abs(result.val)

	if k >= maxk {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Legendre_H3d_0_e(lambda, eta float64, result *Result) err.GSLError {
	if eta < 0.0 {
		return DomainError(result)
	} else if eta == 0.0 || lambda == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else {
		lam_eta := lambda * eta
		var s Result
		Sin_err_e(lam_eta, 2.0*gsl.Float64Eps*math.Abs(lam_eta), &s)
		if eta > -0.5*gsl.LnFloat64Eps {
			f := 2.0 / lambda * math.Exp(-eta)
			result.val = f * s.val
			result.err = math.Abs(f*s.val) * (math.Abs(eta) + 1.0) * gsl.Float64Eps
			result.err += math.Abs(f) * s.err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		} else {
			f := 1.0 / (lambda * math.Sinh(eta))
			result.val = f * s.val
			result.err = math.Abs(f*s.val) * (math.Abs(eta) + 1.0) * gsl.Float64Eps
			result.err += math.Abs(f) * s.err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		}
		return nil
	}
}

func Legendre_H3d_1_e(lambda, eta float64, result *Result) err.GSLError {
	var (
		xi    = math.Abs(eta * lambda)
		lsq   = lambda * lambda
		lsqp1 = lsq + 1.0
	)

	if eta < 0.0 {
		return DomainError(result)
	} else if eta == 0.0 || lambda == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if xi < gsl.Root5Float64Eps && eta < gsl.Root5Float64Eps {
		etasq := eta * eta
		xisq := xi * xi
		term1 := (etasq + xisq) / 3.0
		term2 := -(2.0*etasq*etasq + 5.0*etasq*xisq + 3.0*xisq*xisq) / 90.0
		sinh_term := 1.0 - eta*eta/6.0*(1.0-7.0/60.0*eta*eta)
		pre := sinh_term / math.Sqrt(lsqp1) / eta
		result.val = pre * (term1 + term2)
		result.err = pre * gsl.Float64Eps * (math.Abs(term1) + math.Abs(term2))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		var (
			sin_term     float64 /*  Sin(xi)/xi     */
			cos_term     float64 /*  Cos(xi)        */
			coth_term    float64 /*  eta/Tanh(eta)  */
			sinh_term    float64 /*  eta/Sinh(eta)  */
			sin_term_err float64
			cos_term_err float64
			t1           float64
			pre_val      float64
			pre_err      float64
			term1        float64
			term2        float64
		)
		if xi < gsl.Root5Float64Eps {
			sin_term = 1.0 - xi*xi/6.0*(1.0-xi*xi/20.0)
			cos_term = 1.0 - 0.5*xi*xi*(1.0-xi*xi/12.0)
			sin_term_err = gsl.Float64Eps
			cos_term_err = gsl.Float64Eps
		} else {
			var sin_xi_result, cos_xi_result Result
			Sin_e(xi, &sin_xi_result)
			Cos_e(xi, &cos_xi_result)
			sin_term = sin_xi_result.val / xi
			cos_term = cos_xi_result.val
			sin_term_err = sin_xi_result.err / math.Abs(xi)
			cos_term_err = cos_xi_result.err
		}
		if eta < gsl.Root5Float64Eps {
			coth_term = 1.0 + eta*eta/3.0*(1.0-eta*eta/15.0)
			sinh_term = 1.0 - eta*eta/6.0*(1.0-7.0/60.0*eta*eta)
		} else {
			coth_term = eta / math.Tanh(eta)
			sinh_term = eta / math.Sinh(eta)
		}
		t1 = math.Sqrt(lsqp1) * eta
		pre_val = sinh_term / t1
		pre_err = 2.0 * gsl.Float64Eps * math.Abs(pre_val)
		term1 = sin_term * coth_term
		term2 = cos_term
		result.val = pre_val * (term1 - term2)
		result.err = pre_err * math.Abs(term1-term2)
		result.err += pre_val * (sin_term_err*coth_term + cos_term_err)
		result.err += pre_val * math.Abs(term1-term2) * (math.Abs(eta) + 1.0) * gsl.Float64Eps
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Legendre_H3d_e(ell int, lambda, eta float64, result *Result) err.GSLError {
	var (
		fell     = float64(ell)
		abs_lam  = math.Abs(lambda)
		lsq      = abs_lam * abs_lam
		xi       = abs_lam * eta
		cosh_eta = math.Cosh(eta)
	)

	if eta < 0.0 {
		return DomainError(result)
	} else if eta > gsl.LnMaxFloat64 {
		/* cosh(eta) is too big. */
		return OverflowError(result)
	} else if ell == 0 {
		return Legendre_H3d_0_e(lambda, eta, result)
	} else if ell == 1 {
		return Legendre_H3d_1_e(lambda, eta, result)
	} else if eta == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if xi < 1.0 {

		return legendre_H3d_series(ell, lambda, eta, result)
	} else if (fell*fell+lsq)/math.Sqrt(1.0+lsq)/(cosh_eta*cosh_eta) < 5.0*gsl.Root3Float64Eps {
		/* Large argument.
		 */
		var P Result
		var lm float64
		stat_P := ConicalP_large_x_e(-fell-0.5, lambda, cosh_eta, &P, &lm)
		if P.val == 0.0 {
			result.val = 0.0
			result.err = 0.0
			return stat_P
		} else {
			var lnN float64
			var lnsh Result
			var ln_abslam, lnpre_val, lnpre_err float64
			var stat_e err.GSLError
			Lnsinh_e(eta, &lnsh)
			legendre_H3d_lnnorm(ell, lambda, &lnN)
			ln_abslam = math.Log(abs_lam)
			lnpre_val = 0.5*(gsl.LnPi+lnN-gsl.Ln2-lnsh.val) - ln_abslam
			lnpre_err = lnsh.err
			lnpre_err += 2.0 * gsl.Float64Eps * (0.5*(gsl.LnPi+gsl.Ln2+math.Abs(lnN)) + math.Abs(ln_abslam))
			lnpre_err += 2.0 * gsl.Float64Eps * math.Abs(lnpre_val)
			stat_e = Exp_mult_err_e(lnpre_val+lm, lnpre_err, P.val, P.err, result)
			return err.ErrorSelect(stat_e, stat_P)
		}
	} else if abs_lam > 1000.0*fell*fell {
		/* Large degree.
		 */
		var P Result
		var lm float64
		stat_P := ConicalP_xgt1_neg_mu_largetau_e(fell+0.5, lambda, cosh_eta, eta, &P, &lm)
		if P.val == 0.0 {
			result.val = 0.0
			result.err = 0.0
			return stat_P
		} else {
			var (
				lnN                  float64
				lnsh                 Result
				ln_abslam            float64
				lnpre_val, lnpre_err float64
				stat_e               err.GSLError
			)
			Lnsinh_e(eta, &lnsh)
			legendre_H3d_lnnorm(ell, lambda, &lnN)
			ln_abslam = math.Log(abs_lam)
			lnpre_val = 0.5*(gsl.LnPi+lnN-gsl.Ln2-lnsh.val) - ln_abslam
			lnpre_err = lnsh.err
			lnpre_err += gsl.Float64Eps * (0.5*(gsl.LnPi+gsl.Ln2+math.Abs(lnN)) + math.Abs(ln_abslam))
			lnpre_err += 2.0 * gsl.Float64Eps * math.Abs(lnpre_val)
			stat_e = Exp_mult_err_e(lnpre_val+lm, lnpre_err, P.val, P.err, result)
			return err.ErrorSelect(stat_e, stat_P)
		}
	} else {
		/* Backward recurrence.
		 */
		coth_eta := 1.0 / math.Tanh(eta)
		coth_err_mult := math.Abs(eta) + 1.0
		var rH Result
		stat_CF1 := legendre_H3d_CF1_ser(ell, lambda, coth_eta, &rH)
		var Hlm1 float64
		Hl := gsl.SqrtMinFloat64
		Hlp1 := rH.val * Hl
		var lp int
		for lp = ell; lp > 0; lp-- {
			flp := float64(lp)
			root_term_0 := math.Hypot(lambda, flp)
			root_term_1 := math.Hypot(lambda, flp+1.0)
			Hlm1 = ((2.0*flp+1.0)*coth_eta*Hl - root_term_1*Hlp1) / root_term_0
			Hlp1 = Hl
			Hl = Hlm1
		}

		if math.Abs(Hl) > math.Abs(Hlp1) {
			var H0 Result
			stat_H0 := Legendre_H3d_0_e(lambda, eta, &H0)
			result.val = gsl.SqrtMinFloat64 / Hl * H0.val
			result.err = gsl.SqrtMinFloat64 / math.Abs(Hl) * H0.err
			result.err += math.Abs(rH.err/rH.val) * (fell + 1.0) * coth_err_mult * math.Abs(result.val)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return err.ErrorSelect(stat_H0, stat_CF1)
		} else {
			var H1 Result
			stat_H1 := Legendre_H3d_1_e(lambda, eta, &H1)
			result.val = gsl.SqrtMinFloat64 / Hlp1 * H1.val
			result.err = gsl.SqrtMinFloat64 / math.Abs(Hlp1) * H1.err
			result.err += math.Abs(rH.err/rH.val) * (fell + 1.0) * coth_err_mult * math.Abs(result.val)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return err.ErrorSelect(stat_H1, stat_CF1)
		}
	}
}

func Legendre_H3d_array(lmax int, lambda, eta float64, result_array []float64) err.GSLError {
	if eta < 0.0 || lmax < 0 {
		var ell int
		for ell = 0; ell <= lmax; ell++ {
			result_array[ell] = 0.0
		}
		return err.ERROR("domain error", err.EDOM)
	} else if eta > gsl.LnMaxFloat64 {
		/* cosh(eta) is too big. */
		var ell int
		for ell = 0; ell <= lmax; ell++ {
			result_array[ell] = 0.0
		}
		return err.ERROR("overflow", err.EOVRFLW)
	} else if lmax == 0 {
		var H0 Result
		stat := Legendre_H3d_e(0, lambda, eta, &H0)
		result_array[0] = H0.val
		return stat
	} else {
		/* Not the most efficient method. But what the hell... it's simple.
		 */
		var r_Hlp1, r_Hl Result
		stat_lmax := Legendre_H3d_e(lmax, lambda, eta, &r_Hlp1)
		stat_lmaxm1 := Legendre_H3d_e(lmax-1, lambda, eta, &r_Hl)
		stat_max := err.ErrorSelect(stat_lmax, stat_lmaxm1)

		coth_eta := 1.0 / math.Tanh(eta)
		var stat_recursion err.GSLError
		Hlp1 := r_Hlp1.val
		Hl := r_Hl.val
		var Hlm1 float64
		var ell int

		result_array[lmax] = Hlp1
		result_array[lmax-1] = Hl

		for ell = lmax - 1; ell > 0; ell-- {
			fell := float64(ell)
			root_term_0 := math.Hypot(lambda, fell)
			root_term_1 := math.Hypot(lambda, fell+1.0)
			Hlm1 = ((2.0*fell+1.0)*coth_eta*Hl - root_term_1*Hlp1) / root_term_0
			result_array[ell-1] = Hlm1
			if !(Hlm1 < gsl.MaxFloat64) {
				stat_recursion = err.Overflow()
			}
			Hlp1 = Hl
			Hl = Hlm1
		}

		return err.ErrorSelect(stat_recursion, stat_max)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Legendre_H3d_0(lambda, eta float64) float64 {
	result := new(Result)
	status := Legendre_H3d_0_e(lambda, eta, result)
	return EvalResult(result, status)
}

func Legendre_H3d_1(lambda, eta float64) float64 {
	result := new(Result)
	status := Legendre_H3d_1_e(lambda, eta, result)
	return EvalResult(result, status)
}

func Legendre_H3d(l int, lambda, eta float64) float64 {
	result := new(Result)
	status := Legendre_H3d_e(l, lambda, eta, result)
	return EvalResult(result, status)
}
