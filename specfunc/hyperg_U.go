/* specfunc/hyperg_U.c
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
	"github.com/jtejido/ggsl/sys"
	"math"
)

var (
	INT_THRESHOLD = 1000.0 * gsl.Float64Eps
)

func seriesEvalOk(a, b, x float64) bool {
	return ((math.Abs(a) < 5 && b < 5 && x < 2.0) || (math.Abs(a) < 10 && b < 10 && x < 1.0))
}

func AsympEvalOk(a, b, x float64) bool {
	return (math.Max(math.Abs(a), 1.0)*math.Max(math.Abs(1.0+a-b), 1.0) < 0.99*math.Abs(x))
}

/* Log[U(a,2a,x)]
 * [Abramowitz+stegun, 13.6.21]
 * Assumes x > 0, a > 1/2.
 */
func hyperg_lnU_beq2a(a, x float64, result *Result) err.GSLError {
	lx := math.Log(x)
	nu := a - 0.5
	lnpre := 0.5*(x-gsl.LnPi) - nu*lx
	lnK := new(Result)
	Bessel_lnKnu_e(nu, 0.5*x, lnK)
	result.val = lnpre + lnK.val
	result.err = 2.0 * gsl.Float64Eps * (math.Abs(0.5*x) + 0.5*gsl.LnPi + math.Abs(nu*lx))
	result.err += lnK.err
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return nil
}

/* Evaluate u_{N+1}/u_N by Steed's continued fraction method.
 *
 * u_N := Gamma[a+N]/Gamma[a] U(a + N, b, x)
 *
 * u_{N+1}/u_N = (a+N) U(a+N+1,b,x)/U(a+N,b,x)
 */
func hyperg_U_CF1(a, b float64, N int, x float64, result *float64, count *int) err.GSLError {
	var (
		RECUR_BIG = gsl.SqrtMaxFloat64
		maxiter   = 20000
		n         = 1
		Anm2      = 1.0
		Bnm2      = 0.0
		Anm1      = 0.0
		Bnm1      = 1.0
		a1        = -(a + float64(N))
		b1        = (b - 2.0*a - x - 2.0*(float64(N)+1))
		An        = b1*Anm1 + a1*Anm2
		Bn        = b1*Bnm1 + a1*Bnm2
		an, bn    float64
		fn        = An / Bn
	)

	for n < maxiter {
		var old_fn, del float64
		n++
		fN := float64(N)
		ffn := float64(n)
		Anm2 = Anm1
		Bnm2 = Bnm1
		Anm1 = An
		Bnm1 = Bn
		an = -(a + fN + ffn - b) * (a + fN + ffn - 1.0)
		bn = (b - 2.0*a - x - 2.0*(fN+ffn))
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

		if math.Abs(del-1.0) < 10.0*gsl.Float64Eps {
			break
		}
	}

	*result = fn
	*count = n

	if n == maxiter {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

/* Large x asymptotic for  x^a U(a,b,x)
 * Based on SLATEC D9CHU() [W. Fullerton]
 *
 * Uses a rational approximation due to Luke.
 * See [Luke, Algorithms for the Computation of Special Functions, p. 252]
 *     [Luke, Utilitas Math. (1977)]
 *
 * z^a U(a,b,z) ~ 2F0(a,1+a-b,-1/z)
 *
 * This assumes that a is not a negative integer and
 * that 1+a-b is not a negative integer. If one of them
 * is, then the 2F0 actually terminates, the above
 * relation is an equality, and the sum should be
 * evaluated directly [see below].
 */
func d9chu(a, b, x float64, result *Result) err.GSLError {
	var (
		EPS     = 8.0 * gsl.Float64Eps /* EPS = 4.0D0*D1MACH(4)   */
		maxiter = 500
		aa, bb  [4]float64
		i       int

		bp  = 1.0 + a - b
		ab  = a * bp
		ct2 = 2.0 * (x - ab)
		sab = a + bp

		ct3  = sab + 1.0 + ab
		anbn = ct3 + sab + 3.0
		ct1  = 1.0 + 2.0*x/anbn
	)
	bb[0] = 1.0
	aa[0] = 1.0

	bb[1] = 1.0 + 2.0*x/ct3
	aa[1] = 1.0 + ct2/ct3

	bb[2] = 1.0 + 6.0*ct1*x/ct3
	aa[2] = 1.0 + 6.0*ab/anbn + 3.0*ct1*ct2/ct3

	for i = 4; i < maxiter; i++ {
		var (
			j          int
			c2         float64
			d1z        float64
			g1, g2, g3 float64
			x2i1       = 2*float64(i) - 3
		)
		ct1 = x2i1 / (x2i1 - 2.0)
		anbn += x2i1 + sab
		ct2 = (x2i1 - 1.0) / anbn
		c2 = x2i1*ct2 - 1.0
		d1z = 2.0 * x2i1 * x / anbn

		ct3 = sab * ct2
		g1 = d1z + ct1*(c2+ct3)
		g2 = d1z - c2
		g3 = ct1 * (1.0 - ct3 - 2.0*ct2)

		bb[3] = g1*bb[2] + g2*bb[1] + g3*bb[0]
		aa[3] = g1*aa[2] + g2*aa[1] + g3*aa[0]

		if math.Abs(aa[3]*bb[0]-aa[0]*bb[3]) < EPS*math.Abs(bb[3]*bb[0]) {
			break
		}

		for j = 0; j < 3; j++ {
			aa[j] = aa[j+1]
			bb[j] = bb[j+1]
		}
	}

	result.val = aa[3] / bb[3]
	result.err = 8.0 * gsl.Float64Eps * math.Abs(result.val)

	if i == maxiter {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil

}

/* Evaluate asymptotic for z^a U(a,b,z) ~ 2F0(a,1+a-b,-1/z)
 * We check for termination of the 2F0 as a special case.
 * Assumes x > 0.
 * Also assumes a,b are not too large compared to x.
 */
func hyperg_zaU_asymp(a, b, x float64, result *Result) err.GSLError {
	var (
		ap         = a
		bp         = 1.0 + a - b
		rintap     = math.Floor(ap + 0.5)
		rintbp     = math.Floor(bp + 0.5)
		ap_neg_int = (ap < 0.0 && math.Abs(ap-rintap) < INT_THRESHOLD)
		bp_neg_int = (bp < 0.0 && math.Abs(bp-rintbp) < INT_THRESHOLD)
	)

	if ap_neg_int || bp_neg_int {
		/* Evaluate 2F0 polynomial.
		 */
		var (
			mxi     = -1.0 / x
			nmax    = -int(math.Min(ap, bp) - 0.1)
			tn      = 1.0
			sum     = 1.0
			n       = 1.0
			sum_err = 0.0
		)
		for n <= float64(nmax) {
			apn := (ap + n - 1.0)
			bpn := (bp + n - 1.0)
			tn *= ((apn / n) * mxi) * bpn
			sum += tn
			sum_err += 2.0 * gsl.Float64Eps * math.Abs(tn)
			n += 1.0
		}
		result.val = sum
		result.err = sum_err
		result.err += 2.0 * gsl.Float64Eps * (math.Abs(float64(nmax)) + 1.0) * math.Abs(sum)
		return nil
	} else {
		return d9chu(a, b, x, result)
	}
}

/* Evaluate finite sum which appears below.
 */
func hyperg_U_finite_sum(N int, a, b, x, xeps float64, result *Result) err.GSLError {
	var (
		i       int
		sum_val float64
		sum_err float64
	)

	if N <= 0 {
		t_val := 1.0
		t_err := 0.0
		poch := new(Result)
		var stat_poch err.GSLError

		sum_val = 1.0
		sum_err = 0.0
		for i = 1; i <= -N; i++ {
			xi1 := float64(i) - 1
			mult := (a + xi1) * x / ((b + xi1) * (xi1 + 1.0))
			t_val *= mult
			t_err += math.Abs(mult)*t_err + math.Abs(t_val)*8.0*2.0*gsl.Float64Eps
			sum_val += t_val
			sum_err += t_err
		}

		stat_poch = Poch_e(1.0+a-b, -a, poch)

		result.val = sum_val * poch.val
		result.err = math.Abs(sum_val)*poch.err + sum_err*math.Abs(poch.val)
		result.err += math.Abs(poch.val) * (math.Abs(float64(N)) + 2.0) * gsl.Float64Eps * math.Abs(sum_val)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		result.err *= 2.0 /* FIXME: fudge factor... why is the error estimate too small? */
		return stat_poch
	} else {
		M := N - 2
		if M < 0 {
			result.val = 0.0
			result.err = 0.0
			return nil
		} else {
			gbm1, gamr := new(Result), new(Result)
			var stat_gbm1 err.GSLError
			t_val := 1.0
			t_err := 0.0

			sum_val = 1.0
			sum_err = 0.0
			for i = 1; i <= M; i++ {
				fi := float64(i)
				mult := (a - b + fi) * x / ((1.0 - b + fi) * fi)
				t_val *= mult
				t_err += t_err*math.Abs(mult) + math.Abs(t_val)*8.0*2.0*gsl.Float64Eps
				sum_val += t_val
				sum_err += t_err
			}

			stat_gbm1 = Gamma_e(b-1.0, gbm1)
			Gammainv_e(a, gamr) // stat_gamr

			if stat_gbm1 == nil {
				powx1N := new(Result)
				stat_p := Pow_int_e(x, 1-N, powx1N)
				pe_val := powx1N.val * xeps
				pe_err := powx1N.err*math.Abs(xeps) + 2.0*gsl.Float64Eps*math.Abs(pe_val)
				coeff_val := gbm1.val * gamr.val * pe_val
				coeff_err := gbm1.err*math.Abs(gamr.val*pe_val) + gamr.err*math.Abs(gbm1.val*pe_val) + math.Abs(gbm1.val*gamr.val)*pe_err + 2.0*gsl.Float64Eps*math.Abs(coeff_val)

				result.val = sum_val * coeff_val
				result.err = math.Abs(sum_val)*coeff_err + sum_err*math.Abs(coeff_val)
				result.err += 2.0 * gsl.Float64Eps * (float64(M) + 2.0) * math.Abs(result.val)
				result.err *= 2.0 /* FIXME: fudge factor... why is the error estimate too small? */
				return stat_p
			} else {
				result.val = 0.0
				result.err = 0.0
				return stat_gbm1
			}
		}
	}
}

/* Evaluate infinite sum which appears below.
 */
func hyperg_U_infinite_sum_stable(N int, a, bint, b, beps, x, xeps float64, sum Result, result *Result) err.GSLError {

	EPS := 2.0 * gsl.Float64Eps /* EPS = D1MACH(3) */

	var istrt int
	if N < 1 {
		istrt = 1 - N
	}
	xi := float64(istrt)

	gamr, powx := new(Result), new(Result)
	stat_gamr := Gammainv_e(1.0+a-b, gamr)
	stat_powx := Pow_int_e(x, istrt, powx)
	sarg := beps * gsl.Pi
	sfact := 1.0
	if sarg != 0.0 {
		sfact = sarg / math.Sin(sarg)
	}
	factor_val := sfact
	if gsl.IsOdd(N) {
		factor_val *= -1
	}
	factor_val *= gamr.val * powx.val
	factor_err := math.Abs(gamr.val)*powx.err + math.Abs(powx.val)*gamr.err + 2.0*gsl.Float64Eps*math.Abs(factor_val)

	pochai, gamri1, gamrni := new(Result), new(Result), new(Result)
	stat_pochai := Poch_e(a, xi, pochai)
	stat_gamri1 := Gammainv_e(xi+1.0, gamri1)
	stat_gamrni := Gammainv_e(bint+xi, gamrni)
	stat_gam123 := err.ErrorSelect(stat_gamr, stat_gamri1, stat_gamrni)
	stat_gamall := err.ErrorSelect(stat_gam123, stat_pochai, stat_powx)

	pochaxibeps, gamrxi1beps := new(Result), new(Result)
	stat_pochaxibeps := Poch_e(a, xi-beps, pochaxibeps)
	stat_gamrxi1beps := Gammainv_e(xi+1.0-beps, gamrxi1beps)

	stat_all := err.ErrorSelect(stat_gamall, stat_pochaxibeps, stat_gamrxi1beps)

	b0_val := factor_val * pochaxibeps.val * gamrni.val * gamrxi1beps.val
	b0_err := math.Abs(factor_val*pochaxibeps.val*gamrni.val)*gamrxi1beps.err + math.Abs(factor_val*pochaxibeps.val*gamrxi1beps.val)*gamrni.err + math.Abs(factor_val*gamrni.val*gamrxi1beps.val)*pochaxibeps.err + math.Abs(pochaxibeps.val*gamrni.val*gamrxi1beps.val)*factor_err + 2.0*gsl.Float64Eps*math.Abs(b0_val)

	/*
	   C  X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE
	   C  STRAIGHTFORWARD FORMULATION IS STABLE.
	*/
	var (
		i                                int
		dchu_val, dchu_err, t_val, t_err float64
		dgamrbxi                         Result
	)

	stat_dgamrbxi := Gammainv_e(b+xi, &dgamrbxi)
	a0_val := factor_val * pochai.val * dgamrbxi.val * gamri1.val / beps
	a0_err := math.Abs(factor_val*pochai.val*dgamrbxi.val/beps)*gamri1.err + math.Abs(factor_val*pochai.val*gamri1.val/beps)*dgamrbxi.err + math.Abs(factor_val*dgamrbxi.val*gamri1.val/beps)*pochai.err + math.Abs(pochai.val*dgamrbxi.val*gamri1.val/beps)*factor_err + 2.0*gsl.Float64Eps*math.Abs(a0_val)
	stat_all = err.ErrorSelect(stat_all, stat_dgamrbxi)

	b0_val = xeps * b0_val / beps
	b0_err = math.Abs(xeps/beps)*b0_err + 4.0*gsl.Float64Eps*math.Abs(b0_val)
	dchu_val = sum.val + a0_val - b0_val
	dchu_err = sum.err + a0_err + b0_err + 2.0*gsl.Float64Eps*(math.Abs(sum.val)+math.Abs(a0_val)+math.Abs(b0_val))

	for i = 1; i < 2000; i++ {
		xi := float64(istrt + i)
		xi1 := float64(istrt + i - 1)
		a0_multiplier := (a + xi1) * x / ((b + xi1) * xi)
		b0_multiplier := (a + xi1 - beps) * x / ((bint + xi1) * (xi - beps))
		a0_val *= a0_multiplier
		a0_err += math.Abs(a0_multiplier) * a0_err
		b0_val *= b0_multiplier
		b0_err += math.Abs(b0_multiplier) * b0_err
		t_val = a0_val - b0_val
		t_err = a0_err + b0_err
		dchu_val += t_val
		dchu_err += t_err
		if math.Abs(t_val) < EPS*math.Abs(dchu_val) {
			break
		}
	}

	result.val = dchu_val
	result.err = 2.0 * dchu_err
	result.err += 2.0 * math.Abs(t_val)
	result.err += 4.0 * gsl.Float64Eps * (float64(i) + 2.0) * math.Abs(dchu_val)
	result.err *= 2.0 /* FIXME: fudge factor */

	if i >= 2000 {
		return err.ERROR("error", err.EMAXITER)
	}

	return stat_all

}

func hyperg_U_infinite_sum_simple(N int, a, bint, b, beps, x, xeps float64, sum Result, result *Result) err.GSLError {
	EPS := 2.0 * gsl.Float64Eps /* EPS = D1MACH(3) */
	var istrt int
	if N < 1 {
		istrt = 1 - N
	}
	xi := float64(istrt)
	powx := new(Result)
	stat_powx := Pow_int_e(x, istrt, powx)
	sarg := beps * gsl.Pi
	sfact := 1.0
	if sarg != 0.0 {
		sfact = sarg / math.Sin(sarg)
	}
	factor_val := sfact
	if gsl.IsOdd(N) {
		factor_val *= -1.0
	}
	factor_val *= powx.val
	factor_err := math.Abs(powx.err) + 2.0*gsl.Float64Eps*math.Abs(factor_val)

	pochai, gamri1, gamrni := new(Result), new(Result), new(Result)
	stat_pochai := Poch_e(a, xi, pochai)
	stat_gamri1 := Gammainv_e(xi+1.0, gamri1)
	stat_gamrni := Gammainv_e(bint+xi, gamrni)
	stat_gam123 := err.ErrorSelect(stat_gamri1, stat_gamrni)
	stat_gamall := err.ErrorSelect(stat_gam123, stat_pochai, stat_powx)

	pochaxibeps, gamrxi1beps := new(Result), new(Result)
	stat_pochaxibeps := Poch_e(a, xi-beps, pochaxibeps)
	stat_gamrxi1beps := Gammainv_e(xi+1.0-beps, gamrxi1beps)

	stat_all := err.ErrorSelect(stat_gamall, stat_pochaxibeps, stat_gamrxi1beps)

	X := sfact
	if gsl.IsOdd(N) {
		X *= -1.0
	}
	X *= powx.val * Poch(1+a-b, xi-1+b-beps) * Gammainv(a)

	b0_val := X * gamrni.val * gamrxi1beps.val
	b0_err := math.Abs(factor_val*pochaxibeps.val*gamrni.val)*gamrxi1beps.err + math.Abs(factor_val*pochaxibeps.val*gamrxi1beps.val)*gamrni.err + math.Abs(factor_val*gamrni.val*gamrxi1beps.val)*pochaxibeps.err + math.Abs(pochaxibeps.val*gamrni.val*gamrxi1beps.val)*factor_err + 2.0*gsl.Float64Eps*math.Abs(b0_val)

	/*
	   C  X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE
	   C  STRAIGHTFORWARD FORMULATION IS STABLE.
	*/
	var (
		i        int
		dchu_val float64
		dchu_err float64
		t_val    float64
		t_err    float64
		gamr     Result
		dgamrbxi Result
	)
	stat_gamr := Gammainv_e(1.0+a-b, &gamr)
	stat_dgamrbxi := Gammainv_e(b+xi, &dgamrbxi)
	a0_val := factor_val * gamr.val * pochai.val * dgamrbxi.val * gamri1.val / beps
	a0_err := math.Abs(factor_val*pochai.val*dgamrbxi.val*gamri1.val/beps)*gamr.err + math.Abs(factor_val*gamr.val*dgamrbxi.val*gamri1.val/beps)*pochai.err + math.Abs(factor_val*gamr.val*pochai.val*gamri1.val/beps)*dgamrbxi.err + math.Abs(factor_val*gamr.val*pochai.val*dgamrbxi.val/beps)*gamri1.err + math.Abs(pochai.val*gamr.val*dgamrbxi.val*gamri1.val/beps)*factor_err + 2.0*gsl.Float64Eps*math.Abs(a0_val)
	stat_all = err.ErrorSelect(stat_all, stat_gamr, stat_dgamrbxi)

	b0_val = xeps * b0_val / beps
	b0_err = math.Abs(xeps/beps)*b0_err + 4.0*gsl.Float64Eps*math.Abs(b0_val)
	dchu_val = sum.val + a0_val - b0_val
	dchu_err = sum.err + a0_err + b0_err + 2.0*gsl.Float64Eps*(math.Abs(sum.val)+math.Abs(a0_val)+math.Abs(b0_val))
	for i = 1; i < 2000; i++ {
		xi := float64(istrt + i)
		xi1 := float64(istrt + i - 1)
		a0_multiplier := (a + xi1) * x / ((b + xi1) * xi)
		b0_multiplier := (a + xi1 - beps) * x / ((bint + xi1) * (xi - beps))
		a0_val *= a0_multiplier
		a0_err += math.Abs(a0_multiplier) * a0_err
		b0_val *= b0_multiplier
		b0_err += math.Abs(b0_multiplier) * b0_err
		t_val = a0_val - b0_val
		t_err = a0_err + b0_err
		dchu_val += t_val
		dchu_err += t_err
		if sys.Finite(t_val) == 0 || math.Abs(t_val) < EPS*math.Abs(dchu_val) {
			break
		}
	}

	result.val = dchu_val
	result.err = 2.0 * dchu_err
	result.err += 2.0 * math.Abs(t_val)
	result.err += 4.0 * gsl.Float64Eps * (float64(i) + 2.0) * math.Abs(dchu_val)
	result.err *= 2.0 /* FIXME: fudge factor */

	if i >= 2000 {
		return err.ERROR("error", err.EMAXITER)
	}

	return stat_all

}

func hyperg_U_infinite_sum_improved(N int, a, bint, b, beps, x, xeps float64, sum Result, result *Result) err.GSLError {
	EPS := 2.0 * gsl.Float64Eps /* EPS = D1MACH(3) */
	lnx := math.Log(x)

	var istrt int
	if N < 1 {
		istrt = 1 - N
	}
	xi := float64(istrt)

	gamr, powx := new(Result), new(Result)
	stat_gamr := Gammainv_e(1.0+a-b, gamr)
	stat_powx := Pow_int_e(x, istrt, powx)
	sarg := beps * gsl.Pi
	sfact := 1.0
	if sarg != 0.0 {
		sfact = sarg / math.Sin(sarg)
	}

	factor_val := sfact
	if gsl.IsOdd(N) {
		factor_val *= -1.0
	}
	factor_val *= gamr.val * powx.val
	factor_err := math.Abs(gamr.val)*powx.err + math.Abs(powx.val)*gamr.err + 2.0*gsl.Float64Eps*math.Abs(factor_val)

	pochai, gamri1, gamrni := new(Result), new(Result), new(Result)
	stat_pochai := Poch_e(a, xi, pochai)
	stat_gamri1 := Gammainv_e(xi+1.0, gamri1)
	stat_gamrni := Gammainv_e(bint+xi, gamrni)
	stat_gam123 := err.ErrorSelect(stat_gamr, stat_gamri1, stat_gamrni)
	stat_gamall := err.ErrorSelect(stat_gam123, stat_pochai, stat_powx)

	pochaxibeps, gamrxi1beps := new(Result), new(Result)
	stat_pochaxibeps := Poch_e(a, xi-beps, pochaxibeps)
	stat_gamrxi1beps := Gammainv_e(xi+1.0-beps, gamrxi1beps)

	stat_all := err.ErrorSelect(stat_gamall, stat_pochaxibeps, stat_gamrxi1beps)

	b0_val := factor_val * pochaxibeps.val * gamrni.val * gamrxi1beps.val
	b0_err := math.Abs(factor_val*pochaxibeps.val*gamrni.val)*gamrxi1beps.err + math.Abs(factor_val*pochaxibeps.val*gamrxi1beps.val)*gamrni.err + math.Abs(factor_val*gamrni.val*gamrxi1beps.val)*pochaxibeps.err + math.Abs(pochaxibeps.val*gamrni.val*gamrxi1beps.val)*factor_err + 2.0*gsl.Float64Eps*math.Abs(b0_val)

	/*
	   C  X**(-BEPS) IS CLOSE TO 1.0D0, SO WE MUST BE
	   C  CAREFUL IN EVALUATING THE DIFFERENCES.
	*/
	var (
		i                           int
		pch1ai, pch1i, poch1bxibeps Result
	)

	stat_pch1ai := Pochrel_e(a+xi, -beps, &pch1ai)
	stat_pch1i := Pochrel_e(xi+1.0-beps, beps, &pch1i)
	stat_poch1bxibeps := Pochrel_e(b+xi, -beps, &poch1bxibeps)
	c0_t1_val := beps * pch1ai.val * pch1i.val
	c0_t1_err := math.Abs(beps)*math.Abs(pch1ai.val)*pch1i.err + math.Abs(beps)*math.Abs(pch1i.val)*pch1ai.err + 2.0*gsl.Float64Eps*math.Abs(c0_t1_val)
	c0_t2_val := -poch1bxibeps.val + pch1ai.val - pch1i.val + c0_t1_val
	c0_t2_err := poch1bxibeps.err + pch1ai.err + pch1i.err + c0_t1_err + 2.0*gsl.Float64Eps*math.Abs(c0_t2_val)
	c0_val := factor_val * pochai.val * gamrni.val * gamri1.val * c0_t2_val
	c0_err := math.Abs(factor_val*pochai.val*gamrni.val*gamri1.val)*c0_t2_err + math.Abs(factor_val*pochai.val*gamrni.val*c0_t2_val)*gamri1.err + math.Abs(factor_val*pochai.val*gamri1.val*c0_t2_val)*gamrni.err + math.Abs(factor_val*gamrni.val*gamri1.val*c0_t2_val)*pochai.err + math.Abs(pochai.val*gamrni.val*gamri1.val*c0_t2_val)*factor_err + 2.0*gsl.Float64Eps*math.Abs(c0_val)
	/*
	   C  XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS)
	*/
	dexprl := new(Result)
	stat_dexprl := Exprel_e(-beps*lnx, dexprl)
	xeps1_val := lnx * dexprl.val
	xeps1_err := 2.0*gsl.Float64Eps*(1.0+math.Abs(beps*lnx))*math.Abs(dexprl.val) + math.Abs(lnx)*dexprl.err + 2.0*gsl.Float64Eps*math.Abs(xeps1_val)
	dchu_val := sum.val + c0_val + xeps1_val*b0_val
	dchu_err := sum.err + c0_err + math.Abs(xeps1_val)*b0_err + xeps1_err*math.Abs(b0_val) + math.Abs(b0_val*lnx)*dexprl.err + 2.0*gsl.Float64Eps*(math.Abs(sum.val)+math.Abs(c0_val)+math.Abs(xeps1_val*b0_val))
	xn := float64(N)
	var t_val, t_err float64

	stat_all = err.ErrorSelect(stat_all, stat_dexprl, stat_poch1bxibeps, stat_pch1i, stat_pch1ai)

	for i = 1; i < 2000; i++ {
		xi := float64(istrt + i)
		xi1 := float64(istrt + i - 1)
		tmp := (a-1.0)*(xn+2.0*xi-1.0) + xi*(xi-beps)
		b0_multiplier := (a + xi1 - beps) * x / ((xn + xi1) * (xi - beps))
		c0_multiplier_1 := (a + xi1) * x / ((b + xi1) * xi)
		c0_multiplier_2 := tmp / (xi * (b + xi1) * (a + xi1 - beps))
		b0_val *= b0_multiplier
		b0_err += math.Abs(b0_multiplier)*b0_err + math.Abs(b0_val)*8.0*2.0*gsl.Float64Eps
		c0_val = c0_multiplier_1*c0_val - c0_multiplier_2*b0_val
		c0_err = math.Abs(c0_multiplier_1)*c0_err + math.Abs(c0_multiplier_2)*b0_err + math.Abs(c0_val)*8.0*2.0*gsl.Float64Eps + math.Abs(b0_val*c0_multiplier_2)*16.0*2.0*gsl.Float64Eps
		t_val = c0_val + xeps1_val*b0_val
		t_err = c0_err + math.Abs(xeps1_val)*b0_err
		t_err += math.Abs(b0_val*lnx) * dexprl.err
		t_err += math.Abs(b0_val) * xeps1_err
		dchu_val += t_val
		dchu_err += t_err
		if math.Abs(t_val) < EPS*math.Abs(dchu_val) {
			break
		}
	}

	result.val = dchu_val
	result.err = 2.0 * dchu_err
	result.err += 2.0 * math.Abs(t_val)
	result.err += 4.0 * gsl.Float64Eps * (float64(i) + 2.0) * math.Abs(dchu_val)
	result.err *= 2.0 /* FIXME: fudge factor */

	if i >= 2000 {
		return err.ERROR("error", err.EMAXITER)
	}

	return stat_all

}

/* Based on SLATEC DCHU() [W. Fullerton]
 * Assumes x > 0.
 * This is just a series summation method, and
 * it is not good for large a.
 *
 * I patched up the window for 1+a-b near zero. [GJ]
 */
func hyperg_U_series(a, b, x float64, result *Result) err.GSLError {
	SQRT_EPS := math.Sqrt2 * gsl.SqrtFloat64Eps
	bint := math.Floor(b + 0.5)
	if b < 0.0 {
		bint = math.Ceil(b - 0.5)
	}
	beps := b - bint
	a_beps := a - beps
	r_a_beps := math.Floor(a_beps + 0.5)
	a_beps_int := (math.Abs(a_beps-r_a_beps) < INT_THRESHOLD)
	/*  double a_b_1 = a-b+1;
	    double r_a_b_1 = math.Floor(a_b_1+0.5);
	    double r_a_b_1_int = (math.Abs(a_b_1-r_a_b_1)< INT_THRESHOLD);
	    Check for (a-beps) being a member of -N; N being 0,1,... */
	if a_beps_int && a_beps <= 0 {
		beps = beps - 1 + math.Floor(a_beps)
		bint = bint + 1 - math.Floor(a_beps)
	}

	if math.Abs(1.0+a-b) < SQRT_EPS {
		/* Original Comment: ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X
		 */
		/* We can however do the following:
		 * U(a,b,x) = U(a,a+1,x) when 1+a-b=0
		 * and U(a,a+1,x) = x^(-a).
		 */
		lnr := -a * math.Log(x)
		stat_e := Exp_e(lnr, result)
		result.err += 2.0 * SQRT_EPS * math.Abs(result.val)
		return stat_e
	} else {
		N := int(bint)

		lnx := math.Log(x)
		xeps := math.Exp(-beps * lnx)

		/* Evaluate finite sum.
		 */
		var sum Result
		stat_sum := hyperg_U_finite_sum(N, a, b, x, xeps, &sum)
		var stat_inf err.GSLError

		/* Evaluate infinite sum. */
		if math.Abs(xeps-1.0) > 0.5 {
			stat_inf = hyperg_U_infinite_sum_stable(N, a, bint, b, beps, x, xeps, sum, result)
		} else if 1+a-b < 0 && 1+a-b == math.Floor(1+a-b) && beps != 0 {
			stat_inf = hyperg_U_infinite_sum_simple(N, a, bint, b, beps, x, xeps, sum, result)
		} else {
			stat_inf = hyperg_U_infinite_sum_improved(N, a, bint, b, beps, x, xeps, sum, result)
		}

		return err.ErrorSelect(stat_sum, stat_inf)

	}
}

/* Assumes b > 0 and x > 0.
 */
func hyperg_U_small_ab(a, b, x float64, result *Result) err.GSLError {
	if a == -1.0 {
		/* U(-1,c+1,x) = Laguerre[c,0,x] = -b + x
		 */
		result.val = -b + x
		result.err = 2.0 * gsl.Float64Eps * (math.Abs(b) + math.Abs(x))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if a == 0.0 {
		/* U(0,c+1,x) = Laguerre[c,0,x] = 1
		 */
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if AsympEvalOk(a, b, x) {
		p := math.Pow(x, -a)
		asymp := new(Result)
		stat_asymp := hyperg_zaU_asymp(a, b, x, asymp)
		result.val = asymp.val * p
		result.err = asymp.err * p
		result.err += math.Abs(asymp.val) * gsl.Float64Eps * math.Abs(a) * p
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_asymp
	} else {
		return hyperg_U_series(a, b, x, result)
	}
}

/* Assumes b > 0 and x > 0.
 */
func hyperg_U_small_a_bgt0(a, b, x float64, result *Result, ln_multiplier *float64) err.GSLError {
	if a == 0.0 {
		result.val = 1.0
		result.err = 0.0
		*ln_multiplier = 0.0
		return nil
	} else if (b > 5000.0 && x < 0.90*math.Abs(b)) || (b > 500.0 && x < 0.50*math.Abs(b)) {
		stat := Hyperg_U_large_b_e(a, b, x, result, ln_multiplier)
		if stat != nil {
			if stat.Status() == err.EOVRFLW {
				return nil
			}
		}

		return stat

	} else if b > 15.0 {
		/* Recurse up from b near 1.
		 */
		eps := b - math.Floor(b)
		b0 := 1.0 + eps
		var r_Ubm1, r_Ub Result
		stat_0 := hyperg_U_small_ab(a, b0, x, &r_Ubm1)
		stat_1 := hyperg_U_small_ab(a, b0+1.0, x, &r_Ub)
		Ubm1 := r_Ubm1.val
		Ub := r_Ub.val

		for bp := b0 + 1.0; bp < b-0.1; bp += 1.0 {
			Ubp1 := ((1.0+a-bp)*Ubm1 + (bp+x-1.0)*Ub) / x
			Ubm1 = Ub
			Ub = Ubp1
		}
		result.val = Ub
		result.err = (math.Abs(r_Ubm1.err/r_Ubm1.val) + math.Abs(r_Ub.err/r_Ub.val)) * math.Abs(Ub)
		result.err += 2.0 * gsl.Float64Eps * (math.Abs(b-b0) + 1.0) * math.Abs(Ub)
		*ln_multiplier = 0.0
		return err.ErrorSelect(stat_0, stat_1)
	} else {
		*ln_multiplier = 0.0
		return hyperg_U_small_ab(a, b, x, result)
	}
}

/* We use this to keep track of large
 * dynamic ranges in the recursions.
 * This can be important because sometimes
 * we want to calculate a very large and
 * a very small number and the answer is
 * the product, of order 1. This happens,
 * for instance, when we apply a Kummer
 * transform to make b positive and
 * both x and b are large.
 */
func reScale2(u0, u1, factor *float64, count *int) {
	au0 := math.Abs(*u0)
	if au0 > *factor {
		*u0 /= *factor
		*u1 /= *factor
		*count++
	} else if au0 < 1.0/(*factor) {
		*u0 *= *factor
		*u1 *= *factor
		*count--
	}
}

/* Specialization to b >= 1, for integer parameters.
 * Assumes x > 0.
 */
func hyperg_U_int_bge1(a, b int, x float64, result *Result_e10) err.GSLError {
	if a == 0 {
		result.val = 1.0
		result.err = 0.0
		result.e10 = 0
		return nil
	} else if a == -1 {
		result.val = -float64(b) + x
		result.err = 2.0 * gsl.Float64Eps * (math.Abs(float64(b)) + math.Abs(x))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		result.e10 = 0
		return nil
	} else if b == a+1 {
		/* U(a,a+1,x) = x^(-a)
		 */
		return Exp_e10_e(-float64(a)*math.Log(x), result)
	} else if AsympEvalOk(float64(a), float64(b), x) {
		ln_pre_val := -float64(a) * math.Log(x)
		ln_pre_err := 2.0 * gsl.Float64Eps * math.Abs(ln_pre_val)
		var asymp Result
		stat_asymp := hyperg_zaU_asymp(float64(a), float64(b), x, &asymp)
		stat_e := Exp_mult_err_e10_e(ln_pre_val, ln_pre_err, asymp.val, asymp.err, result)
		return err.ErrorSelect(stat_e, stat_asymp)
	} else if seriesEvalOk(float64(a), float64(b), x) && 1+a-b > 0 {
		var ser Result
		stat_ser := hyperg_U_series(float64(a), float64(b), x, &ser)
		result.val = ser.val
		result.err = ser.err
		result.e10 = 0
		return stat_ser
	} else if a < 0 {
		/* Recurse backward from a = -1,0.
		 */
		scale_count := 0
		scale_factor := gsl.SqrtMaxFloat64
		var lnm, y Result
		var lnscale float64
		Uap1 := 1.0           /* U(0,b,x)  */
		Ua := -float64(b) + x /* U(-1,b,x) */

		for ap := -1; ap > a; ap-- {
			fb := float64(b)
			fap := float64(ap)
			Uam1 := fap*(fb-fap-1.0)*Uap1 + (x+2.0*fap-fb)*Ua
			Uap1 = Ua
			Ua = Uam1
			reScale2(&Ua, &Uap1, &scale_factor, &scale_count)
		}

		lnscale = math.Log(scale_factor)
		lnm.val = float64(scale_count) * lnscale
		lnm.err = 2.0 * gsl.Float64Eps * math.Abs(lnm.val)
		y.val = Ua
		y.err = 4.0 * gsl.Float64Eps * (math.Abs(float64(a)) + 1.0) * math.Abs(Ua)
		return Exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result)
	} else if float64(b) >= 2.0*float64(a)+x {
		/* Recurse forward from a = 0,1.
		 */
		scale_count := 0
		scale_factor := gsl.SqrtMaxFloat64
		var r_Ua, lnm, y Result
		var lnscale, lm float64
		stat_1 := hyperg_U_small_a_bgt0(1.0, float64(b), x, &r_Ua, &lm) /* U(1,b,x) */
		var stat_e err.GSLError
		Uam1 := 1.0 /* U(0,b,x) */
		Ua := r_Ua.val
		Uam1 *= math.Exp(-lm)

		for ap := 1; ap < a; ap++ {
			fb := float64(b)
			fap := float64(ap)
			Uap1 := -(Uam1 + (fb-2.0*fap-x)*Ua) / (fap * (1.0 + fap - fb))
			Uam1 = Ua
			Ua = Uap1
			reScale2(&Ua, &Uam1, &scale_factor, &scale_count)
		}

		lnscale = math.Log(scale_factor)
		lnm.val = lm + float64(scale_count)*lnscale
		lnm.err = 2.0 * gsl.Float64Eps * (math.Abs(lm) + math.Abs(float64(scale_count)*lnscale))
		y.val = Ua
		y.err = math.Abs(r_Ua.err/r_Ua.val) * math.Abs(Ua)
		y.err += 2.0 * gsl.Float64Eps * (math.Abs(float64(a)) + 1.0) * math.Abs(Ua)
		stat_e = Exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result)
		return err.ErrorSelect(stat_e, stat_1)
	} else {
		if float64(b) <= x {
			/* Recurse backward either to the b=a+1 line
			 * or to a=0, whichever we hit.
			 */
			scale_factor := gsl.SqrtMaxFloat64
			scale_count := 0
			var stat_CF1 err.GSLError
			var ru float64
			var CF1_count int
			var a_target int
			var lnU_target, Ua, Uap1, Uam1 float64

			if b < a+1 {
				a_target = b - 1
				lnU_target = -float64(a_target) * math.Log(x)
			} else {
				a_target = 0
				lnU_target = 0.0
			}

			stat_CF1 = hyperg_U_CF1(float64(a), float64(b), 0, x, &ru, &CF1_count)

			Ua = 1.0
			Uap1 = ru / float64(a) * Ua
			for ap := a; ap > a_target; ap-- {
				fb := float64(b)
				fap := float64(ap)
				Uam1 = -((fb-2.0*fap-x)*Ua + fap*(1.0+fap-fb)*Uap1)
				Uap1 = Ua
				Ua = Uam1
				reScale2(&Ua, &Uap1, &scale_factor, &scale_count)
			}

			if Ua == 0.0 {
				result.val = 0.0
				result.err = 0.0
				result.e10 = 0
				return err.ERROR("error", err.EZERODIV)
			} else {
				lnscl := -float64(scale_count) * math.Log(scale_factor)
				lnpre_val := lnU_target + lnscl
				lnpre_err := 2.0 * gsl.Float64Eps * (math.Abs(lnU_target) + math.Abs(lnscl))
				oUa_err := 2.0 * (math.Abs(float64(a_target-a)) + float64(CF1_count) + 1.0) * gsl.Float64Eps * math.Abs(1.0/Ua)
				stat_e := Exp_mult_err_e10_e(lnpre_val, lnpre_err, 1.0/Ua, oUa_err, result)
				return err.ErrorSelect(stat_e, stat_CF1)
			}
		} else {
			/* Recurse backward to near the b=2a+x line, then
			 * determine normalization by either direct evaluation
			 * or by a forward recursion. The direct evaluation
			 * is needed when x is small (which is precisely
			 * when it is easy to do).
			 */
			scale_factor := gsl.SqrtMaxFloat64
			scale_count_for := 0
			scale_count_bck := 0
			a0 := 1
			a1 := a0 + int(math.Ceil(0.5*(float64(b)-x)-float64(a0)))
			var Ua1_bck_val, Ua1_bck_err, Ua1_for_val, Ua1_for_err float64
			var stat_for, stat_bck err.GSLError
			var lm_for Result

			{
				/* Recurse back to determine U(a1,b), sans normalization.
				 */
				var (
					ru        float64
					CF1_count int
					stat_CF1  = hyperg_U_CF1(float64(a), float64(b), 0, x, &ru, &CF1_count)
					Ua        = 1.0
					Uap1      = ru / float64(a) * Ua
				)
				for ap := a; ap > a1; ap-- {
					fb := float64(b)
					fap := float64(ap)
					Uam1 := -((fb-2.0*fap-x)*Ua + fap*(1.0+fap-fb)*Uap1)
					Uap1 = Ua
					Ua = Uam1
					reScale2(&Ua, &Uap1, &scale_factor, &scale_count_bck)
				}
				Ua1_bck_val = Ua
				Ua1_bck_err = 2.0 * gsl.Float64Eps * (math.Abs(float64(a1-a)) + float64(CF1_count) + 1.0) * math.Abs(Ua)
				stat_bck = stat_CF1
			}

			if b == 2*a1 && a1 > 1 {
				/* This can happen when x is small, which is
				 * precisely when we need to be careful with
				 * this evaluation.
				 */
				hyperg_lnU_beq2a(float64(a1), x, &lm_for)
				Ua1_for_val = 1.0
				Ua1_for_err = 0.0
				stat_for = nil
			} else if b == 2*a1-1 && a1 > 1 {
				/* Similar to the above. Happens when x is small.
				 * Use
				 *   U(a,2a-1) = (x U(a,2a) - U(a-1,2(a-1))) / (2a - 2)
				 */
				var lnU00, lnU12, U00, U12 Result
				hyperg_lnU_beq2a(float64(a1)-1.0, x, &lnU00)
				hyperg_lnU_beq2a(float64(a1), x, &lnU12)
				if lnU00.val > lnU12.val {
					lm_for.val = lnU00.val
					lm_for.err = lnU00.err
					U00.val = 1.0
					U00.err = 0.0
					Exp_err_e(lnU12.val-lm_for.val, lnU12.err+lm_for.err, &U12)
				} else {
					lm_for.val = lnU12.val
					lm_for.err = lnU12.err
					U12.val = 1.0
					U12.err = 0.0
					Exp_err_e(lnU00.val-lm_for.val, lnU00.err+lm_for.err, &U00)
				}
				Ua1_for_val = (x*U12.val - U00.val) / (2.0*float64(a1) - 2.0)
				Ua1_for_err = (math.Abs(x)*U12.err + U00.err) / math.Abs(2.0*float64(a1)-2.0)
				Ua1_for_err += 2.0 * gsl.Float64Eps * math.Abs(Ua1_for_val)
				stat_for = nil
			} else {
				/* Recurse forward to determine U(a1,b) with
				 * absolute normalization.
				 */
				r_Ua := new(Result)
				Uam1 := 1.0 /* U(a0-1,b,x) = U(0,b,x) */
				var lm_for_local float64
				stat_for = hyperg_U_small_a_bgt0(float64(a0), float64(b), x, r_Ua, &lm_for_local) /* U(1,b,x) */
				Ua := r_Ua.val
				Uam1 *= math.Exp(-lm_for_local)
				lm_for.val = lm_for_local
				lm_for.err = 0.0

				for ap := a0; ap < a1; ap++ {
					fb := float64(b)
					fap := float64(ap)
					Uap1 := -(Uam1 + (fb-2.0*fap-x)*Ua) / (fap * (1.0 + fap - fb))
					Uam1 = Ua
					Ua = Uap1
					reScale2(&Ua, &Uam1, &scale_factor, &scale_count_for)
				}
				Ua1_for_val = Ua
				Ua1_for_err = math.Abs(Ua) * math.Abs(r_Ua.err/r_Ua.val)
				Ua1_for_err += 2.0 * gsl.Float64Eps * (math.Abs(float64(a1-a0)) + 1.0) * math.Abs(Ua1_for_val)
			}

			/* Now do the matching to produce the final result.
			 */
			if Ua1_bck_val == 0.0 {
				result.val = 0.0
				result.err = 0.0
				result.e10 = 0
				return err.ERROR("error", err.EZERODIV)
			} else if Ua1_for_val == 0.0 {
				/* Should never happen. */
				return UnderflowError_e10(result)
			} else {
				lns := float64(scale_count_for-scale_count_bck) * math.Log(scale_factor)
				ln_for_val := math.Log(math.Abs(Ua1_for_val))
				ln_for_err := gsl.Float64Eps + math.Abs(Ua1_for_err/Ua1_for_val)
				ln_bck_val := math.Log(math.Abs(Ua1_bck_val))
				ln_bck_err := gsl.Float64Eps + math.Abs(Ua1_bck_err/Ua1_bck_val)
				lnr_val := lm_for.val + ln_for_val - ln_bck_val + lns
				lnr_err := lm_for.err + ln_for_err + ln_bck_err + 2.0*gsl.Float64Eps*(math.Abs(lm_for.val)+math.Abs(ln_for_val)+math.Abs(ln_bck_val)+math.Abs(lns))
				sgn := gsl.Sign(Ua1_for_val) * gsl.Sign(Ua1_bck_val)
				stat_e := Exp_err_e10_e(lnr_val, lnr_err, result)
				result.val *= sgn
				return err.ErrorSelect(stat_e, stat_bck, stat_for)
			}
		}
	}
}

/* Handle b >= 1 for generic a,b values.
 */
func hyperg_U_bge1(a, b, x float64, result *Result_e10) err.GSLError {
	rinta := math.Floor(a + 0.5)
	a_neg_integer := (a < 0.0 && math.Abs(a-rinta) < INT_THRESHOLD)

	if a == 0.0 {
		result.val = 1.0
		result.err = 0.0
		result.e10 = 0
		return nil
	} else if a_neg_integer && math.Abs(rinta) < math.MaxInt32 {
		/* U(-n,b,x) = (-1)^n n! Laguerre[n,b-1,x]
		 */
		n := -int(rinta)
		sgn := 1.0
		if gsl.IsOdd(n) {
			sgn = -1.0
		}
		var lnfact, L Result
		stat_L := Laguerre_n_e(n, b-1.0, x, &L)
		Lnfact_e(uint(n), &lnfact)
		{
			stat_e := Exp_mult_err_e10_e(lnfact.val, lnfact.err, sgn*L.val, L.err, result)
			return err.ErrorSelect(stat_e, stat_L)
		}
	} else if AsympEvalOk(a, b, x) {
		ln_pre_val := -a * math.Log(x)
		ln_pre_err := 2.0 * gsl.Float64Eps * math.Abs(ln_pre_val)
		var asymp Result
		stat_asymp := hyperg_zaU_asymp(a, b, x, &asymp)
		stat_e := Exp_mult_err_e10_e(ln_pre_val, ln_pre_err, asymp.val, asymp.err, result)
		return err.ErrorSelect(stat_e, stat_asymp)
	} else if math.Abs(a) <= 1.0 {
		var rU Result
		var ln_multiplier float64
		stat_U := hyperg_U_small_a_bgt0(a, b, x, &rU, &ln_multiplier)
		stat_e := Exp_mult_err_e10_e(ln_multiplier, 2.0*gsl.Float64Eps*math.Abs(ln_multiplier), rU.val, rU.err, result)
		return err.ErrorSelect(stat_U, stat_e)
	} else if seriesEvalOk(a, b, x) {
		var ser Result
		stat_ser := hyperg_U_series(a, b, x, &ser)
		result.val = ser.val
		result.err = ser.err
		result.e10 = 0
		return stat_ser
	} else if a < 0.0 {
		/* Recurse backward on a and then upward on b.
		 */
		scale_factor := gsl.SqrtMaxFloat64
		a0 := a - math.Floor(a) - 1.0
		b0 := b - math.Floor(b) + 1.0
		scale_count := 0
		var lm_0, lm_1, lm_max float64
		var r_Uap1, r_Ua Result
		stat_0 := hyperg_U_small_a_bgt0(a0+1.0, b0, x, &r_Uap1, &lm_0)
		stat_1 := hyperg_U_small_a_bgt0(a0, b0, x, &r_Ua, &lm_1)
		var stat_e err.GSLError
		Uap1 := r_Uap1.val
		Ua := r_Ua.val
		lm_max = math.Max(lm_0, lm_1)
		Uap1 *= math.Exp(lm_0 - lm_max)
		Ua *= math.Exp(lm_1 - lm_max)

		/* Downward recursion on a.
		 */
		for ap := a0; ap > a+0.1; ap -= 1.0 {
			Uam1 := ap*(b0-ap-1.0)*Uap1 + (x+2.0*ap-b0)*Ua
			Uap1 = Ua
			Ua = Uam1
			reScale2(&Ua, &Uap1, &scale_factor, &scale_count)
		}

		if b < 2.0 {
			/* b == b0, so no recursion necessary
			 */
			lnscale := math.Log(scale_factor)
			var lnm, y Result
			lnm.val = lm_max + float64(scale_count)*lnscale
			lnm.err = 2.0 * gsl.Float64Eps * (math.Abs(lm_max) + float64(scale_count)*math.Abs(lnscale))
			y.val = Ua
			y.err = math.Abs(r_Uap1.err/r_Uap1.val) * math.Abs(Ua)
			y.err += math.Abs(r_Ua.err/r_Ua.val) * math.Abs(Ua)
			y.err += 2.0 * gsl.Float64Eps * (math.Abs(a-a0) + 1.0) * math.Abs(Ua)
			y.err *= math.Abs(lm_0-lm_max) + 1.0
			y.err *= math.Abs(lm_1-lm_max) + 1.0
			stat_e = Exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result)
		} else {
			/* Upward recursion on b.
			 */
			err_mult := math.Abs(b-b0) + math.Abs(a-a0) + 1.0
			lnscale := math.Log(scale_factor)
			var lnm, y Result

			Ubm1 := Ua                               /* U(a,b0)   */
			Ub := (a*(b0-a-1.0)*Uap1 + (a+x)*Ua) / x /* U(a,b0+1) */

			for bp := b0 + 1.0; bp < b-0.1; bp += 1.0 {
				Ubp1 := ((1.0+a-bp)*Ubm1 + (bp+x-1.0)*Ub) / x
				Ubm1 = Ub
				Ub = Ubp1
				reScale2(&Ub, &Ubm1, &scale_factor, &scale_count)
			}

			lnm.val = lm_max + float64(scale_count)*lnscale
			lnm.err = 2.0 * gsl.Float64Eps * (math.Abs(lm_max) + math.Abs(float64(scale_count)*lnscale))
			y.val = Ub
			y.err = 2.0 * err_mult * math.Abs(r_Uap1.err/r_Uap1.val) * math.Abs(Ub)
			y.err += 2.0 * err_mult * math.Abs(r_Ua.err/r_Ua.val) * math.Abs(Ub)
			y.err += 2.0 * gsl.Float64Eps * err_mult * math.Abs(Ub)
			y.err *= math.Abs(lm_0-lm_max) + 1.0
			y.err *= math.Abs(lm_1-lm_max) + 1.0
			stat_e = Exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result)
		}
		return err.ErrorSelect(stat_e, stat_0, stat_1)
	} else if b >= 2*a+x {
		/* Recurse forward from a near zero.
		 * Note that we cannot cross the singularity at
		 * the line b=a+1, because the only way we could
		 * be in that little wedge is if a < 1. But we
		 * have already dealt with the small a case.
		 */
		scale_count := 0
		a0 := a - math.Floor(a)
		scale_factor := gsl.SqrtMaxFloat64
		var lnscale, lm_0, lm_1, lm_max float64
		var r_Uam1, r_Ua Result
		stat_0 := hyperg_U_small_a_bgt0(a0-1.0, b, x, &r_Uam1, &lm_0)
		stat_1 := hyperg_U_small_a_bgt0(a0, b, x, &r_Ua, &lm_1)
		var stat_e err.GSLError
		var lnm, y Result
		Uam1 := r_Uam1.val
		Ua := r_Ua.val
		lm_max = math.Max(lm_0, lm_1)
		Uam1 *= math.Exp(lm_0 - lm_max)
		Ua *= math.Exp(lm_1 - lm_max)

		for ap := a0; ap < a-0.1; ap += 1.0 {
			Uap1 := -(Uam1 + (b-2.0*ap-x)*Ua) / (ap * (1.0 + ap - b))
			Uam1 = Ua
			Ua = Uap1
			reScale2(&Ua, &Uam1, &scale_factor, &scale_count)
		}

		lnscale = math.Log(scale_factor)
		lnm.val = lm_max + float64(scale_count)*lnscale
		lnm.err = 2.0 * gsl.Float64Eps * (math.Abs(lm_max) + math.Abs(float64(scale_count)*lnscale))
		y.val = Ua
		y.err = math.Abs(r_Uam1.err/r_Uam1.val) * math.Abs(Ua)
		y.err += math.Abs(r_Ua.err/r_Ua.val) * math.Abs(Ua)
		y.err += 2.0 * gsl.Float64Eps * (math.Abs(a-a0) + 1.0) * math.Abs(Ua)
		y.err *= math.Abs(lm_0-lm_max) + 1.0
		y.err *= math.Abs(lm_1-lm_max) + 1.0
		stat_e = Exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result)
		return err.ErrorSelect(stat_e, stat_0, stat_1)
	} else {
		if b <= x {
			/* Recurse backward to a near zero.
			 */
			a0 := a - math.Floor(a)
			scale_factor := gsl.SqrtMaxFloat64
			scale_count := 0
			var lnm, y Result
			var lnscale, lm_0, Uap1, Ua, Uam1 float64
			var U0 Result
			var ap, ru, r float64
			var CF1_count int
			stat_CF1 := hyperg_U_CF1(a, b, 0, x, &ru, &CF1_count)
			var stat_U0, stat_e err.GSLError
			r = ru / a
			Ua = gsl.SqrtMinFloat64
			Uap1 = r * Ua
			for ap = a; ap > a0+0.1; ap -= 1.0 {
				Uam1 = -((b-2.0*ap-x)*Ua + ap*(1.0+ap-b)*Uap1)
				Uap1 = Ua
				Ua = Uam1
				reScale2(&Ua, &Uap1, &scale_factor, &scale_count)
			}

			stat_U0 = hyperg_U_small_a_bgt0(a0, b, x, &U0, &lm_0)

			lnscale = math.Log(scale_factor)
			lnm.val = lm_0 - float64(scale_count)*lnscale
			lnm.err = 2.0 * gsl.Float64Eps * (math.Abs(lm_0) + math.Abs(float64(scale_count)*lnscale))
			y.val = gsl.SqrtMinFloat64 * (U0.val / Ua)
			y.err = gsl.SqrtMinFloat64 * (U0.err / math.Abs(Ua))
			y.err += 2.0 * gsl.Float64Eps * (math.Abs(a0-a) + float64(CF1_count) + 1.0) * math.Abs(y.val)
			stat_e = Exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result)
			return err.ErrorSelect(stat_e, stat_U0, stat_CF1)
		} else {
			/* Recurse backward to near the b=2a+x line, then
			 * forward from a near zero to get the normalization.
			 */
			scale_count_for := 0
			scale_count_bck := 0
			scale_factor := gsl.SqrtMaxFloat64
			eps := a - math.Floor(a)
			a0 := eps
			if eps == 0.0 {
				a0 = 1.0
			}
			a1 := a0 + math.Ceil(0.5*(b-x)-a0)
			var lnm, y Result
			var lm_for, lnscale, Ua1_bck, Ua1_for float64
			var stat_for, stat_bck, stat_e err.GSLError
			var CF1_count int

			{
				/* Recurse back to determine U(a1,b), sans normalization.
				 */
				var Uap1, Ua, Uam1, ap, ru, r float64
				stat_CF1 := hyperg_U_CF1(a, b, 0, x, &ru, &CF1_count)
				r = ru / a
				Ua = gsl.SqrtMinFloat64
				Uap1 = r * Ua
				for ap = a; ap > a1+0.1; ap -= 1.0 {
					Uam1 = -((b-2.0*ap-x)*Ua + ap*(1.0+ap-b)*Uap1)
					Uap1 = Ua
					Ua = Uam1
					reScale2(&Ua, &Uap1, &scale_factor, &scale_count_bck)
				}
				Ua1_bck = Ua
				stat_bck = stat_CF1
			}
			{
				/* Recurse forward to determine U(a1,b) with
				 * absolute normalization.
				 */
				var r_Uam1, r_Ua Result
				var lm_0, lm_1 float64
				stat_0 := hyperg_U_small_a_bgt0(a0-1.0, b, x, &r_Uam1, &lm_0)
				stat_1 := hyperg_U_small_a_bgt0(a0, b, x, &r_Ua, &lm_1)
				Uam1 := r_Uam1.val
				Ua := r_Ua.val

				lm_for = math.Max(lm_0, lm_1)
				Uam1 *= math.Exp(lm_0 - lm_for)
				Ua *= math.Exp(lm_1 - lm_for)

				for ap := a0; ap < a1-0.1; ap += 1.0 {
					Uap1 := -(Uam1 + (b-2.0*ap-x)*Ua) / (ap * (1.0 + ap - b))
					Uam1 = Ua
					Ua = Uap1
					reScale2(&Ua, &Uam1, &scale_factor, &scale_count_for)
				}
				Ua1_for = Ua
				stat_for = err.ErrorSelect(stat_0, stat_1)
			}

			lnscale = math.Log(scale_factor)
			lnm.val = lm_for + float64(scale_count_for-scale_count_bck)*lnscale
			lnm.err = 2.0 * gsl.Float64Eps * (math.Abs(lm_for) + math.Abs(float64(scale_count_for-scale_count_bck))*math.Abs(lnscale))
			y.val = gsl.SqrtMinFloat64 * Ua1_for / Ua1_bck
			y.err = 2.0 * gsl.Float64Eps * (math.Abs(a-a0) + float64(CF1_count) + 1.0) * math.Abs(y.val)
			stat_e = Exp_mult_err_e10_e(lnm.val, lnm.err, y.val, y.err, result)
			return err.ErrorSelect(stat_e, stat_bck, stat_for)
		}
	}
}

/* Handle U(a,b,x) using AMS 13.1.3 when x = 0.  This requires b<1 to
   avoid a singularity from the term z^(1-b)

   U(a,b,z=0) = (pi/sin(pi*b)) * 1/(gamma(1+a-b)*gamma(b))

   There are a lot of special limiting cases here

   b = 0, positive integer, negative integer
   1+a-b = 0, negative integer

   I haven't implemented these yet - BJG
*/
func hyperg_U_origin(a, b float64, result *Result_e10) err.GSLError {
	var r1, r2 Result
	stat_1 := Gammainv_e(1+a-b, &r1)
	stat_2 := Gammainv_e(b, &r2)
	factor := gsl.Pi / math.Sin(gsl.Pi*b)

	result.val = factor * r1.val * r2.val
	result.err = math.Abs(factor) * (r1.err + r2.err)
	result.e10 = 0

	return err.ErrorSelect(stat_1, stat_2)
}

func hyperg_U_int_origin(a, b int, result *Result_e10) err.GSLError {
	return hyperg_U_origin(float64(a), float64(b), result)
}

/* Calculate U(a,b,x) for x < 0

   Abramowitz and Stegun formula 13.1.3

     U(a,b,x) = (gamma(1-b)/gamma(1+a-b)) M(a,b,x)
                 - z^(1-b) (gamma(1-b)/gamma(a)) M(1+a-b,2-b,x)

   can be transformed into

     U(a,b,x) = poch(1+a-b,-a) M(a,b,x)
                 + z^(1-b) poch(a,-(1+a-b)) M(1+a-b,2-b,x)

   using the reflection formula 6.1.17 and the definition of
   Poch(a,b)=gamma(a+b)/gamma(a).  Our poch function already handles
   the special cases of ratios of gamma functions with negative
   integer argument.

   Note that U(a,b,x) is complex in general for x<0 due to the term
   x^(1-b), but is real when

   1) b is an integer

   4) a is zero or a negative integer so x^(1-b)/gamma(a) is zero.

   For integer b U(a,b,x) is defined as the limit beta.b U(a,beta,x).
   This makes the situation slightly more complicated.

*/
func hyperg_U_negx(a, b, x float64, result *Result_e10) err.GSLError {
	var r1, r2 Result
	var stat_1, stat_2, status err.GSLError
	a_int := (a == math.Floor(a))
	b_int := (b == math.Floor(b))

	var T1, T1_err, T2, T2_err float64

	/* Compute the first term poch(1+a-b) M(a,b,x) */

	if b_int && b <= 0 && !(a_int && a <= 0 && a >= b) {
		/* Need to handle first term as

		   lim_{beta.b} poch(1+a-beta,-a) M(a,beta,x)

		   due to pole in M(a,b,x) for b == 0 or -ve integer

		   We skip this case when a is zero or a negative integer and
		   a>=b because the hypergeometric series terminates before any
		   singular terms
		*/

		/* FIXME: TO BE IMPLEMENTED ! */
		result.val = math.NaN()
		result.err = math.NaN()
		return err.ERROR("limit case integer b <= 0 unimplemented", err.EUNIMPL)
	} else {
		stat_1 = Poch_e(1+a-b, -a, &r1)
		status = stat_1

		if r1.val != 0.0 {
			var Mr1 Result
			stat_Mr1 := Hyperg_1F1_e(a, b, x, &Mr1)
			status = err.ErrorSelect(status, stat_Mr1)

			T1 = Mr1.val * r1.val
			T1_err = 2.0*gsl.Float64Eps*math.Abs(T1) + math.Abs(Mr1.err*r1.val) + math.Abs(Mr1.val*r1.err)
		}
	}

	/* Compute the second term z^(1-b) poch(a,-(1+a-b)) M(1+a-b,2-b,x) */

	if b_int && b >= 2 && !(a_int && a <= (b-2)) {
		/* Need to handle second term as a limit due to pole in
		   M(1+a-b,2-b,x).

		   We skip this case when a is integer and a <= b-2 because the
		   hypergeometric series terminates before any singular terms
		*/

		/* FIXME: TO BE IMPLEMENTED ! */
		result.val = math.NaN()
		result.err = math.NaN()
		return err.ERROR("limit case integer b >= 2 unimplemented", err.EUNIMPL)
	} else {
		if a_int && a <= 0 && (b >= 1) {
			r2.val = 0
			r2.err = 0
		} else {
			stat_2 = Poch_e(a, -(1 + a - b), &r2)
			status = err.ErrorSelect(status, stat_2)
		}

		if r2.val != 0.0 {
			var Mr2 Result
			stat_Mr2 := Hyperg_1F1_e(1+a-b, 2-b, x, &Mr2)
			T2 = Mr2.val * r2.val
			T2_err = 2.0*gsl.Float64Eps*math.Abs(T2) + math.Abs(Mr2.err*r2.val) + math.Abs(Mr2.val*r2.err)
			status = err.ErrorSelect(status, stat_Mr2)

			if T2 != 0.0 {
				x1mb := math.Pow(x, 1-b)
				T2 = x1mb * T2
				T2_err = math.Abs(x1mb) * T2_err
			}
		}
	}

	result.val = (T1 + T2)
	result.err = 2.0*gsl.Float64Eps*math.Abs(result.val) + (T1_err + T2_err)
	result.e10 = 0

	return status
}

func hyperg_U_int_negx(a, b int, x float64, result *Result_e10) err.GSLError {
	/* Looking at the tests it seems that everything is handled correctly by hyperg_U_negx
	   except a<b<=0.  The b=1 case seems strange since the poch(1+a-b,-a) should blow up.
	   In order to do no (little) harm I fix up only a<b<=0 using DLMF 13.2.11.  These are the failing conditions
	*/
	if a < b && b <= 0 {
		var r1 Result_e10
		z21_z := math.Pow(x, 1-float64(b))
		status := hyperg_U_negx(float64(1+a-b), float64(2-b), x, &r1)
		res_tem := z21_z * r1.val
		res_tem_err := math.Abs(z21_z) * r1.err
		result.val = res_tem
		result.err = res_tem_err
		return status
	} else {
		return hyperg_U_negx(float64(a), float64(b), x, result)
	}
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Hyperg_U_int_e10_e(a, b int, x float64, result *Result_e10) err.GSLError {

	if x == 0.0 && b >= 1 {
		return DomainError_e10(result)
	} else if x == 0.0 {

		return hyperg_U_int_origin(a, b, result)
	} else if x < 0.0 {
		return hyperg_U_int_negx(a, b, x, result)
	} else {
		if b >= 1 {
			return hyperg_U_int_bge1(a, b, x, result)
		} else {
			/* Use the reflection formula
			 * U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)
			 */

			var U Result_e10
			ln_x := math.Log(x)
			ap := 1 + a - b
			bp := 2 - b
			stat_U := hyperg_U_int_bge1(ap, bp, x, &U)
			ln_pre_val := (1.0 - float64(b)) * ln_x
			ln_pre_err := 2.0 * gsl.Float64Eps * (math.Abs(float64(b)) + 1.0) * math.Abs(ln_x)
			ln_pre_err += 2.0 * gsl.Float64Eps * math.Abs(1.0-float64(b)) /* error in math.Log(x) */
			stat_e := Exp_mult_err_e10_e(ln_pre_val+float64(U.e10)*gsl.Ln10, ln_pre_err, U.val, U.err, result)
			return err.ErrorSelect(stat_e, stat_U)
		}
	}
}

func Hyperg_U_e10_e(a, b, x float64, result *Result_e10) err.GSLError {
	rinta := math.Floor(a + 0.5)
	rintb := math.Floor(b + 0.5)
	a_integer := (math.Abs(a-rinta) < INT_THRESHOLD)
	b_integer := (math.Abs(b-rintb) < INT_THRESHOLD)

	/* CHECK_POINTER(result) */

	if x == 0.0 && b >= 1 {
		return DomainError_e10(result)
	} else if a == 0.0 {
		result.val = 1.0
		result.err = 0.0
		result.e10 = 0
		return nil
	} else if x == 0.0 {
		return hyperg_U_origin(a, b, result)
	} else if a_integer && b == a+1 /* This is DLMF 13.6.4 */ {
		var powx1N_1 Result
		Pow_int_e(x, -int(a), &powx1N_1)
		result.val = powx1N_1.val
		result.err = powx1N_1.err
		result.e10 = 0
		return nil

	} else if a_integer && b_integer {
		return Hyperg_U_int_e10_e(int(rinta), int(rintb), x, result)
	} else if x < 0.0 {
		return hyperg_U_negx(a, b, x, result)
	} else {
		if b >= 1.0 {
			/* Use b >= 1 function.
			 */
			return hyperg_U_bge1(a, b, x, result)
		} else {
			/* Use the reflection formula
			 * U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)
			 */
			lnx := math.Log(x)
			ln_pre_val := (1.0 - b) * lnx
			ln_pre_err := math.Abs(lnx) * 2.0 * gsl.Float64Eps * (1.0 + math.Abs(b))
			ap := 1.0 + a - b
			bp := 2.0 - b
			var U Result_e10
			stat_U := hyperg_U_bge1(ap, bp, x, &U)
			stat_e := Exp_mult_err_e10_e(ln_pre_val+float64(U.e10)*gsl.Ln10, ln_pre_err, U.val, U.err, result)
			return err.ErrorSelect(stat_e, stat_U)
		}
	}
}

func Hyperg_U_int_e(a, b int, x float64, result *Result) err.GSLError {
	re := new(Result_e10)
	stat_U := Hyperg_U_int_e10_e(a, b, x, re)
	stat_c := Result_smash_e(re, result)
	return err.ErrorSelect(stat_c, stat_U)
}

func Hyperg_U_e(a, b, x float64, result *Result) err.GSLError {
	re := new(Result_e10)
	stat_U := Hyperg_U_e10_e(a, b, x, re)
	stat_c := Result_smash_e(re, result)
	return err.ErrorSelect(stat_c, stat_U)
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Hyperg_U_int(a, b int, x float64) float64 {
	result := new(Result)
	status := Hyperg_U_int_e(a, b, x, result)
	return EvalResult(result, status)
}

func Hyperg_U(a, b, x float64) float64 {
	result := new(Result)
	status := Hyperg_U_e(a, b, x, result)
	return EvalResult(result, status)
}
