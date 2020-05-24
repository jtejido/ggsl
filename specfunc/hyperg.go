/* specfunc/hyperg.c
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

var (
	SUM_LARGE = 1.0e-5 * gsl.MaxFloat64
)

func Hyperg_1F1_series_e(a, b, x float64, result *Result) err.GSLError {
	var (
		an          = a
		bn          = b
		n           = 1.0
		del         = 1.0
		abs_del     = 1.0
		max_abs_del = 1.0
		sum_val     = 1.0
		sum_err     = 0.0
	)

	for abs_del/math.Abs(sum_val) > 0.25*gsl.Float64Eps {
		var u, abs_u float64

		if bn == 0.0 {
			return DomainError(result)
		}

		if an == 0.0 {
			result.val = sum_val
			result.err = sum_err
			result.err += 2.0 * gsl.Float64Eps * n * math.Abs(sum_val)
			return nil
		}

		if n > 10000.0 {
			result.val = sum_val
			result.err = sum_err
			return err.ERROR("hypergeometric series failed to converge", err.EFAILED)
		}

		u = x * (an / (bn * n))
		abs_u = math.Abs(u)
		if abs_u > 1.0 && max_abs_del > gsl.MaxFloat64/abs_u {
			result.val = sum_val
			result.err = math.Abs(sum_val)
			return err.ERROR("overflow", err.EOVRFLW)
		}
		del *= u
		sum_val += del
		if math.Abs(sum_val) > SUM_LARGE {
			result.val = sum_val
			result.err = math.Abs(sum_val)
			return err.ERROR("overflow", err.EOVRFLW)
		}

		abs_del = math.Abs(del)
		max_abs_del = math.Max(abs_del, max_abs_del)
		sum_err += 2.0 * gsl.Float64Eps * abs_del

		an += 1.0
		bn += 1.0
		n += 1.0
	}

	result.val = sum_val
	result.err = sum_err
	result.err += abs_del
	result.err += 2.0 * gsl.Float64Eps * n * math.Abs(sum_val)

	return nil
}

func Hyperg_1F1_large_b_e(a, b, x float64, result *Result) err.GSLError {
	if math.Abs(x/b) < 1.0 {
		var (
			u   = x / b
			v   = 1.0 / (1.0 - u)
			pre = math.Pow(v, a)
			uv  = u * v
			uv2 = uv * uv
			t1  = a * (a + 1.0) / (2.0 * b) * uv2
			t2a = a * (a + 1.0) / (24.0 * b * b) * uv2
			t2b = 12.0 + 16.0*(a+2.0)*uv + 3.0*(a+2.0)*(a+3.0)*uv2
			t2  = t2a * t2b
		)
		result.val = pre * (1.0 - t1 + t2)
		result.err = pre * gsl.Float64Eps * (1.0 + math.Abs(t1) + math.Abs(t2))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		return DomainError(result)
	}
}

func Hyperg_U_large_b_e(a, b, x float64, result *Result, ln_multiplier *float64) err.GSLError {
	N := math.Floor(b) /* b = N + eps */
	eps := b - N

	if math.Abs(eps) < gsl.SqrtFloat64Eps {
		var lnpre_val, lnpre_err float64
		M := new(Result)
		if b > 1.0 {
			tmp := (1.0 - b) * math.Log(x)
			lg_bm1, lg_a := new(Result), new(Result)
			Lngamma_e(b-1.0, lg_bm1)
			Lngamma_e(a, lg_a)
			lnpre_val = tmp + x + lg_bm1.val - lg_a.val
			lnpre_err = lg_bm1.err + lg_a.err + gsl.Float64Eps*(math.Abs(x)+math.Abs(tmp))
			Hyperg_1F1_large_b_e(1.0-a, 2.0-b, -x, M)
		} else {
			lg_1mb, lg_1pamb := new(Result), new(Result)
			Lngamma_e(1.0-b, lg_1mb)
			Lngamma_e(1.0+a-b, lg_1pamb)
			lnpre_val = lg_1mb.val - lg_1pamb.val
			lnpre_err = lg_1mb.err + lg_1pamb.err
			Hyperg_1F1_large_b_e(a, b, x, M)
		}

		if lnpre_val > gsl.LnMaxFloat64-10.0 {
			result.val = M.val
			result.err = M.err
			*ln_multiplier = lnpre_val
			return err.ERROR("overflow", err.EOVRFLW)
		} else {
			epre := new(Result)
			stat_e := Exp_err_e(lnpre_val, lnpre_err, epre)
			result.val = epre.val * M.val
			result.err = epre.val*M.err + epre.err*math.Abs(M.val)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			*ln_multiplier = 0.0
			return stat_e
		}
	} else {
		omb_lnx := (1.0 - b) * math.Log(x)
		var (
			sgn_1mb, sgn_1pamb, sgn_bm1, sgn_a, lnpre1_val, lnpre2_val, lnpre1_err, lnpre2_err, sgpre1, sgpre2 float64
			lg_1mb, lg_1pamb, lg_bm1, lg_a, M1, M2                                                             Result
		)

		Hyperg_1F1_large_b_e(a, b, x, &M1)
		Hyperg_1F1_large_b_e(1.0-a, 2.0-b, x, &M2)

		Lngamma_sgn_e(1.0-b, &lg_1mb, &sgn_1mb)
		Lngamma_sgn_e(1.0+a-b, &lg_1pamb, &sgn_1pamb)

		Lngamma_sgn_e(b-1.0, &lg_bm1, &sgn_bm1)
		Lngamma_sgn_e(a, &lg_a, &sgn_a)

		lnpre1_val = lg_1mb.val - lg_1pamb.val
		lnpre1_err = lg_1mb.err + lg_1pamb.err
		lnpre2_val = lg_bm1.val - lg_a.val - omb_lnx - x
		lnpre2_err = lg_bm1.err + lg_a.err + gsl.Float64Eps*(math.Abs(omb_lnx)+math.Abs(x))
		sgpre1 = sgn_1mb * sgn_1pamb
		sgpre2 = sgn_bm1 * sgn_a

		if lnpre1_val > gsl.LnMaxFloat64-10.0 || lnpre2_val > gsl.LnMaxFloat64-10.0 {
			max_lnpre_val := math.Max(lnpre1_val, lnpre2_val)
			max_lnpre_err := math.Max(lnpre1_err, lnpre2_err)
			lp1 := lnpre1_val - max_lnpre_val
			lp2 := lnpre2_val - max_lnpre_val
			t1 := sgpre1 * math.Exp(lp1)
			t2 := sgpre2 * math.Exp(lp2)
			result.val = t1*M1.val + t2*M2.val
			result.err = math.Abs(t1)*M1.err + math.Abs(t2)*M2.err
			result.err += gsl.Float64Eps * math.Exp(max_lnpre_err) * (math.Abs(t1*M1.val) + math.Abs(t2*M2.val))
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			*ln_multiplier = max_lnpre_val
			return err.ERROR("overflow", err.EOVRFLW)
		} else {
			t1 := sgpre1 * math.Exp(lnpre1_val)
			t2 := sgpre2 * math.Exp(lnpre2_val)
			result.val = t1*M1.val + t2*M2.val
			result.err = math.Abs(t1)*M1.err + math.Abs(t2)*M2.err
			result.err += gsl.Float64Eps * (math.Exp(lnpre1_err)*math.Abs(t1*M1.val) + math.Exp(lnpre2_err)*math.Abs(t2*M2.val))
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			*ln_multiplier = 0.0
			return nil
		}
	}
}

/* [Carlson, p.109] says the error in truncating this asymptotic series
 * is less than the absolute value of the first neglected term.
 *
 * A termination argument is provided, so that the series will
 * be summed at most up to n=n_trunc. If n_trunc is set negative,
 * then the series is summed until it appears to start diverging.
 */
func Hyperg_2F0_series_e(a, b, x float64, n_trunc int, result *Result) err.GSLError {
	var (
		maxiter      = 2000
		an           = a
		bn           = b
		n            = 1.0
		sum          = 1.0
		del          = 1.0
		abs_del      = 1.0
		max_abs_del  = 1.0
		last_abs_del = 1.0
	)

	for abs_del/math.Abs(sum) > gsl.Float64Eps && n < float64(maxiter) {

		u := an * (bn / n * x)
		abs_u := math.Abs(u)

		if abs_u > 1.0 && (max_abs_del > gsl.MaxFloat64/abs_u) {
			result.val = sum
			result.err = math.Abs(sum)
			return err.ERROR("overflow", err.EOVRFLW)
		}

		del *= u
		sum += del

		abs_del = math.Abs(del)

		if abs_del > last_abs_del {
			break /* series is probably starting to grow */
		}

		last_abs_del = abs_del
		max_abs_del = math.Max(abs_del, max_abs_del)

		an += 1.0
		bn += 1.0
		n += 1.0

		if an == 0.0 || bn == 0.0 {
			break /* series terminated */
		}

		if n_trunc >= 0 && n >= float64(n_trunc) {
			break /* reached requested timeout */
		}
	}

	result.val = sum
	result.err = gsl.Float64Eps*n + abs_del
	if n >= float64(maxiter) {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}
