/* specfunc/beta_inc.c
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

func beta_cont_frac(a, b, x float64, result *Result) err.GSLError {
	max_iter := 512                /* control iterations      */
	cutoff := 2.0 * gsl.MinFloat64 /* control the zero cutoff */
	iter_count := 0
	var cf float64

	/* standard initialization for continued fraction */
	num_term := 1.0
	den_term := 1.0 - (a+b)*x/(a+1.0)
	if math.Abs(den_term) < cutoff {
		den_term = cutoff
	}
	den_term = 1.0 / den_term
	cf = den_term

	for iter_count < max_iter {
		k := iter_count + 1
		fk := float64(k)
		coeff := fk * (b - fk) * x / (((a - 1.0) + 2*fk) * (a + 2*fk))
		var delta_frac float64

		/* first step */
		den_term = 1.0 + coeff*den_term
		num_term = 1.0 + coeff/num_term
		if math.Abs(den_term) < cutoff {
			den_term = cutoff
		}
		if math.Abs(num_term) < cutoff {
			num_term = cutoff
		}
		den_term = 1.0 / den_term

		delta_frac = den_term * num_term
		cf *= delta_frac

		coeff = -(a + fk) * (a + b + fk) * x / ((a + 2*fk) * (a + 2*fk + 1.0))

		/* second step */
		den_term = 1.0 + coeff*den_term
		num_term = 1.0 + coeff/num_term
		if math.Abs(den_term) < cutoff {
			den_term = cutoff
		}
		if math.Abs(num_term) < cutoff {
			num_term = cutoff
		}
		den_term = 1.0 / den_term

		delta_frac = den_term * num_term
		cf *= delta_frac

		if math.Abs(delta_frac-1.0) < 2.0*gsl.Float64Eps {
			break
		}

		iter_count++
	}

	result.val = cf
	result.err = float64(iter_count) * 4.0 * gsl.Float64Eps * math.Abs(cf)

	if iter_count >= max_iter {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Beta_inc_e(a, b, x float64, result *Result) err.GSLError {
	if x < 0.0 || x > 1.0 {
		return DomainError(result)
	} else if isnegint(a) || isnegint(b) {
		return DomainError(result)
	} else if isnegint(a + b) {
		return DomainError(result)
	} else if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if x == 1.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if a <= 0 || b <= 0 {
		var f, beta Result
		var stat err.GSLError
		stat_f := Hyperg_2F1_e(a, 1-b, a+1, x, &f)
		stat_beta := Beta_e(a, b, &beta)
		prefactor := (math.Pow(x, a) / a)
		result.val = prefactor * f.val / beta.val
		result.err = math.Abs(prefactor)*f.err/math.Abs(beta.val) + math.Abs(result.val/beta.val)*beta.err

		stat = err.ErrorSelect(stat_f, stat_beta)
		if stat == nil {
			CheckUnderflow(result)
		}

		return stat
	} else {
		var ln_beta, ln_x, ln_1mx, prefactor Result
		stat_ln_beta := Lnbeta_e(a, b, &ln_beta)
		stat_ln_1mx := Log_1plusx_e(-x, &ln_1mx)
		stat_ln_x := Log_e(x, &ln_x)
		stat_ln := err.ErrorSelect(stat_ln_beta, stat_ln_1mx, stat_ln_x)

		ln_pre_val := -ln_beta.val + a*ln_x.val + b*ln_1mx.val
		ln_pre_err := ln_beta.err + math.Abs(a*ln_x.err) + math.Abs(b*ln_1mx.err)
		stat_exp := Exp_err_e(ln_pre_val, ln_pre_err, &prefactor)
		if stat_ln != nil {
			result.val = 0.0
			result.err = 0.0
			return err.ERROR("error", err.ESANITY)
		}

		if x < (a+1.0)/(a+b+2.0) {
			/* Apply continued fraction directly. */
			var cf Result
			stat_cf := beta_cont_frac(a, b, x, &cf)
			var stat err.GSLError
			result.val = prefactor.val * cf.val / a
			result.err = (math.Abs(prefactor.err*cf.val) + math.Abs(prefactor.val*cf.err)) / a

			stat = err.ErrorSelect(stat_exp, stat_cf)
			if stat == nil {
				CheckUnderflow(result)
			}

			return stat
		} else {
			/* Apply continued fraction after hypergeometric transformation. */
			var cf Result
			stat_cf := beta_cont_frac(b, a, 1.0-x, &cf)
			var stat err.GSLError
			term := prefactor.val * cf.val / b
			result.val = 1.0 - term
			result.err = math.Abs(prefactor.err*cf.val) / b
			result.err += math.Abs(prefactor.val*cf.err) / b
			result.err += 2.0 * gsl.Float64Eps * (1.0 + math.Abs(term))
			/* since the prefactor term is subtracted from 1 we need to
			   ignore underflow */
			if stat_exp != nil && stat_exp.Status() != err.EUNDRFLW {
				stat = err.ErrorSelect(stat_exp, stat_cf)
			} else {
				stat = stat_cf
			}

			if stat == nil {
				CheckUnderflow(result)
			}
			return stat
		}
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Beta_inc(a, b, x float64) float64 {
	result := new(Result)
	status := Beta_inc_e(a, b, x, result)
	return EvalResult(result, status)
}
