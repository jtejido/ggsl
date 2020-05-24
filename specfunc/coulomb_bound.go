/* specfunc/coulomb_bound.c
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
 */package specfunc

import (
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"math"
)

/* normalization for hydrogenic wave functions */
func rNorm(n, l int, Z float64, result *Result) err.GSLError {
	A := 2.0 * Z / float64(n)
	pre := math.Sqrt(A * A * A / (2.0 * float64(n)))
	var ln_a, ln_b, ex Result
	stat_a := Lnfact_e(uint(n+l), &ln_a)
	stat_b := Lnfact_e(uint(n-l-1), &ln_b)
	diff_val := 0.5 * (ln_b.val - ln_a.val)
	diff_err := 0.5*(ln_b.err+ln_a.err) + gsl.Float64Eps*math.Abs(diff_val)
	stat_e := Exp_err_e(diff_val, diff_err, &ex)
	result.val = pre * ex.val
	result.err = pre * ex.err
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return err.ErrorSelect(stat_e, stat_a, stat_b)
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func HydrogenicR_1_e(Z, r float64, result *Result) err.GSLError {
	if Z > 0.0 && r >= 0.0 {
		A := 2.0 * Z
		norm := A * math.Sqrt(Z)
		ea := math.Exp(-Z * r)
		result.val = norm * ea
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val) * math.Abs(Z*r)
		CheckUnderflow(result)
		return nil
	} else {
		return DomainError(result)
	}
}

func HydrogenicR_e(n, l int, Z, r float64, result *Result) err.GSLError {
	if n < 1 || l > n-1 || Z <= 0.0 || r < 0.0 {
		return DomainError(result)
	} else {
		A := 2.0 * Z / float64(n)
		var norm Result
		stat_norm := rNorm(n, l, Z, &norm)
		rho := A * r
		ea := math.Exp(-0.5 * rho)
		pp := Pow_int(rho, l)
		var lag Result
		stat_lag := Laguerre_n_e(n-l-1, float64(2*l+1), rho, &lag)
		W_val := norm.val * ea * pp
		W_err := norm.err * ea * pp
		W_err += norm.val * ((0.5*rho + 1.0) * gsl.Float64Eps) * ea * pp
		W_err += norm.val * ea * ((float64(l) + 1.0) * gsl.Float64Eps) * pp
		result.val = W_val * lag.val
		result.err = W_val*lag.err + W_err*math.Abs(lag.val)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		if (l == 0 || (r > 0 && l > 0)) && lag.val != 0.0 && stat_lag == nil && stat_norm == nil {
			CheckUnderflow(result)
		}

		return err.ErrorSelect(stat_lag, stat_norm)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func HydrogenicR_1(Z, r float64) float64 {
	result := new(Result)
	status := HydrogenicR_1_e(Z, r, result)
	return EvalResult(result, status)
}

func HydrogenicR(n, l int, Z, r float64) float64 {
	result := new(Result)
	status := HydrogenicR_e(n, l, Z, r, result)
	return EvalResult(result, status)
}
