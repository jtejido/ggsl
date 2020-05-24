/* specfunc/beta.c
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

func isnegint(x float64) bool {
	return (x < 0) && (x == math.Floor(x))
}

func Lnbeta_e(x, y float64, result *Result) err.GSLError {
	var sgn float64
	status := Lnbeta_sgn_e(x, y, result, &sgn)
	if sgn == -1 {
		return DomainError(result)
	}

	return status
}

func Lnbeta_sgn_e(x, y float64, result *Result, sgn *float64) err.GSLError {
	if x == 0.0 || y == 0.0 {
		*sgn = 0.0
		return DomainError(result)
	} else if isnegint(x) || isnegint(y) {
		*sgn = 0.0
		return DomainError(result)
	}

	/* See if we can handle the postive case with min/max < 0.2 */
	if x > 0 && y > 0 {
		max := math.Max(x, y)
		min := math.Min(x, y)
		rat := min / max

		if rat < 0.2 {
			var (
				lnpre_val  float64
				lnpre_err  float64
				lnpow_val  float64
				lnpow_err  float64
				t1, t2, t3 float64
			)
			lnopr := new(Result)
			gsx, gsy, gsxy := new(Result), new(Result), new(Result)
			Gammastar_e(x, gsx)
			Gammastar_e(y, gsy)
			Gammastar_e(x+y, gsxy)
			Log_1plusx_e(rat, lnopr)
			lnpre_val = math.Log(gsx.val * gsy.val / gsxy.val * math.Sqrt2 * math.SqrtPi)
			lnpre_err = gsx.err/gsx.val + gsy.err/gsy.val + gsxy.err/gsxy.val
			t1 = min * math.Log(rat)
			t2 = 0.5 * math.Log(min)
			t3 = (x + y - 0.5) * lnopr.val
			lnpow_val = t1 - t2 - t3
			lnpow_err = gsl.Float64Eps * (math.Abs(t1) + math.Abs(t2) + math.Abs(t3))
			lnpow_err += math.Abs(x+y-0.5) * lnopr.err
			result.val = lnpre_val + lnpow_val
			result.err = lnpre_err + lnpow_err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			*sgn = 1.0
			return nil
		}
	}

	lgx, lgy, lgxy := new(Result), new(Result), new(Result)
	var sgx, sgy, sgxy float64
	xy := x + y
	stat_gx := Lngamma_sgn_e(x, lgx, &sgx)
	stat_gy := Lngamma_sgn_e(y, lgy, &sgy)
	stat_gxy := Lngamma_sgn_e(xy, lgxy, &sgxy)

	*sgn = sgx * sgy * sgxy
	result.val = lgx.val + lgy.val - lgxy.val
	result.err = lgx.err + lgy.err + lgxy.err
	result.err += 2.0 * gsl.Float64Eps * (math.Abs(lgx.val) + math.Abs(lgy.val) + math.Abs(lgxy.val))
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return err.ErrorSelect(stat_gx, stat_gy, stat_gxy)

}

func Beta_e(x, y float64, result *Result) err.GSLError {
	if (x > 0 && y > 0) && x < 50.0 && y < 50.0 {
		gx, gy, gxy := new(Result), new(Result), new(Result)
		Gamma_e(x, gx)
		Gamma_e(y, gy)
		Gamma_e(x+y, gxy)
		result.val = (gx.val * gy.val) / gxy.val
		result.err = gx.err * math.Abs(gy.val/gxy.val)
		result.err += gy.err * math.Abs(gx.val/gxy.val)
		result.err += math.Abs((gx.val*gy.val)/(gxy.val*gxy.val)) * gxy.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if isnegint(x) || isnegint(y) {
		return DomainError(result)
	} else if isnegint(x + y) {
		result.val = 0.0
		result.err = 0.0
		return nil
	}

	lb := new(Result)
	var sgn float64
	stat_lb := Lnbeta_sgn_e(x, y, lb, &sgn)
	if stat_lb == nil {
		status := Exp_err_e(lb.val, lb.err, result)
		result.val *= sgn
		return status
	}

	result.val = 0.0
	result.err = 0.0
	return stat_lb
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

func Lnbeta(x, y float64) float64 {
	result := new(Result)
	status := Lnbeta_e(x, y, result)
	return EvalResult(result, status)
}

func Beta(x, y float64) float64 {
	result := new(Result)
	status := Beta_e(x, y, result)
	return EvalResult(result, status)
}
