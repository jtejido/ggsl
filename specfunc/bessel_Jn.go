/* specfunc/bessel_Jn.c
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

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_Jn_e(n int, x float64, result *Result) err.GSLError {
	sign := 1

	if n < 0 {
		/* reduce to case n >= 0 */
		n = -n
		if gsl.IsOdd(n) {
			sign = -sign
		}
	}

	if x < 0.0 {
		/* reduce to case x >= 0. */
		x = -x
		if gsl.IsOdd(n) {
			sign = -sign
		}
	}

	if n == 0 {
		b0 := new(Result)
		stat_J0 := Bessel_J0_e(x, b0)
		result.val = float64(sign) * b0.val
		result.err = b0.err
		return stat_J0
	} else if n == 1 {
		b1 := new(Result)
		stat_J1 := Bessel_J1_e(x, b1)
		result.val = float64(sign) * b1.val
		result.err = b1.err
		return stat_J1
	} else {
		if x == 0.0 {
			result.val = 0.0
			result.err = 0.0
			return nil
		} else if x*x < 10.0*(float64(n)+1.0)*gsl.Root5Float64Eps {
			b := new(Result)
			status := Bessel_IJ_taylor_e(float64(n), x, -1, 50, gsl.Float64Eps, b)
			result.val = float64(sign) * b.val
			result.err = b.err
			result.err += gsl.Float64Eps * math.Abs(result.val)
			return status
		} else if gsl.Root4Float64Eps*x > (float64(n)*float64(n) + 1.0) {
			status := bessel_Jnu_asympx_e(float64(n), x, result)
			result.val *= float64(sign)
			return status
		} else if n > 50 {
			status := bessel_Jnu_asymp_Olver_e(float64(n), x, result)
			result.val *= float64(sign)
			return status
		} else if x > 1000.0 {
			/* We need this to avoid feeding large x to CF1; note that
			 * due to the above check, we know that n <= 50.
			 */
			status := bessel_Jnu_asympx_e(float64(n), x, result)
			result.val *= float64(sign)
			return status
		} else {
			var ans, errs, ratio, sgn float64
			var stat_b err.GSLError
			stat_CF1 := bessel_J_CF1(float64(n), x, &ratio, &sgn)

			/* backward recurrence */
			Jkp1 := gsl.SqrtMinFloat64 * ratio
			Jk := gsl.SqrtMinFloat64

			for k := n; k > 0; k-- {
				Jkm1 := 2.0*float64(k)/x*Jk - Jkp1
				Jkp1 = Jk
				Jk = Jkm1
			}

			if math.Abs(Jkp1) > math.Abs(Jk) {
				b1 := new(Result)
				stat_b = Bessel_J1_e(x, b1)
				ans = b1.val / Jkp1 * gsl.SqrtMinFloat64
				errs = b1.err / Jkp1 * gsl.SqrtMinFloat64
			} else {
				b0 := new(Result)
				stat_b = Bessel_J0_e(x, b0)
				ans = b0.val / Jk * gsl.SqrtMinFloat64
				errs = b0.err / Jk * gsl.SqrtMinFloat64
			}

			result.val = float64(sign) * ans
			result.err = math.Abs(errs)
			return err.ErrorSelect(stat_CF1, stat_b)
		}
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_Jn(n int, x float64) float64 {
	result := new(Result)
	status := Bessel_Jn_e(n, x, result)
	return EvalResult(result, status)
}
