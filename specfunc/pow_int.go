/* specfunc/pow_int.c
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

/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error handling *-*-*-*-*-*-*-*-*-*-*-*/

func Pow_int_e(x float64, n int, result *Result) err.GSLError {
	value := 1.0
	count := 0

	if n < 0 {
		n = -n

		if x == 0.0 {
			u := 1.0 / x
			if (n % 2) != 0 { /* correct sign of infinity */
				result.val = u
			} else {
				result.val = (u * u)
			}

			result.err = math.Inf(1)
			return err.ERROR("overflow", err.EOVRFLW)
		}

		x = 1.0 / x
	}

	/* repeated squaring method
	 * returns 0.0^0 = 1.0, so continuous in x
	 */
	for ok := true; ok; ok = n != 0 {
		if gsl.IsOdd(n) {
			value *= x
		}
		n >>= 1
		x *= x
		count++
	}

	result.val = value
	result.err = 2.0 * gsl.Float64Eps * (float64(count) + 1.0) * math.Abs(value)

	return nil
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Pow_int(x float64, n int) float64 {
	result := new(Result)
	status := Pow_int_e(x, n, result)
	return EvalResult(result, status)
}
