/* specfunc/sincos_pi.c
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

const TWOBIG float64 = 2.0 / gsl.Float64Eps

func sin_pi_taylor(x float64, result *Result) err.GSLError {
	result.val = 0.0
	result.err = 0.0
	if 16.0*math.Abs(x) < 1.0 {
		y := gsl.Pi * x
		a := y * y
		result.val = y * (1.0 - a*(1.0-a*(1.0-a*(1.0-a*(1.0-a/110.0)/72.0)/42.0)/20.0)/6.0)
	} else {
		result.val = math.Sin(gsl.Pi * x)
	}

	result.err = gsl.Float64Eps * math.Abs(result.val)
	return nil
}

func cos_pi_taylor(x float64, result *Result) err.GSLError {
	result.val = 0.0
	result.err = 0.0
	if 20.0*math.Abs(x) < 1.0 {
		y := gsl.Pi * x
		a := y * y
		result.val = 1.0 - 0.5*a*(1.0-a*(1.0-a*(1.0-a*(1.0-a/90.0)/56.0)/30.0)/12.0)
	} else {
		result.val = math.Cos(gsl.Pi * x)
	}

	result.err = gsl.Float64Eps * math.Abs(result.val)
	return nil
}

func Sin_pi_e(x float64, result *Result) err.GSLError {
	var q float64

	sign := 1
	result.val = 0.0
	result.err = 0.0
	intx, fracx := math.Modf(x)
	if fracx == 0.0 {
		return nil
	}

	if math.Abs(intx) >= TWOBIG {
		return nil
	}

	if (intx >= math.MinInt64) && (intx <= math.MaxInt64) {
		q = intx
	} else {
		q = math.Mod(intx, 2.)
	}

	if int64(q)%2 != 0 {
		sign = -1
	}

	if math.Abs(fracx) == 0.5 {
		if fracx < 0.0 {
			sign = -sign
		}

		if sign != 1 {
			result.val = -1.
		} else {
			result.val = 1.
		}

		return nil
	}

	if math.Abs(fracx) > 0.5 {
		sign = -sign

		if fracx > 0.0 {
			fracx = fracx - 1
		} else {
			fracx = fracx + 1
		}
	}
	var status err.GSLError
	if fracx > 0.25 {
		status = cos_pi_taylor((fracx - 0.5), result)
	} else if fracx < -0.25 {
		status = cos_pi_taylor((fracx + 0.5), result)
		sign = -sign
	} else {
		status = sin_pi_taylor(fracx, result)
	}

	if sign != 1 {
		result.val = -result.val
	}

	return status
}

func Cos_pi_e(x float64, result *Result) err.GSLError {
	var q float64
	sign := 1
	result.val = 0.0
	result.err = 0.0
	intx, fracx := math.Modf(x)
	if math.Abs(fracx) == 0.5 {
		return nil
	}

	if math.Abs(intx) >= TWOBIG {
		result.val = 1.
		return nil
	}
	if (intx >= math.MinInt64) && (intx <= math.MaxInt64) {
		q = intx
	} else {
		q = math.Mod(intx, 2.)
	}

	if int64(q)%2 != 0 {
		sign = -1
	}

	if fracx == 0.0 {
		if sign != 1 {
			result.val = -1
		} else {
			result.val = 1
		}
		return nil
	}
	if math.Abs(fracx) > 0.5 {
		sign = -sign

		if fracx > 0.0 {
			fracx -= 1
		} else {
			fracx += 1
		}
	}

	var status err.GSLError
	if fracx > 0.25 {
		status = sin_pi_taylor((fracx - 0.5), result)
		sign = -sign
	} else if fracx < -0.25 {
		status = sin_pi_taylor((fracx + 0.5), result)
	} else {
		status = cos_pi_taylor(fracx, result)
	}

	if sign != 1 {
		result.val = -result.val
	}

	return status
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Sin_pi(x float64) float64 {
	result := new(Result)
	status := Sin_pi_e(x, result)
	return EvalResult(result, status)
}

func Cos_pi(x float64) float64 {
	result := new(Result)
	status := Cos_pi_e(x, result)
	return EvalResult(result, status)
}
