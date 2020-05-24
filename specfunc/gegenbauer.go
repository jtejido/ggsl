/* specfunc/gegenbauer.c
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

/* See: [Thompson, Atlas for Computing Mathematical Functions] */
func Gegenpoly_1_e(lambda, x float64, result *Result) err.GSLError {
	if lambda == 0.0 {
		result.val = 2.0 * x
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = 2.0 * lambda * x
		result.err = 4.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Gegenpoly_2_e(lambda, x float64, result *Result) err.GSLError {
	if lambda == 0.0 {
		txx := 2.0 * x * x
		result.val = -1.0 + txx
		result.err = 2.0 * gsl.Float64Eps * math.Abs(txx)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = lambda * (-1.0 + 2.0*(1.0+lambda)*x*x)
		result.err = gsl.Float64Eps * (2.0*math.Abs(result.val) + math.Abs(lambda))
		return nil
	}
}

func Gegenpoly_3_e(lambda, x float64, result *Result) err.GSLError {

	if lambda == 0.0 {
		result.val = x * (-2.0 + 4.0/3.0*x*x)
		result.err = gsl.Float64Eps * (2.0*math.Abs(result.val) + math.Abs(x))
		return nil
	} else {
		c := 4.0 + lambda*(6.0+2.0*lambda)
		result.val = 2.0 * lambda * x * (-1.0 - lambda + c*x*x/3.0)
		result.err = gsl.Float64Eps * (2.0*math.Abs(result.val) + math.Abs(lambda*x))
		return nil
	}
}

func Gegenpoly_n_e(n int, lambda, x float64, result *Result) err.GSLError {

	if lambda <= -0.5 || n < 0 {
		return DomainError(result)
	} else if n == 0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if n == 1 {
		return Gegenpoly_1_e(lambda, x, result)
	} else if n == 2 {
		return Gegenpoly_2_e(lambda, x, result)
	} else if n == 3 {
		return Gegenpoly_3_e(lambda, x, result)
	} else {
		if lambda == 0.0 && (x >= -1.0 && x <= 1.0) {
			/* 2 T_n(x)/n */
			z := float64(n) * math.Acos(x)
			result.val = 2.0 * math.Cos(z) / float64(n)
			result.err = 2.0 * gsl.Float64Eps * math.Abs(z*result.val)
			return nil
		} else {
			var k int
			var g2, g3 Result
			stat_g2 := Gegenpoly_2_e(lambda, x, &g2)
			stat_g3 := Gegenpoly_3_e(lambda, x, &g3)
			stat_g := err.ErrorSelect(stat_g2, stat_g3)
			gkm2 := g2.val
			gkm1 := g3.val
			gk := 0.0
			for k = 4; k <= n; k++ {
				fk := float64(k)
				gk = (2.0*(fk+lambda-1.0)*x*gkm1 - (fk+2.0*lambda-2.0)*gkm2) / fk
				gkm2 = gkm1
				gkm1 = gk
			}
			result.val = gk
			result.err = 2.0 * gsl.Float64Eps * 0.5 * float64(n) * math.Abs(gk)
			return stat_g
		}
	}
}

func Gegenpoly_array(nmax int, lambda, x float64, result_array []float64) err.GSLError {
	var k int
	if lambda <= -0.5 || nmax < 0 {
		return err.ERROR("domain error", err.EDOM)
	}

	/* n == 0 */
	result_array[0] = 1.0
	if nmax == 0 {
		return nil
	}

	/* n == 1 */
	if lambda == 0.0 {
		result_array[1] = 2.0 * x
	} else {
		result_array[1] = 2.0 * lambda * x
	}

	/* n <= nmax */
	for k = 2; k <= nmax; k++ {
		fk := float64(k)
		term1 := 2.0 * (fk + lambda - 1.0) * x * result_array[k-1]
		term2 := (fk + 2.0*lambda - 2.0) * result_array[k-2]
		result_array[k] = (term1 - term2) / fk
	}

	return nil
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Gegenpoly_1(lambda, x float64) float64 {
	result := new(Result)
	status := Gegenpoly_1_e(lambda, x, result)
	return EvalResult(result, status)
}

func Gegenpoly_2(lambda, x float64) float64 {
	result := new(Result)
	status := Gegenpoly_2_e(lambda, x, result)
	return EvalResult(result, status)
}

func Gegenpoly_3(lambda, x float64) float64 {
	result := new(Result)
	status := Gegenpoly_3_e(lambda, x, result)
	return EvalResult(result, status)
}

func Gegenpoly_n(n int, lambda, x float64) float64 {
	result := new(Result)
	status := Gegenpoly_n_e(n, lambda, x, result)
	return EvalResult(result, status)
}
