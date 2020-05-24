/* specfunc/bessel_y.c
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

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* [Abramowitz+Stegun, 10.1.3]
 * with lmax=15, precision ~ 15D for x < 3
 *
 * checked OK [GJ] Wed May 13 15:41:25 MDT 1998
 */
func bessel_yl_small_x(l int, x float64, result *Result) err.GSLError {
	var num_fact Result
	den := Pow_int(x, l+1)
	stat_df := Doublefact_e(uint(2*l-1), &num_fact)

	if stat_df != nil || den == 0.0 {
		return OverflowError(result)
	} else {
		lmax := 200
		t := -0.5 * x * x
		sum := 1.0
		t_coeff := 1.0
		t_power := 1.0
		var delta float64
		var i int
		for i = 1; i <= lmax; i++ {
			t_coeff /= float64(i * (2*(i-l) - 1))
			t_power *= t
			delta = t_power * t_coeff
			sum += delta
			if math.Abs(delta/sum) < 0.5*gsl.Float64Eps {
				break
			}
		}
		result.val = -num_fact.val / den * sum
		result.err = gsl.Float64Eps * math.Abs(result.val)

		return nil
	}
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_y0_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if 1.0/gsl.MaxFloat64 > 0.0 && x < 1.0/gsl.MaxFloat64 {
		return OverflowError(result)
	} else {
		var cos_result Result
		stat := Cos_e(x, &cos_result)
		result.val = -cos_result.val / x
		result.err = math.Abs(cos_result.err / x)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat
	}
}

func Bessel_y1_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if x < 1.0/gsl.SqrtMaxFloat64 {
		return OverflowError(result)
	} else if x < 0.25 {
		var (
			y   = x * x
			c1  = 1.0 / 2.0
			c2  = -1.0 / 8.0
			c3  = 1.0 / 144.0
			c4  = -1.0 / 5760.0
			c5  = 1.0 / 403200.0
			c6  = -1.0 / 43545600.0
			sum = 1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*c6)))))
		)
		result.val = -sum / y
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		var cos_result, sin_result Result
		stat_cos := Cos_e(x, &cos_result)
		stat_sin := Sin_e(x, &sin_result)
		cx := cos_result.val
		sx := sin_result.val
		result.val = -(cx/x + sx) / x
		result.err = (math.Abs(cos_result.err/x) + sin_result.err) / math.Abs(x)
		result.err += gsl.Float64Eps * (math.Abs(sx/x) + math.Abs(cx/(x*x)))
		return err.ErrorSelect(stat_cos, stat_sin)
	}
}

func Bessel_y2_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if x < 1.0/gsl.Root3MaxFloat64 {
		return OverflowError(result)
	} else if x < 0.5 {
		var (
			y   = x * x
			c1  = 1.0 / 6.0
			c2  = 1.0 / 24.0
			c3  = -1.0 / 144.0
			c4  = 1.0 / 3456.0
			c5  = -1.0 / 172800.0
			c6  = 1.0 / 14515200.0
			c7  = -1.0 / 1828915200.0
			sum = 1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))
		)
		result.val = -3.0 / (x * x * x) * sum
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		var cos_result, sin_result Result
		stat_cos := Cos_e(x, &cos_result)
		stat_sin := Sin_e(x, &sin_result)
		sx := sin_result.val
		cx := cos_result.val
		a := 3.0 / (x * x)
		result.val = (1.0-a)/x*cx - a*sx
		result.err = cos_result.err*math.Abs((1.0-a)/x) + sin_result.err*math.Abs(a)
		result.err += gsl.Float64Eps * (math.Abs(cx/x) + math.Abs(sx/(x*x)))
		return err.ErrorSelect(stat_cos, stat_sin)
	}
}

func Bessel_yl_e(l int, x float64, result *Result) err.GSLError {
	fl := float64(l)
	if l < 0 || x <= 0.0 {
		return DomainError(result)
	} else if l == 0 {
		return Bessel_y0_e(x, result)
	} else if l == 1 {
		return Bessel_y1_e(x, result)
	} else if l == 2 {
		return Bessel_y2_e(x, result)
	} else if x < 3.0 {
		return bessel_yl_small_x(l, x, result)
	} else if gsl.Root3Float64Eps*x > (fl*fl + fl + 1.0) {
		status := bessel_Ynu_asympx_e(fl+0.5, x, result)
		pre := math.Sqrt((0.5 * gsl.Pi) / x)
		result.val *= pre
		result.err *= pre
		return status
	} else if l > 40 {
		status := bessel_Ynu_asymp_Olver_e(fl+0.5, x, result)
		pre := math.Sqrt((0.5 * gsl.Pi) / x)
		result.val *= pre
		result.err *= pre
		return status
	} else {
		/* recurse upward */
		var r_by, r_bym Result
		stat_1 := Bessel_y1_e(x, &r_by)
		stat_0 := Bessel_y0_e(x, &r_bym)
		bym := r_bym.val
		by := r_by.val
		var byp float64
		var j int
		for j = 1; j < l; j++ {
			byp = float64(2*j+1)/x*by - bym
			bym = by
			by = byp
		}
		result.val = by
		result.err = math.Abs(result.val) * (gsl.Float64Eps + math.Abs(r_by.err/r_by.val) + math.Abs(r_bym.err/r_bym.val))

		return err.ErrorSelect(stat_1, stat_0)
	}
}

func Bessel_yl_array(lmax int, x float64, result_array []float64) err.GSLError {

	if lmax < 0 || x <= 0.0 {
		return err.ERROR("error", err.EDOM)
	} else if lmax == 0 {
		var result Result
		stat := Bessel_y0_e(x, &result)
		result_array[0] = result.val
		return stat
	} else {
		var r_yell, r_yellm1 Result
		stat_1 := Bessel_y1_e(x, &r_yell)
		stat_0 := Bessel_y0_e(x, &r_yellm1)
		var yellp1 float64
		yell := r_yell.val
		yellm1 := r_yellm1.val
		var ell int

		result_array[0] = yellm1
		result_array[1] = yell

		for ell = 1; ell < lmax; ell++ {
			yellp1 = (2*float64(ell)+1)/x*yell - yellm1
			result_array[ell+1] = yellp1
			yellm1 = yell
			yell = yellp1
		}

		return err.ErrorSelect(stat_0, stat_1)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_y0(x float64) float64 {
	result := new(Result)
	status := Bessel_y0_e(x, result)
	return EvalResult(result, status)
}

func Bessel_y1(x float64) float64 {
	result := new(Result)
	status := Bessel_y1_e(x, result)
	return EvalResult(result, status)
}

func Bessel_y2(x float64) float64 {
	result := new(Result)
	status := Bessel_y2_e(x, result)
	return EvalResult(result, status)
}

func Bessel_yl(l int, x float64) float64 {
	result := new(Result)
	status := Bessel_yl_e(l, x, result)
	return EvalResult(result, status)
}
