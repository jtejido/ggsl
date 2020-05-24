/* specfunc/bessel_i.c
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

/* i_{l+1}/i_l
 */
func bessel_il_CF1(l int, x, threshold float64, ratio *float64) err.GSLError {
	kmax := 2000
	tk := 1.0
	sum := 1.0
	rhok := 0.0
	fl := float64(l)
	var k int

	for k = 1; k <= kmax; k++ {
		fk := float64(k)
		ak := (x / (2.0*fl + 1.0 + 2.0*fk)) * (x / (2.0*fl + 3.0 + 2.0*fk))
		rhok = -ak * (1.0 + rhok) / (1.0 + ak*(1.0+rhok))
		tk *= rhok
		sum += tk
		if math.Abs(tk/sum) < threshold {
			break
		}
	}

	*ratio = x / (2.0*fl + 3.0) * sum

	if k == kmax {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_i0_scaled_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	if x == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if ax < 0.2 {
		var (
			eax = math.Exp(-ax)
			y   = ax * ax
			c1  = 1.0 / 6.0
			c2  = 1.0 / 120.0
			c3  = 1.0 / 5040.0
			c4  = 1.0 / 362880.0
			c5  = 1.0 / 39916800.0
			sum = 1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*c5))))
		)
		result.val = eax * sum
		result.err = 2.0 * gsl.Float64Eps * result.val
	} else if ax < -0.5*gsl.LnFloat64Eps {
		result.val = (1.0 - math.Exp(-2.0*ax)) / (2.0 * ax)
		result.err = 2.0 * gsl.Float64Eps * result.val
	} else {
		result.val = 1.0 / (2.0 * ax)
		result.err = 2.0 * gsl.Float64Eps * result.val
	}
	return nil
}

func Bessel_i1_scaled_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if ax < 3.0*gsl.MinFloat64 {
		return UnderflowError(result)
	} else if ax < 0.25 {
		var (
			eax = math.Exp(-ax)
			y   = x * x
			c1  = 1.0 / 10.0
			c2  = 1.0 / 280.0
			c3  = 1.0 / 15120.0
			c4  = 1.0 / 1330560.0
			c5  = 1.0 / 172972800.0
			sum = 1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*c5))))
		)
		result.val = eax * x / 3.0 * sum
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		ex := math.Exp(-2.0 * ax)
		result.val = 0.5 * (ax*(1.0+ex) - (1.0 - ex)) / (ax * ax)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		if x < 0.0 {
			result.val = -result.val
		}
		return nil
	}
}

func Bessel_i2_scaled_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if ax < 4.0*gsl.SqrtMinFloat64 {
		return UnderflowError(result)
	} else if ax < 0.25 {
		var (
			y   = x * x
			c1  = 1.0 / 14.0
			c2  = 1.0 / 504.0
			c3  = 1.0 / 33264.0
			c4  = 1.0 / 3459456.0
			c5  = 1.0 / 518918400.0
			sum = 1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*c5))))
			pre = math.Exp(-ax) * x * x / 15.0
		)
		result.val = pre * sum
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		ex := math.Exp(-2.0 * ax)
		x2 := x * x
		result.val = 0.5 * ((3.0+x2)*(1.0-ex) - 3.0*ax*(1.0+ex)) / (ax * ax * ax)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Bessel_il_scaled_e(l int, x float64, result *Result) err.GSLError {
	sgn := 1.0
	ax := math.Abs(x)
	fl := float64(l)
	if x < 0.0 {
		/* i_l(-x) = (-1)^l i_l(x) */
		sgn = 1.
		if gsl.IsOdd(l) {
			sgn = -1.
		}
		x = -x
	}

	if l < 0 {
		return DomainError(result)
	} else if x == 0.0 {
		result.val = 0.
		if l == 0 {
			result.val = 1.
		}
		result.err = 0.0
		return nil
	} else if l == 0 {
		var il Result
		stat_il := Bessel_i0_scaled_e(x, &il)
		result.val = sgn * il.val
		result.err = il.err
		return stat_il
	} else if l == 1 {
		var il Result
		stat_il := Bessel_i1_scaled_e(x, &il)
		result.val = sgn * il.val
		result.err = il.err
		return stat_il
	} else if l == 2 {
		var il Result
		stat_il := Bessel_i2_scaled_e(x, &il)
		result.val = sgn * il.val
		result.err = il.err
		return stat_il
	} else if x*x < 10.0*(fl+1.5)/math.E {
		var b Result
		stat := Bessel_IJ_taylor_e(fl+0.5, x, 1, 50, gsl.Float64Eps, &b)
		pre := math.Exp(-ax) * math.Sqrt((0.5*gsl.Pi)/x)
		result.val = sgn * pre * b.val
		result.err = pre * b.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat
	} else if l < 150 {
		var i0_scaled Result
		stat_i0 := Bessel_i0_scaled_e(ax, &i0_scaled)
		var rat float64
		stat_CF1 := bessel_il_CF1(l, ax, gsl.Float64Eps, &rat)
		iellp1 := rat * gsl.SqrtMinFloat64
		iell := gsl.SqrtMinFloat64
		var iellm1 float64
		var ell int
		for ell = l; ell >= 1; ell-- {
			iellm1 = iellp1 + (2*float64(ell)+1)/x*iell
			iellp1 = iell
			iell = iellm1
		}
		result.val = sgn * i0_scaled.val * (gsl.SqrtMinFloat64 / iell)
		result.err = i0_scaled.err * (gsl.SqrtMinFloat64 / iell)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_i0, stat_CF1)
	} else if gsl.Min(0.29/(fl*fl+1.0), 0.5/(fl*fl+1.0+x*x)) < 0.5*gsl.Root3Float64Eps {
		status := bessel_Inu_scaled_asymp_unif_e(fl+0.5, x, result)
		pre := math.Sqrt((0.5 * gsl.Pi) / x)
		result.val *= sgn * pre
		result.err *= pre
		return status
	} else {
		/* recurse down from safe values */
		rt_term := math.Sqrt((0.5 * gsl.Pi) / x)
		LMAX := 2 + (1.2 / gsl.Root6Float64Eps)
		var r_iellp1, r_iell Result
		stat_a1 := bessel_Inu_scaled_asymp_unif_e(LMAX+1+0.5, x, &r_iellp1)
		stat_a2 := bessel_Inu_scaled_asymp_unif_e(LMAX+0.5, x, &r_iell)
		iellp1 := r_iellp1.val
		iell := r_iell.val
		iellm1 := 0.0
		var ell int
		iellp1 *= rt_term
		iell *= rt_term
		for ell = int(LMAX); ell >= l+1; ell-- {
			iellm1 = iellp1 + (2*float64(ell)+1)/x*iell
			iellp1 = iell
			iell = iellm1
		}
		result.val = sgn * iellm1
		result.err = math.Abs(result.val) * (gsl.Float64Eps + math.Abs(r_iellp1.err/r_iellp1.val) + math.Abs(r_iell.err/r_iell.val))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

		return err.ErrorSelect(stat_a1, stat_a2)
	}
}

func Bessel_il_scaled_array(lmax int, x float64, result_array []float64) err.GSLError {
	if x == 0.0 {
		var ell int
		result_array[0] = 1.0
		for ell = lmax; ell >= 1; ell-- {
			result_array[ell] = 0.0
		}
		return nil
	} else {
		var ell int
		var r_iellp1, r_iell Result
		stat_0 := Bessel_il_scaled_e(lmax+1, x, &r_iellp1)
		stat_1 := Bessel_il_scaled_e(lmax, x, &r_iell)
		iellp1 := r_iellp1.val
		iell := r_iell.val
		var iellm1 float64
		result_array[lmax] = iell
		for ell = lmax; ell >= 1; ell-- {
			iellm1 = iellp1 + (2*float64(ell)+1)/x*iell
			iellp1 = iell
			iell = iellm1
			result_array[ell-1] = iellm1
		}
		return err.ErrorSelect(stat_0, stat_1)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_i0_scaled(x float64) float64 {
	result := new(Result)
	status := Bessel_i0_scaled_e(x, result)
	return EvalResult(result, status)
}

func Bessel_i1_scaled(x float64) float64 {
	result := new(Result)
	status := Bessel_i1_scaled_e(x, result)
	return EvalResult(result, status)
}

func Bessel_i2_scaled(x float64) float64 {
	result := new(Result)
	status := Bessel_i2_scaled_e(x, result)
	return EvalResult(result, status)
}

func Bessel_il_scaled(l int, x float64) float64 {
	result := new(Result)
	status := Bessel_il_scaled_e(l, x, result)
	return EvalResult(result, status)
}
