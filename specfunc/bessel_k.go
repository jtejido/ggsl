/* specfunc/bessel_k.c
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

/* [Abramowitz+Stegun, 10.2.4 + 10.2.6]
 * with lmax=15, precision ~ 15D for x < 3
 *
 * assumes l >= 1
 */
func bessel_kl_scaled_small_x(l int, x float64, result *Result) err.GSLError {
	var num_fact Result
	den := Pow_int(x, l+1)
	stat_df := Doublefact_e(uint(2*l-1), &num_fact)

	if stat_df != nil || den == 0.0 {
		return OverflowError(result)
	} else {
		lmax := 50
		var ipos_term Result
		var ineg_term float64
		sgn := 1.0
		if gsl.IsOdd(l) {
			sgn = -1.0
		}
		ex := math.Exp(x)
		t := 0.5 * x * x
		sum := 1.0
		t_coeff := 1.0
		t_power := 1.0
		var delta float64
		var stat_il err.GSLError
		var i int

		for i = 1; i < lmax; i++ {
			t_coeff /= float64(i * (2*(i-l) - 1))
			t_power *= t
			delta = t_power * t_coeff
			sum += delta
			if math.Abs(delta/sum) < gsl.Float64Eps {
				break
			}
		}

		stat_il = Bessel_il_scaled_e(l, x, &ipos_term)
		ineg_term = sgn * num_fact.val / den * sum
		result.val = -sgn * 0.5 * gsl.Pi * (ex*ipos_term.val - ineg_term)
		result.val *= ex
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_il
	}
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_k0_scaled_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else {
		result.val = gsl.Pi / (2.0 * x)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		CheckUnderflow(result)
		return nil
	}
}

func Bessel_k1_scaled_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if x < (math.SqrtPi+1.0)/(gsl.Sqrt2*gsl.SqrtMaxFloat64) {
		return OverflowError(result)
	} else {
		result.val = gsl.Pi / (2.0 * x) * (1.0 + 1.0/x)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		CheckUnderflow(result)
		return nil
	}
}

func Bessel_k2_scaled_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if x < 2.0/gsl.Root3MaxFloat64 {
		return OverflowError(result)
	} else {
		result.val = gsl.Pi / (2.0 * x) * (1.0 + 3.0/x*(1.0+1.0/x))
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		CheckUnderflow(result)
		return nil
	}
}

func Bessel_kl_scaled_e(l int, x float64, result *Result) err.GSLError {
	fl := float64(l)
	if l < 0 || x <= 0.0 {
		return DomainError(result)
	} else if l == 0 {
		return Bessel_k0_scaled_e(x, result)
	} else if l == 1 {
		return Bessel_k1_scaled_e(x, result)
	} else if l == 2 {
		return Bessel_k2_scaled_e(x, result)
	} else if x < 3.0 {
		return bessel_kl_scaled_small_x(l, x, result)
	} else if gsl.Root3Float64Eps*x > (fl*fl + fl + 1) {
		status := bessel_Knu_scaled_asympx_e(fl+0.5, x, result)
		pre := math.Sqrt((0.5 * gsl.Pi) / x)
		result.val *= pre
		result.err *= pre
		return status
	} else if gsl.Min(0.29/(fl*fl+1.0), 0.5/(fl*fl+1.0+x*x)) < gsl.Root3Float64Eps {
		status := bessel_Knu_scaled_asymp_unif_e(fl+0.5, x, result)
		pre := math.Sqrt((0.5 * gsl.Pi) / x)
		result.val *= pre
		result.err *= pre
		return status
	} else {
		/* recurse upward */
		var r_bk, r_bkm Result
		stat_1 := Bessel_k1_scaled_e(x, &r_bk)
		stat_0 := Bessel_k0_scaled_e(x, &r_bkm)
		var bkp float64
		bk := r_bk.val
		bkm := r_bkm.val
		var j int
		for j = 1; j < l; j++ {
			bkp = (2*float64(j)+1)/x*bk + bkm
			bkm = bk
			bk = bkp
		}
		result.val = bk
		result.err = math.Abs(bk) * (math.Abs(r_bk.err/r_bk.val) + math.Abs(r_bkm.err/r_bkm.val))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

		return err.ErrorSelect(stat_1, stat_0)
	}
}

func Bessel_kl_scaled_array(lmax int, x float64, result_array []float64) err.GSLError {
	if lmax < 0 || x <= 0.0 {
		return err.ERROR("domain error", err.EDOM)
	} else if lmax == 0 {
		var result Result
		stat := Bessel_k0_scaled_e(x, &result)
		result_array[0] = result.val
		return stat
	} else {
		var ell int
		var kellp1, kell, kellm1 float64
		var r_kell, r_kellm1 Result
		Bessel_k1_scaled_e(x, &r_kell)
		Bessel_k0_scaled_e(x, &r_kellm1)
		kell = r_kell.val
		kellm1 = r_kellm1.val
		result_array[0] = kellm1
		result_array[1] = kell
		for ell = 1; ell < lmax; ell++ {
			kellp1 = (2*float64(ell)+1)/x*kell + kellm1
			result_array[ell+1] = kellp1
			kellm1 = kell
			kell = kellp1
		}
		return nil
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_k0_scaled(x float64) float64 {
	result := new(Result)
	status := Bessel_k0_scaled_e(x, result)
	return EvalResult(result, status)
}

func Bessel_k1_scaled(x float64) float64 {
	result := new(Result)
	status := Bessel_k1_scaled_e(x, result)
	return EvalResult(result, status)
}

func Bessel_k2_scaled(x float64) float64 {
	result := new(Result)
	status := Bessel_k2_scaled_e(x, result)
	return EvalResult(result, status)
}

func Bessel_kl_scaled(l int, x float64) float64 {
	result := new(Result)
	status := Bessel_kl_scaled_e(l, x, result)
	return EvalResult(result, status)
}
