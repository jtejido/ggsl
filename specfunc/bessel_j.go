/* specfunc/bessel_j.c
 *
 * Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003 Gerard Jungman
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
func Bessel_j0_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	/* CHECK_POINTER(result) */

	if ax < 0.5 {
		var (
			y  = x * x
			c1 = -1.0 / 6.0
			c2 = 1.0 / 120.0
			c3 = -1.0 / 5040.0
			c4 = 1.0 / 362880.0
			c5 = -1.0 / 39916800.0
			c6 = 1.0 / 6227020800.0
		)
		result.val = 1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*c6)))))
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = math.Sin(x) / x
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Bessel_j1_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if ax < 3.1*gsl.MinFloat64 {
		return UnderflowError(result)
	} else if ax < 0.25 {
		var (
			y   = x * x
			c1  = -1.0 / 10.0
			c2  = 1.0 / 280.0
			c3  = -1.0 / 15120.0
			c4  = 1.0 / 1330560.0
			c5  = -1.0 / 172972800.0
			sum = 1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*c5))))
		)
		result.val = x / 3.0 * sum
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		cos_x := math.Cos(x)
		sin_x := math.Sin(x)
		result.val = (sin_x/x - cos_x) / x
		result.err = 2.0 * gsl.Float64Eps * (math.Abs(sin_x/(x*x)) + math.Abs(cos_x/x))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Bessel_j2_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if ax < 4.0*gsl.SqrtMinFloat64 {
		return UnderflowError(result)
	} else if ax < 1.3 {
		var (
			y   = x * x
			c1  = -1.0 / 14.0
			c2  = 1.0 / 504.0
			c3  = -1.0 / 33264.0
			c4  = 1.0 / 3459456.0
			c5  = -1.0 / 518918400
			c6  = 1.0 / 105859353600.0
			c7  = -1.0 / 28158588057600.0
			c8  = 1.0 / 9461285587353600.0
			c9  = -1.0 / 3916972233164390400.0
			sum = 1.0 + y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*(c8+y*c9))))))))
		)
		result.val = y / 15.0 * sum
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		/* bug #45730: switch from gsl_sf_{cos,sin} to cos/sin to fix large inputs */
		// #if 0
		//     gsl_sf_result cos_result;
		//     gsl_sf_result sin_result;
		//     const int stat_cos = gsl_sf_cos_e(x, &cos_result);
		//     const int stat_sin = gsl_sf_sin_e(x, &sin_result);
		//     const double cos_x = cos_result.val;
		//     const double sin_x = sin_result.val;
		//     const double f = (3.0/(x*x) - 1.0);
		//     result.val  = (f * sin_x - 3.0*cos_x/x)/x;
		//     result.err  = math.Abs(f * sin_result.err/x) + math.Abs((3.0*cos_result.err/x)/x);
		//     result.err += 2.0 * gsl.Float64Eps * (math.Abs(f*sin_x/x) + 3.0*math.Abs(cos_x/(x*x)));
		//     result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val);
		//     return GSL_ERROR_SELECT_2(stat_cos, stat_sin);
		// #else
		cos_x := math.Cos(x)
		sin_x := math.Sin(x)
		f := (3.0/(x*x) - 1.0)
		result.val = (f*sin_x - 3.0*cos_x/x) / x
		result.err = 2.0 * gsl.Float64Eps * (math.Abs(f*sin_x/x) + 3.0*math.Abs(cos_x/(x*x)))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		/*return GSL_ERROR_SELECT_2(stat_cos, stat_sin);*/
		return nil
	}
}

func Bessel_jl_e(l int, x float64, result *Result) err.GSLError {
	fl := float64(l)
	if l < 0 || x < 0.0 {
		return DomainError(result)
	} else if x == 0.0 {
		result.val = 1.0
		if l > 0 {
			result.val = 0.0
		}
		result.err = 0.0
		return nil
	} else if l == 0 {
		return Bessel_j0_e(x, result)
	} else if l == 1 {
		return Bessel_j1_e(x, result)
	} else if l == 2 {
		return Bessel_j2_e(x, result)
	} else if x*x < 10.0*(fl+0.5)/math.E {
		var b Result
		status := Bessel_IJ_taylor_e(fl+0.5, x, -1, 50, gsl.Float64Eps, &b)
		pre := math.Sqrt((0.5 * gsl.Pi) / x)
		result.val = pre * b.val
		result.err = pre * b.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return status
	} else if gsl.Root4Float64Eps*x > (fl*fl + fl + 1.0) {
		var b Result
		status := bessel_Jnu_asympx_e(fl+0.5, x, &b)
		pre := math.Sqrt((0.5 * gsl.Pi) / x)
		result.val = pre * b.val
		result.err = 2.0*gsl.Float64Eps*math.Abs(result.val) + pre*b.err
		return status
	} else if fl > 1.0/gsl.Root6Float64Eps {
		var b Result
		status := bessel_Jnu_asymp_Olver_e(fl+0.5, x, &b)
		pre := math.Sqrt((0.5 * gsl.Pi) / x)
		result.val = pre * b.val
		result.err = 2.0*gsl.Float64Eps*math.Abs(result.val) + pre*b.err
		return status
	} else if x > 1000.0 && x > fl*fl {
		/* We need this path to avoid feeding large x to CF1 below; */
		var b Result
		status := bessel_Jnu_asympx_e(fl+0.5, x, &b)
		pre := math.Sqrt((0.5 * gsl.Pi) / x)
		result.val = pre * b.val
		result.err = 2.0*gsl.Float64Eps*math.Abs(result.val) + pre*b.err
		return status
	} else {
		var sgn, ratio float64
		/* The CF1 call will hit 10000 iterations for x > 10000 + l */
		stat_CF1 := bessel_J_CF1(fl+0.5, x, &ratio, &sgn)
		BESSEL_J_SMALL := gsl.MinFloat64 / gsl.Float64Eps
		jellp1 := BESSEL_J_SMALL * ratio
		jell := BESSEL_J_SMALL
		var jellm1 float64
		var ell int
		for ell = l; ell > 0; ell-- {
			jellm1 = -jellp1 + (2*float64(ell)+1)/x*jell
			jellp1 = jell
			jell = jellm1
		}

		if math.Abs(jell) > math.Abs(jellp1) {
			var j0_result Result
			stat_j0 := Bessel_j0_e(x, &j0_result)
			pre := BESSEL_J_SMALL / jell
			result.val = j0_result.val * pre
			result.err = j0_result.err * math.Abs(pre)
			result.err += 4.0 * gsl.Float64Eps * (0.5*fl + 1.0) * math.Abs(result.val)
			return err.ErrorSelect(stat_j0, stat_CF1)
		} else {
			var j1_result Result
			stat_j1 := Bessel_j1_e(x, &j1_result)
			pre := BESSEL_J_SMALL / jellp1
			result.val = j1_result.val * pre
			result.err = j1_result.err * math.Abs(pre)
			result.err += 4.0 * gsl.Float64Eps * (0.5*fl + 1.0) * math.Abs(result.val)
			return err.ErrorSelect(stat_j1, stat_CF1)
		}
	}
}

func Bessel_jl_array(lmax int, x float64, result_array []float64) err.GSLError {

	if lmax < 0 || x < 0.0 {
		var j int
		for j = 0; j <= lmax; j++ {
			result_array[j] = 0.0
		}
		return err.ERROR("error", err.EDOM)
	} else if x == 0.0 {
		var j int
		for j = 1; j <= lmax; j++ {
			result_array[j] = 0.0
		}
		result_array[0] = 1.0
		return nil
	} else {
		var r_jellp1, r_jell Result
		stat_0 := Bessel_jl_e(lmax+1, x, &r_jellp1)
		stat_1 := Bessel_jl_e(lmax, x, &r_jell)
		jellp1 := r_jellp1.val
		jell := r_jell.val
		var jellm1 float64
		var ell int

		result_array[lmax] = jell
		for ell = lmax; ell >= 1; ell-- {
			jellm1 = -jellp1 + (2*float64(ell)+1)/x*jell
			jellp1 = jell
			jell = jellm1
			result_array[ell-1] = jellm1
		}

		return err.ErrorSelect(stat_0, stat_1)
	}
}

func Bessel_jl_steed_array(lmax int, x float64, jl_x []float64) err.GSLError {

	if lmax < 0 || x < 0.0 {
		var j int
		for j = 0; j <= lmax; j++ {
			jl_x[j] = 0.0
		}
		return err.ERROR("error", err.EDOM)
	} else if x == 0.0 {
		var j int
		for j = 1; j <= lmax; j++ {
			jl_x[j] = 0.0
		}
		jl_x[0] = 1.0
		return nil
	} else if x < 2.0*gsl.Root4Float64Eps {
		/* first two terms of Taylor series */
		inv_fact := 1.0 /* 1/(1 3 5 ... (2l+1)) */
		x_l := 1.0      /* x^l */
		var l int
		for l = 0; l <= lmax; l++ {
			jl_x[l] = x_l * inv_fact
			jl_x[l] *= 1.0 - 0.5*x*x/(2.0*float64(l)+3.0)
			inv_fact /= 2.0*float64(l) + 3.0
			x_l *= x
		}
		return nil
	} else {
		/* Steed/Barnett algorithm [Comp. Phys. Comm. 21, 297 (1981)] */
		x_inv := 1.0 / x
		W := 2.0 * x_inv
		F := 1.0
		FP := (float64(lmax) + 1.0) * x_inv
		B := 2.0*FP + x_inv
		end := B + 20000.0*W
		D := 1.0 / B
		del := -D

		FP += del

		/* continued fraction */
		for ok := true; ok; ok = math.Abs(del) >= math.Abs(FP)*gsl.Float64Eps {
			B += W
			D = 1.0 / (B - D)
			del *= (B*D - 1.)
			FP += del
			if D < 0.0 {
				F = -F
			}
			if B > end {
				return err.ERROR("error", err.EMAXITER)
			}
		}

		FP *= F

		if lmax > 0 {
			/* downward recursion */
			XP2 := FP
			PL := float64(lmax) * x_inv
			L := lmax
			var LP int
			jl_x[lmax] = F
			for LP = 1; LP <= lmax; LP++ {
				jl_x[L-1] = PL*jl_x[L] + XP2
				FP = PL*jl_x[L-1] - jl_x[L]
				XP2 = FP
				PL -= x_inv
				L--
			}
			F = jl_x[0]
		}

		/* normalization */
		W = x_inv / math.Hypot(FP, F)
		jl_x[0] = W * F
		if lmax > 0 {
			var L int
			for L = 1; L <= lmax; L++ {
				jl_x[L] *= W
			}
		}

		return nil
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_j0(x float64) float64 {
	result := new(Result)
	status := Bessel_j0_e(x, result)
	return EvalResult(result, status)
}

func Bessel_j1(x float64) float64 {
	result := new(Result)
	status := Bessel_j1_e(x, result)
	return EvalResult(result, status)
}

func Bessel_j2(x float64) float64 {
	result := new(Result)
	status := Bessel_j2_e(x, result)
	return EvalResult(result, status)
}

func Bessel_jl(l int, x float64) float64 {
	result := new(Result)
	status := Bessel_jl_e(l, x, result)
	return EvalResult(result, status)
}
