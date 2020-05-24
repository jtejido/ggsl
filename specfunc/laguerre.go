/* specfunc/laguerre.c
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

/* based on the large 2b-4a asymptotic for 1F1
 * [Abramowitz+Stegun, 13.5.21]
 * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
 */
func laguerre_large_n(n int, alpha, x float64, result *Result) err.GSLError {
	a := -n
	b := alpha + 1.0
	eta := 2.0*b - 4.0*float64(a)
	cos2th := x / eta
	sin2th := 1.0 - cos2th
	eps := math.Asin(math.Sqrt(cos2th))
	pre_h := 0.25 * gsl.Pi * gsl.Pi * eta * eta * cos2th * sin2th
	lg_b, lnfact := new(Result), new(Result)
	stat_lg := Lngamma_e(b+float64(n), lg_b)
	stat_lf := Lnfact_e(uint(n), lnfact)
	pre_term1 := 0.5 * (1.0 - b) * math.Log(0.25*x*eta)
	pre_term2 := 0.25 * math.Log(pre_h)
	lnpre_val := lg_b.val - lnfact.val + 0.5*x + pre_term1 - pre_term2
	lnpre_err := lg_b.err + lnfact.err + gsl.Float64Eps*(math.Abs(pre_term1)+math.Abs(pre_term2))

	phi1 := 0.25 * eta * (2*eps + math.Sin(2.0*eps))
	ser_term1 := -math.Sin(phi1)

	A1 := (1.0 / 12.0) * (5.0/(4.0*sin2th) + (3.0*b*b-6.0*b+2.0)*sin2th - 1.0)
	ser_term2 := -A1 * math.Cos(phi1) / (0.25 * eta * math.Sin(2.0*eps))

	ser_val := ser_term1 + ser_term2
	ser_err := ser_term2*ser_term2 + gsl.Float64Eps*(math.Abs(ser_term1)+math.Abs(ser_term2))
	stat_e := Exp_mult_err_e(lnpre_val, lnpre_err, ser_val, ser_err, result)
	result.err += 2.0 * gsl.SqrtFloat64Eps * math.Abs(result.val)
	return err.ErrorSelect(stat_e, stat_lf, stat_lg)
}

/* Evaluate polynomial based on confluent hypergeometric representation.
 *
 * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
 *
 * assumes n > 0 and a != negative integer greater than -n
 */
func laguerre_n_cp(n int, a, x float64, result *Result) err.GSLError {
	lnfact, lg1, lg2 := new(Result), new(Result), new(Result)
	var s1, s2 float64
	stat_f := Lnfact_e(uint(n), lnfact)
	stat_g1 := Lngamma_sgn_e(a+1.0+float64(n), lg1, &s1)
	stat_g2 := Lngamma_sgn_e(a+1.0, lg2, &s2)
	poly_1F1_val := 1.0
	poly_1F1_err := 0.0
	var stat_e err.GSLError

	lnpre_val := (lg1.val - lg2.val) - lnfact.val
	lnpre_err := lg1.err + lg2.err + lnfact.err + 2.0*gsl.Float64Eps*math.Abs(lnpre_val)

	for k := n - 1; k >= 0; k-- {
		kk := float64(k)
		nn := float64(n)
		t := (-nn + kk) / (a + 1.0 + kk) * (x / (kk + 1))
		r := t + 1.0/poly_1F1_val
		if r > 0.9*gsl.MaxFloat64/poly_1F1_val {
			/* internal error only, don't call the error handler */
			return InternalOverflowError(result)
		} else {
			/* Collect the Horner terms. */
			poly_1F1_val = 1.0 + t*poly_1F1_val
			poly_1F1_err += gsl.Float64Eps + math.Abs(t)*poly_1F1_err
		}
	}

	stat_e = Exp_mult_err_e(lnpre_val, lnpre_err, poly_1F1_val, poly_1F1_err, result)
	return err.ErrorSelect(stat_e, stat_f, stat_g1, stat_g2)
}

/* Evaluate the polynomial based on the confluent hypergeometric
 * function in a safe way, with no restriction on the arguments.
 *
 * assumes x != 0
 */
func laguerre_n_poly_safe(n int, a, x float64, result *Result) err.GSLError {
	b := a + 1.0
	mx := -x
	tc_sgn := 1.

	if x >= 0. && gsl.IsOdd(n) {
		tc_sgn = -1
	}
	tc := new(Result)
	stat_tc := Taylorcoeff_e(n, math.Abs(x), tc)

	if stat_tc == nil {
		term := tc.val * tc_sgn
		sum_val := term
		sum_err := tc.err

		for k := n - 1; k >= 0; k-- {
			term *= ((b + float64(k)) / float64(n-k)) * (float64(k) + 1.0) / mx
			sum_val += term
			sum_err += 4.0 * gsl.Float64Eps * math.Abs(term)
		}
		result.val = sum_val
		result.err = sum_err + 2.0*gsl.Float64Eps*math.Abs(result.val)
		return nil
	} else if stat_tc.Status() == err.EOVRFLW {
		result.val = 0.0 /* FIXME: should be Inf */
		result.err = 0.0
		return stat_tc
	} else {
		result.val = 0.0
		result.err = 0.0
		return stat_tc
	}

}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*/
func Laguerre_1_e(a, x float64, result *Result) err.GSLError {
	result.val = 1.0 + a - x
	result.err = 2.0 * gsl.Float64Eps * (1.0 + math.Abs(a) + math.Abs(x))
	return nil
}

func Laguerre_2_e(a, x float64, result *Result) err.GSLError {

	if a == -2.0 {
		result.val = 0.5 * x * x
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}

	c0 := 0.5 * (2.0 + a) * (1.0 + a)
	c1 := -(2.0 + a)
	c2 := -0.5 / (2.0 + a)
	result.val = c0 + c1*x*(1.0+c2*x)
	result.err = 2.0 * gsl.Float64Eps * (math.Abs(c0) + 2.0*math.Abs(c1*x)*(1.0+2.0*math.Abs(c2*x)))
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return nil

}

func Laguerre_3_e(a, x float64, result *Result) err.GSLError {
	if a == -2.0 {
		x2_6 := x * x / 6.0
		result.val = x2_6 * (3.0 - x)
		result.err = x2_6 * (3.0 + math.Abs(x)) * 2.0 * gsl.Float64Eps
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if a == -3.0 {
		result.val = -x * x / 6.0
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}

	c0 := (3.0 + a) * (2.0 + a) * (1.0 + a) / 6.0
	c1 := -c0 * 3.0 / (1.0 + a)
	c2 := -1.0 / (2.0 + a)
	c3 := -1.0 / (3.0 * (3.0 + a))
	result.val = c0 + c1*x*(1.0+c2*x*(1.0+c3*x))
	result.err = 1.0 + 2.0*math.Abs(c3*x)
	result.err = 1.0 + 2.0*math.Abs(c2*x)*result.err
	result.err = 2.0 * gsl.Float64Eps * (math.Abs(c0) + 2.0*math.Abs(c1*x)*result.err)
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return nil

}

func Laguerre_n_e(n int, a, x float64, result *Result) err.GSLError {
	if n < 0 {
		return DomainError(result)
	} else if n == 0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if n == 1 {
		result.val = 1.0 + a - x
		result.err = 2.0 * gsl.Float64Eps * (1.0 + math.Abs(a) + math.Abs(x))
		return nil
	} else if x == 0.0 {
		product := a + 1.0

		for k := 2; k <= n; k++ {
			product *= (a + float64(k)) / float64(k)
		}
		result.val = product
		result.err = 2.0*(float64(n)+1.0)*gsl.Float64Eps*math.Abs(product) + gsl.Float64Eps
		return nil
	} else if x < 0.0 && a > -1.0 {
		/* In this case all the terms in the polynomial
		 * are of the same sign. Note that this also
		 * catches overflows correctly.
		 */

		return laguerre_n_cp(n, a, x, result)
	} else if n < 5 || (x > 0.0 && a < -float64(n)-1) {
		/* Either the polynomial will not lose too much accuracy
		 * or all the terms are negative. In any case,
		 * the error estimate here is good. We try both
		 * explicit summation methods, as they have different
		 * characteristics. One may underflow/overflow while the
		 * other does not.
		 */
		if laguerre_n_cp(n, a, x, result) == nil {
			return nil
		} else {
			return laguerre_n_poly_safe(n, a, x, result)
		}
	} else if float64(n) > 1.0e+07 && x > 0.0 && a > -1.0 && x < 2.0*(a+1.0)+4.0*float64(n) {
		return laguerre_large_n(n, a, x, result)
	} else if a >= 0.0 || (x > 0.0 && a < -float64(n)-1) {

		lg2 := new(Result)
		stat_lg2 := Laguerre_2_e(a, x, lg2)
		Lkm1 := 1.0 + a - x
		Lk := lg2.val

		for k := 2; k < n; k++ {
			kk := float64(k)
			Lkp1 := (-(kk+a)*Lkm1 + (2.0*kk+a+1.0-x)*Lk) / (kk + 1.0)
			Lkm1 = Lk
			Lk = Lkp1
		}
		result.val = Lk
		result.err = (math.Abs(lg2.err/lg2.val) + gsl.Float64Eps) * float64(n) * math.Abs(Lk)
		return stat_lg2
	}

	/* Despair... or magic? */
	return laguerre_n_poly_safe(n, a, x, result)

}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Laguerre_1(a, x float64) float64 {
	result := new(Result)
	status := Laguerre_1_e(a, x, result)
	return EvalResult(result, status)
}

func Laguerre_2(a, x float64) float64 {
	result := new(Result)
	status := Laguerre_2_e(a, x, result)
	return EvalResult(result, status)
}

func Laguerre_3(a, x float64) float64 {
	result := new(Result)
	status := Laguerre_3_e(a, x, result)
	return EvalResult(result, status)
}

func Laguerre_n(n int, a, x float64) float64 {
	result := new(Result)
	status := Laguerre_n_e(n, a, x, result)
	return EvalResult(result, status)
}
