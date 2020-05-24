/* specfunc/dilog.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
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

/* Evaluate series for real dimath.Log(x)
 * Sum[ x^k / k^2, {k,1,Infinity}]
 *
 * Converges rapidly for |x| < 1/2.
 */
func dilog_series_1(x float64, result *Result) err.GSLError {
	kmax := 1000
	sum := x
	term := x
	var k int

	for k = 2; k < kmax; k++ {
		rk := (float64(k) - 1.0) / float64(k)
		term *= x
		term *= rk * rk
		sum += term
		if math.Abs(term/sum) < gsl.Float64Eps {
			break
		}
	}

	result.val = sum
	result.err = 2.0 * math.Abs(term)
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

	if k == kmax {
		return err.ERROR("error", err.EMAXITER)
	}

	return nil
}

/* Compute the associated series
 *
 *   sum_{k=1}{infty} r^k / (k^2 (k+1))
 *
 * This is a series which appears in the one-step accelerated
 * method, which splits out one elementary function from the
 * full definition of Li_2(x). See below.
 */
func series_2(r float64, result *Result) err.GSLError {
	kmax := 100
	rk := r
	sum := 0.5 * r
	var k int

	for k = 2; k < 10; k++ {
		var ds float64
		fk := float64(k)
		rk *= r
		ds = rk / (fk * fk * (fk + 1.0))
		sum += ds
	}
	for ; k < kmax; k++ {
		var ds float64
		fk := float64(k)
		rk *= r
		ds = rk / (fk * fk * (fk + 1.0))
		sum += ds
		if math.Abs(ds/sum) < 0.5*gsl.Float64Eps {
			break
		}
	}

	result.val = sum
	result.err = 2.0 * float64(kmax) * gsl.Float64Eps * math.Abs(sum)

	return nil
}

/* Compute Li_2(x) using the accelerated series representation.
 *
 * Li_2(x) = 1 + (1-x)ln(1-x)/x + series_2(x)
 *
 * assumes: -1 < x < 1
 */
func dilog_series_2(x float64, result *Result) err.GSLError {
	stat_s3 := series_2(x, result)
	var t float64
	if x > 0.01 {
		t = (1.0 - x) * math.Log(1.0-x) / x
	} else {
		c3 := 1.0 / 3.0
		c4 := 1.0 / 4.0
		c5 := 1.0 / 5.0
		c6 := 1.0 / 6.0
		c7 := 1.0 / 7.0
		c8 := 1.0 / 8.0
		t68 := c6 + x*(c7+x*c8)
		t38 := c3 + x*(c4+x*(c5+x*t68))
		t = (x - 1.0) * (1.0 + x*(0.5+x*t38))
	}
	result.val += 1.0 + t
	result.err += 2.0 * gsl.Float64Eps * math.Abs(t)
	return stat_s3
}

/* Calculates Li_2(x) for real x. Assumes x >= 0.0.
 */
func dilog_xge0(x float64, result *Result) err.GSLError {
	if x > 2.0 {
		var ser Result
		stat_ser := dilog_series_2(1.0/x, &ser)
		log_x := math.Log(x)
		t1 := gsl.Pi * gsl.Pi / 3.0
		t2 := ser.val
		t3 := 0.5 * log_x * log_x
		result.val = t1 - t2 - t3
		result.err = gsl.Float64Eps*math.Abs(log_x) + ser.err
		result.err += gsl.Float64Eps * (math.Abs(t1) + math.Abs(t2) + math.Abs(t3))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_ser
	} else if x > 1.01 {
		var ser Result
		stat_ser := dilog_series_2(1.0-1.0/x, &ser)
		log_x := math.Log(x)
		log_term := log_x * (math.Log(1.0-1.0/x) + 0.5*log_x)
		t1 := gsl.Pi * gsl.Pi / 6.0
		t2 := ser.val
		t3 := log_term
		result.val = t1 + t2 - t3
		result.err = gsl.Float64Eps*math.Abs(log_x) + ser.err
		result.err += gsl.Float64Eps * (math.Abs(t1) + math.Abs(t2) + math.Abs(t3))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_ser
	} else if x > 1.0 {
		/* series around x = 1.0 */
		eps := x - 1.0
		lne := math.Log(eps)
		c0 := gsl.Pi * gsl.Pi / 6.0
		c1 := 1.0 - lne
		c2 := -(1.0 - 2.0*lne) / 4.0
		c3 := (1.0 - 3.0*lne) / 9.0
		c4 := -(1.0 - 4.0*lne) / 16.0
		c5 := (1.0 - 5.0*lne) / 25.0
		c6 := -(1.0 - 6.0*lne) / 36.0
		c7 := (1.0 - 7.0*lne) / 49.0
		c8 := -(1.0 - 8.0*lne) / 64.0
		result.val = c0 + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8)))))))
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x == 1.0 {
		result.val = gsl.Pi * gsl.Pi / 6.0
		result.err = 2.0 * gsl.Float64Eps * gsl.Pi * gsl.Pi / 6.0
		return nil
	} else if x > 0.5 {
		var ser Result
		stat_ser := dilog_series_2(1.0-x, &ser)
		log_x := math.Log(x)
		t1 := gsl.Pi * gsl.Pi / 6.0
		t2 := ser.val
		t3 := log_x * math.Log(1.0-x)
		result.val = t1 - t2 - t3
		result.err = gsl.Float64Eps*math.Abs(log_x) + ser.err
		result.err += gsl.Float64Eps * (math.Abs(t1) + math.Abs(t2) + math.Abs(t3))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_ser
	} else if x > 0.25 {
		return dilog_series_2(x, result)
	} else if x > 0.0 {
		return dilog_series_1(x, result)
	} else {
		/* x == 0.0 */
		result.val = 0.0
		result.err = 0.0
		return nil
	}
}

/* Evaluate the series representation for Li2(z):
 *
 *   Li2(z) = Sum[ |z|^k / k^2 Exp[i k arg(z)], {k,1,Infinity}]
 *   |z|    = r
 *   arg(z) = theta
 *
 * Assumes 0 < r < 1.
 * It is used only for small r.
 */
func dilogc_series_1(r, x, y float64, real_result, imag_result *Result) err.GSLError {
	cos_theta := x / r
	sin_theta := y / r
	alpha := 1.0 - cos_theta
	beta := sin_theta
	ck := cos_theta
	sk := sin_theta
	rk := r
	real_sum := r * ck
	imag_sum := r * sk
	kmax := 50 + int(22.0/(-math.Log(r))) /* tuned for double-precision */

	for k := 2; k < kmax; k++ {
		var dr, di float64
		ck_tmp := ck
		ck = ck - (alpha*ck + beta*sk)
		sk = sk - (alpha*sk - beta*ck_tmp)
		rk *= r
		dr = rk / float64(k*k) * ck
		di = rk / float64(k*k) * sk
		real_sum += dr
		imag_sum += di
		if math.Abs((dr*dr+di*di)/(real_sum*real_sum+imag_sum*imag_sum)) < gsl.Float64Eps*gsl.Float64Eps {
			break
		}
	}

	real_result.val = real_sum
	real_result.err = 2.0 * float64(kmax) * gsl.Float64Eps * math.Abs(real_sum)
	imag_result.val = imag_sum
	imag_result.err = 2.0 * float64(kmax) * gsl.Float64Eps * math.Abs(imag_sum)

	return nil
}

/* Compute
 *
 *   sum_{k=1}{infty} z^k / (k^2 (k+1))
 *
 * This is a series which appears in the one-step accelerated
 * method, which splits out one elementary function from the
 * full definition of Li_2.
 */
func series_2_c(r, x, y float64, sum_re, sum_im *Result) err.GSLError {
	cos_theta := x / r
	sin_theta := y / r
	alpha := 1.0 - cos_theta
	beta := sin_theta
	ck := cos_theta
	sk := sin_theta
	rk := r
	real_sum := 0.5 * r * ck
	imag_sum := 0.5 * r * sk
	kmax := 30 + (int)(18.0/(-math.Log(r))) /* tuned for double-precision */

	for k := 2; k < kmax; k++ {
		var dr, di float64
		ck_tmp := ck
		ck = ck - (alpha*ck + beta*sk)
		sk = sk - (alpha*sk - beta*ck_tmp)
		rk *= r
		dr = rk / float64(k*k*(k+1)) * ck
		di = rk / float64(k*k*(k+1)) * sk
		real_sum += dr
		imag_sum += di
		if math.Abs((dr*dr+di*di)/(real_sum*real_sum+imag_sum*imag_sum)) < gsl.Float64Eps*gsl.Float64Eps {
			break
		}
	}

	sum_re.val = real_sum
	sum_re.err = 2.0 * float64(kmax) * gsl.Float64Eps * math.Abs(real_sum)
	sum_im.val = imag_sum
	sum_im.err = 2.0 * float64(kmax) * gsl.Float64Eps * math.Abs(imag_sum)

	return nil
}

/* Compute Li_2(z) using the one-step accelerated series.
 *
 * Li_2(z) = 1 + (1-z)ln(1-z)/z + series_2_c(z)
 *
 * z = r exp(i theta)
 * assumes: r < 1
 * assumes: r > epsilon, so that we take no special care with math.Log(1-z)
 */
func dilogc_series_2(r, x, y float64, real_dl, imag_dl *Result) err.GSLError {
	if r == 0.0 {
		real_dl.val = 0.0
		imag_dl.val = 0.0
		real_dl.err = 0.0
		imag_dl.err = 0.0
		return nil
	} else {
		var sum_re, sum_im Result
		stat_s3 := series_2_c(r, x, y, &sum_re, &sum_im)

		/* t = ln(1-z)/z */
		var ln_omz_r, ln_omz_theta Result
		stat_log := Complex_log_e(1.0-x, -y, &ln_omz_r, &ln_omz_theta)
		t_x := (ln_omz_r.val*x + ln_omz_theta.val*y) / (r * r)
		t_y := (-ln_omz_r.val*y + ln_omz_theta.val*x) / (r * r)

		/* r = (1-z) ln(1-z)/z */
		r_x := (1.0-x)*t_x + y*t_y
		r_y := (1.0-x)*t_y - y*t_x

		real_dl.val = sum_re.val + r_x + 1.0
		imag_dl.val = sum_im.val + r_y
		real_dl.err = sum_re.err + 2.0*gsl.Float64Eps*(math.Abs(real_dl.val)+math.Abs(r_x))
		imag_dl.err = sum_im.err + 2.0*gsl.Float64Eps*(math.Abs(imag_dl.val)+math.Abs(r_y))
		return err.ErrorSelect(stat_s3, stat_log)
	}
}

/* Evaluate a series for Li_2(z) when |z| is near 1.
 * This is uniformly good away from z=1.
 *
 *   Li_2(z) = Sum[ a^n/n! H_n(theta), {n, 0, Infinity}]
 *
 * where
 *   H_n(theta) = Sum[ e^(i m theta) m^n / m^2, {m, 1, Infinity}]
 *   a = ln(r)
 *
 *  H_0(t) = Gl_2(t) + i Cl_2(t)
 *  H_1(t) = 1/2 ln(2(1-c)) + I atan2(-s, 1-c)
 *  H_2(t) = -1/2 + I/2 s/(1-c)
 *  H_3(t) = -1/2 /(1-c)
 *  H_4(t) = -I/2 s/(1-c)^2
 *  H_5(t) = 1/2 (2 + c)/(1-c)^2
 *  H_6(t) = I/2 s/(1-c)^5 (8(1-c) - s^2 (3 + c))
 */
func dilogc_series_3(r, x, y float64, real_result, imag_result *Result) err.GSLError {
	theta := math.Atan2(y, x)
	cos_theta := x / r
	sin_theta := y / r
	a := math.Log(r)
	omc := 1.0 - cos_theta
	omc2 := omc * omc
	var H_re, H_im [7]float64
	var an, nfact, sum_re, sum_im float64
	var Him0 Result

	H_re[0] = gsl.Pi*gsl.Pi/6.0 + 0.25*(theta*theta-2.0*gsl.Pi*math.Abs(theta))
	Clausen_e(theta, &Him0)
	H_im[0] = Him0.val

	H_re[1] = -0.5 * math.Log(2.0*omc)
	H_im[1] = -math.Atan2(-sin_theta, omc)

	H_re[2] = -0.5
	H_im[2] = 0.5 * sin_theta / omc

	H_re[3] = -0.5 / omc
	H_im[3] = 0.0

	H_re[4] = 0.0
	H_im[4] = -0.5 * sin_theta / omc2

	H_re[5] = 0.5 * (2.0 + cos_theta) / omc2
	H_im[5] = 0.0

	H_re[6] = 0.0
	H_im[6] = 0.5 * sin_theta / (omc2 * omc2 * omc) * (8.0*omc - sin_theta*sin_theta*(3.0+cos_theta))

	sum_re = H_re[0]
	sum_im = H_im[0]
	an = 1.0
	nfact = 1.0
	for n := 1; n <= 6; n++ {
		var t float64
		an *= a
		nfact *= float64(n)
		t = an / nfact
		sum_re += t * H_re[n]
		sum_im += t * H_im[n]
	}

	real_result.val = sum_re
	real_result.err = 2.0*6.0*gsl.Float64Eps*math.Abs(sum_re) + math.Abs(an/nfact)
	imag_result.val = sum_im
	imag_result.err = 2.0*6.0*gsl.Float64Eps*math.Abs(sum_im) + Him0.err + math.Abs(an/nfact)

	return nil
}

/* Calculate complex dilogarithm Li_2(z) in the fundamental region,
 * which we take to be the intersection of the unit disk with the
 * half-space x < MAGIC_SPLIT_VALUE. It turns out that 0.732 is a
 * nice choice for MAGIC_SPLIT_VALUE since then points mapped out
 * of the x > MAGIC_SPLIT_VALUE region and into another part of the
 * unit disk are bounded in radius by MAGIC_SPLIT_VALUE itself.
 *
 * If |z| < 0.98 we use a direct series summation. Otherwise z is very
 * near the unit circle, and the series_2 expansion is used; see above.
 * Because the fundamental region is bounded away from z = 1, this
 * works well.
 */
func dilogc_fundamental(r, x, y float64, real_dl, imag_dl *Result) err.GSLError {
	if r > 0.98 {
		return dilogc_series_3(r, x, y, real_dl, imag_dl)
	} else if r > 0.25 {
		return dilogc_series_2(r, x, y, real_dl, imag_dl)
	}

	return dilogc_series_1(r, x, y, real_dl, imag_dl)
}

/* Compute Li_2(z) for z in the unit disk, |z| < 1. If z is outside
 * the fundamental region, which means that it is too close to z = 1,
 * then it is reflected into the fundamental region using the identity
 *
 *   Li2(z) = -Li2(1-z) + zeta(2) - ln(z) ln(1-z).
 */
func dilogc_unitdisk(x, y float64, real_dl, imag_dl *Result) err.GSLError {
	MAGIC_SPLIT_VALUE := 0.732
	zeta2 := gsl.Pi * gsl.Pi / 6.0
	r := math.Hypot(x, y)

	if x > MAGIC_SPLIT_VALUE {
		/* Reflect away from z = 1 if we are too close. The magic value
		 * insures that the reflected value of the radius satisfies the
		 * related inequality r_tmp < MAGIC_SPLIT_VALUE.
		 */
		x_tmp := 1.0 - x
		y_tmp := -y
		r_tmp := math.Hypot(x_tmp, y_tmp)
		/* const double cos_theta_tmp = x_tmp/r_tmp; */
		/* const double sin_theta_tmp = y_tmp/r_tmp; */

		var result_re_tmp, result_im_tmp Result

		stat_dilog := dilogc_fundamental(r_tmp, x_tmp, y_tmp, &result_re_tmp, &result_im_tmp)

		lnz := math.Log(r)                 /*  math.Log(|z|)   */
		lnomz := math.Log(r_tmp)           /*  math.Log(|1-z|) */
		argz := math.Atan2(y, x)           /*  arg(z) assuming principal branch */
		argomz := math.Atan2(y_tmp, x_tmp) /*  arg(1-z)   */
		real_dl.val = -result_re_tmp.val + zeta2 - lnz*lnomz + argz*argomz
		real_dl.err = result_re_tmp.err
		real_dl.err += 2.0 * gsl.Float64Eps * (zeta2 + math.Abs(lnz*lnomz) + math.Abs(argz*argomz))
		imag_dl.val = -result_im_tmp.val - argz*lnomz - argomz*lnz
		imag_dl.err = result_im_tmp.err
		imag_dl.err += 2.0 * gsl.Float64Eps * (math.Abs(argz*lnomz) + math.Abs(argomz*lnz))

		return stat_dilog
	}

	return dilogc_fundamental(r, x, y, real_dl, imag_dl)

}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Dilog_e(x float64, result *Result) err.GSLError {
	if x >= 0.0 {
		return dilog_xge0(x, result)
	} else {
		var d1, d2 Result
		stat_d1 := dilog_xge0(-x, &d1)
		stat_d2 := dilog_xge0(x*x, &d2)
		result.val = -d1.val + 0.5*d2.val
		result.err = d1.err + 0.5*d2.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_d1, stat_d2)
	}
}

func Complex_dilog_xy_e(x, y float64, real_dl, imag_dl *Result) err.GSLError {
	zeta2 := gsl.Pi * gsl.Pi / 6.0
	r2 := x*x + y*y

	if y == 0.0 {
		if x >= 1.0 {
			imag_dl.val = -gsl.Pi * math.Log(x)
			imag_dl.err = 2.0 * gsl.Float64Eps * math.Abs(imag_dl.val)
		} else {
			imag_dl.val = 0.0
			imag_dl.err = 0.0
		}
		return Dilog_e(x, real_dl)
	} else if math.Abs(r2-1.0) < gsl.Float64Eps {
		/* Lewin A.2.4.1 and A.2.4.2 */

		theta := math.Atan2(y, x)
		term1 := theta * theta / 4.0
		term2 := gsl.Pi * math.Abs(theta) / 2.0
		real_dl.val = zeta2 + term1 - term2
		real_dl.err = 2.0 * gsl.Float64Eps * (zeta2 + term1 + term2)
		return Clausen_e(theta, imag_dl)
	} else if r2 < 1.0 {
		return dilogc_unitdisk(x, y, real_dl, imag_dl)
	} else {
		/* Reduce argument to unit disk. */
		r := math.Sqrt(r2)
		x_tmp := x / r2
		y_tmp := -y / r2
		/* const double r_tmp = 1.0/r; */
		var result_re_tmp, result_im_tmp Result

		stat_dilog := dilogc_unitdisk(x_tmp, y_tmp, &result_re_tmp, &result_im_tmp)

		/* Unwind the inversion.
		 *
		 *  Li_2(z) + Li_2(1/z) = -zeta(2) - 1/2 ln(-z)^2
		 */
		theta := math.Atan2(y, x)
		theta_abs := math.Abs(theta)
		theta_sgn := 1.0
		if theta < 0.0 {
			theta_sgn = -1.0
		}
		ln_minusz_re := math.Log(r)
		ln_minusz_im := theta_sgn * (theta_abs - gsl.Pi)
		lmz2_re := ln_minusz_re*ln_minusz_re - ln_minusz_im*ln_minusz_im
		lmz2_im := 2.0 * ln_minusz_re * ln_minusz_im
		real_dl.val = -result_re_tmp.val - 0.5*lmz2_re - zeta2
		real_dl.err = result_re_tmp.err + 2.0*gsl.Float64Eps*(0.5*math.Abs(lmz2_re)+zeta2)
		imag_dl.val = -result_im_tmp.val - 0.5*lmz2_im
		imag_dl.err = result_im_tmp.err + 2.0*gsl.Float64Eps*math.Abs(lmz2_im)
		return stat_dilog
	}
}

func Complex_dilog_e(r, theta float64, real_dl, imag_dl *Result) err.GSLError {
	cos_theta := math.Cos(theta)
	sin_theta := math.Sin(theta)
	x := r * cos_theta
	y := r * sin_theta
	return Complex_dilog_xy_e(x, y, real_dl, imag_dl)
}

func Complex_spence_xy_e(x, y float64, real_sp, imag_sp *Result) err.GSLError {
	oms_x := 1.0 - x
	oms_y := -y
	return Complex_dilog_xy_e(oms_x, oms_y, real_sp, imag_sp)
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Dilog(x float64) float64 {
	result := new(Result)
	status := Dilog_e(x, result)
	return EvalResult(result, status)
}
