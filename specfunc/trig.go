/* specfunc/trig.c
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

var (
	sinc_data = []float64{
		1.133648177811747875422,
		-0.532677564732557348781,
		-0.068293048346633177859,
		0.033403684226353715020,
		0.001485679893925747818,
		-0.000734421305768455295,
		-0.000016837282388837229,
		0.000008359950146618018,
		0.000000117382095601192,
		-0.000000058413665922724,
		-0.000000000554763755743,
		0.000000000276434190426,
		0.000000000001895374892,
		-0.000000000000945237101,
		-0.000000000000004900690,
		0.000000000000002445383,
		0.000000000000000009925,
	}

	sinc_cs = &chebyshevSeries{
		sinc_data,
		16,
		-1, 1,
		10,
	}

	sin_data = []float64{
		-0.3295190160663511504173,
		0.0025374284671667991990,
		0.0006261928782647355874,
		-4.6495547521854042157541e-06,
		-5.6917531549379706526677e-07,
		3.7283335140973803627866e-09,
		3.0267376484747473727186e-10,
		-1.7400875016436622322022e-12,
		-1.0554678305790849834462e-13,
		5.3701981409132410797062e-16,
		2.5984137983099020336115e-17,
		-1.1821555255364833468288e-19,
	}

	sins = &chebyshevSeries{
		sin_data,
		11,
		-1, 1,
		11,
	}

	cos_data = []float64{
		0.165391825637921473505668118136,
		-0.00084852883845000173671196530195,
		-0.000210086507222940730213625768083,
		1.16582269619760204299639757584e-6,
		1.43319375856259870334412701165e-7,
		-7.4770883429007141617951330184e-10,
		-6.0969994944584252706997438007e-11,
		2.90748249201909353949854872638e-13,
		1.77126739876261435667156490461e-14,
		-7.6896421502815579078577263149e-17,
		-3.7363121133079412079201377318e-18,
	}

	coss = &chebyshevSeries{
		cos_data,
		10,
		-1, 1,
		10,
	}
)

/* sinh(x) series
 * double-precision for |x| < 1.0
 */
func sinh_series(x float64, result *float64) err.GSLError {
	y := x * x
	c0 := 1.0 / 6.0
	c1 := 1.0 / 120.0
	c2 := 1.0 / 5040.0
	c3 := 1.0 / 362880.0
	c4 := 1.0 / 39916800.0
	c5 := 1.0 / 6227020800.0
	c6 := 1.0 / 1307674368000.0
	c7 := 1.0 / 355687428096000.0

	*result = x * (1.0 + y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))))
	return nil
}

/* cosh(x)-1 series
 * double-precision for |x| < 1.0
 */
func cosh_m1_series(x float64, result *float64) err.GSLError {
	y := x * x
	c0 := 0.5
	c1 := 1.0 / 24.0
	c2 := 1.0 / 720.0
	c3 := 1.0 / 40320.0
	c4 := 1.0 / 3628800.0
	c5 := 1.0 / 479001600.0
	c6 := 1.0 / 87178291200.0
	c7 := 1.0 / 20922789888000.0
	c8 := 1.0 / 6402373705728000.0

	*result = y * (c0 + y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*c8))))))))
	return nil
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/* I would have prefered just using the library sin() function.
 * But after some experimentation I decided that there was
 * no good way to understand the error; library sin() is just a black box.
 * So we have to roll our own.
 */
func Sin_e(x float64, result *Result) err.GSLError {
	const (
		P1 = 7.85398125648498535156e-1
		P2 = 3.77489470793079817668e-8
		P3 = 2.69515142907905952645e-15
	)

	sgn_x := gsl.Sign(x)
	abs_x := math.Abs(x)

	if abs_x < gsl.Root4Float64Eps {
		x2 := x * x
		result.val = x * (1.0 - x2/6.0)
		result.err = math.Abs(x * x2 * x2 / 100.0)
		return nil
	} else {
		var stat_cs err.GSLError
		sgn_result := sgn_x
		y := math.Floor(abs_x / (0.25 * gsl.Pi))
		octant := int(y - math.Ldexp(math.Floor(math.Ldexp(y, -3)), 3))

		if gsl.IsOdd(octant) {
			octant += 1
			octant &= 07
			y += 1.0
		}

		if octant > 3 {
			octant -= 4
			sgn_result = -sgn_result
		}

		z := ((abs_x - y*P1) - y*P2) - y*P3

		if octant == 0 {
			sin_cs_result := new(Result)
			t := 8.0*math.Abs(z)/gsl.Pi - 1.0
			stat_cs = sins.Evaluate(t, sin_cs_result)
			result.val = z * (1.0 + z*z*sin_cs_result.val)
		} else {
			cos_cs_result := new(Result)
			t := 8.0*math.Abs(z)/gsl.Pi - 1.0
			stat_cs = coss.Evaluate(t, cos_cs_result)
			result.val = 1.0 - 0.5*z*z*(1.0-z*z*cos_cs_result.val)
		}

		result.val *= sgn_result

		if abs_x > 1.0/gsl.Float64Eps {
			result.err = math.Abs(result.val)
		} else if abs_x > 100.0/gsl.SqrtFloat64Eps {
			result.err = 2.0 * abs_x * gsl.Float64Eps * math.Abs(result.val)
		} else if abs_x > 0.1/gsl.SqrtFloat64Eps {
			result.err = 2.0 * gsl.SqrtFloat64Eps * math.Abs(result.val)
		} else {
			result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		}

		return stat_cs
	}

}

func Cos_e(x float64, result *Result) err.GSLError {
	const (
		P1 = 7.85398125648498535156e-1
		P2 = 3.77489470793079817668e-8
		P3 = 2.69515142907905952645e-15
	)

	abs_x := math.Abs(x)

	if abs_x < gsl.Root4Float64Eps {
		x2 := x * x
		result.val = 1.0 - 0.5*x2
		result.err = math.Abs(x2 * x2 / 12.0)
		return nil
	} else {
		sgn_result := 1.0
		y := math.Floor(abs_x / (0.25 * gsl.Pi))
		octant := int(y - math.Ldexp(math.Floor(math.Ldexp(y, -3)), 3))
		var stat_cs err.GSLError
		var z float64

		if gsl.IsOdd(octant) {
			octant += 1
			octant &= 07
			y += 1.0
		}

		if octant > 3 {
			octant -= 4
			sgn_result = -sgn_result
		}

		if octant > 1 {
			sgn_result = -sgn_result
		}

		z = ((abs_x - y*P1) - y*P2) - y*P3

		if octant == 0 {
			cos_cs_result := new(Result)
			t := 8.0*math.Abs(z)/gsl.Pi - 1.0
			stat_cs = coss.Evaluate(t, cos_cs_result)
			result.val = 1.0 - 0.5*z*z*(1.0-z*z*cos_cs_result.val)
		} else { /* octant == 2 */
			sin_cs_result := new(Result)
			t := 8.0*math.Abs(z)/gsl.Pi - 1.0
			stat_cs = sins.Evaluate(t, sin_cs_result)
			result.val = z * (1.0 + z*z*sin_cs_result.val)
		}

		result.val *= sgn_result

		if abs_x > 1.0/gsl.Float64Eps {
			result.err = math.Abs(result.val)
		} else if abs_x > 100.0/gsl.SqrtFloat64Eps {
			result.err = 2.0 * abs_x * gsl.Float64Eps * math.Abs(result.val)
		} else if abs_x > 0.1/gsl.SqrtFloat64Eps {
			result.err = 2.0 * gsl.SqrtFloat64Eps * math.Abs(result.val)
		} else {
			result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		}

		return stat_cs
	}

}

func Hypot_e(x, y float64, result *Result) err.GSLError {
	if x == 0.0 && y == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		a := math.Abs(x)
		b := math.Abs(y)
		min := gsl.Min(a, b)
		max := gsl.Max(a, b)
		rat := min / max
		root_term := math.Sqrt(1.0 + rat*rat)

		if max < gsl.MaxFloat64/root_term {
			result.val = max * root_term
			result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return nil
		} else {
			return OverflowError(result)
		}
	}
}

func Complex_sin_e(zr, zi float64, szr, szi *Result) err.GSLError {
	if math.Abs(zi) < 1.0 {
		var ch_m1, sh float64
		sinh_series(zi, &sh)
		cosh_m1_series(zi, &ch_m1)
		szr.val = math.Sin(zr) * (ch_m1 + 1.0)
		szi.val = math.Cos(zr) * sh
		szr.err = 2.0 * gsl.Float64Eps * math.Abs(szr.val)
		szi.err = 2.0 * gsl.Float64Eps * math.Abs(szi.val)
		return nil
	} else if math.Abs(zi) < gsl.LnMaxFloat64 {
		ex := math.Exp(zi)
		ch := 0.5 * (ex + 1.0/ex)
		sh := 0.5 * (ex - 1.0/ex)
		szr.val = math.Sin(zr) * ch
		szi.val = math.Cos(zr) * sh
		szr.err = 2.0 * gsl.Float64Eps * math.Abs(szr.val)
		szi.err = 2.0 * gsl.Float64Eps * math.Abs(szi.val)
		return nil
	}

	return OverflowError(szr, szi)
}

func Complex_cos_e(zr, zi float64, czr, czi *Result) err.GSLError {

	if math.Abs(zi) < 1.0 {
		var ch_m1, sh float64
		sinh_series(zi, &sh)
		cosh_m1_series(zi, &ch_m1)
		czr.val = math.Cos(zr) * (ch_m1 + 1.0)
		czi.val = -math.Sin(zr) * sh
		czr.err = 2.0 * gsl.Float64Eps * math.Abs(czr.val)
		czi.err = 2.0 * gsl.Float64Eps * math.Abs(czi.val)
		return nil
	} else if math.Abs(zi) < gsl.LnMaxFloat64 {
		ex := math.Exp(zi)
		ch := 0.5 * (ex + 1.0/ex)
		sh := 0.5 * (ex - 1.0/ex)
		czr.val = math.Cos(zr) * ch
		czi.val = -math.Sin(zr) * sh
		czr.err = 2.0 * gsl.Float64Eps * math.Abs(czr.val)
		czi.err = 2.0 * gsl.Float64Eps * math.Abs(czi.val)
		return nil

	}

	return OverflowError(czr, czi)
}

func Complex_logsin_e(zr, zi float64, lszr, lszi *Result) err.GSLError {

	if zi > 60.0 {
		lszr.val = -math.Ln2 + zi
		lszi.val = 0.5*gsl.Pi - zr
		lszr.err = 2.0 * gsl.Float64Eps * math.Abs(lszr.val)
		lszi.err = 2.0 * gsl.Float64Eps * math.Abs(lszi.val)
	} else if zi < -60.0 {
		lszr.val = -math.Ln2 - zi
		lszi.val = -0.5*gsl.Pi + zr
		lszr.err = 2.0 * gsl.Float64Eps * math.Abs(lszr.val)
		lszi.err = 2.0 * gsl.Float64Eps * math.Abs(lszi.val)
	} else {
		sin_r, sin_i := new(Result), new(Result)
		Complex_sin_e(zr, zi, sin_r, sin_i)
		status := Complex_log_e(sin_r.val, sin_i.val, lszr, lszi)
		if status != nil {
			if status_err, ok := status.(err.GSLError); ok {
				if status_err.Status() == err.EDOM {
					return DomainError(lszr, lszi)
				}
			}

		}

	}
	return Angle_restrict_symm_e(&lszi.val)
}

func Lnsinh_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if math.Abs(x) < 1.0 {
		var eps float64
		sinh_series(x, &eps)
		result.val = math.Log(eps)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < -0.5*gsl.LnFloat64Eps {
		result.val = x + math.Log(0.5*(1.0-math.Exp(-2.0*x)))
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = -math.Ln2 + x
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Lncosh_e(x float64, result *Result) err.GSLError {
	if math.Abs(x) < 1.0 {
		var eps float64
		cosh_m1_series(x, &eps)
		return Log_1plusx_e(eps, result)
	} else if math.Abs(x) < -0.5*gsl.LnFloat64Eps {
		result.val = math.Abs(x) + math.Log(0.5*(1.0+math.Exp(-2.0*math.Abs(x))))
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result.val = -math.Ln2 + math.Abs(x)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Polar_to_rect(r, theta float64, x, y *Result) err.GSLError {
	t := theta
	status := Angle_restrict_symm_e(&t)
	c := math.Cos(t)
	s := math.Sin(t)
	x.val = r * math.Cos(t)
	y.val = r * math.Sin(t)
	x.err = r * math.Abs(s*gsl.Float64Eps*t)
	x.err += 2.0 * gsl.Float64Eps * math.Abs(x.val)
	y.err = r * math.Abs(c*gsl.Float64Eps*t)
	y.err += 2.0 * gsl.Float64Eps * math.Abs(y.val)
	return status
}

func Rect_to_polar(x, y float64, r, theta *Result) err.GSLError {
	stat_h := Hypot_e(x, y, r)
	if r.val > 0.0 {
		theta.val = math.Atan2(y, x)
		theta.err = 2.0 * gsl.Float64Eps * math.Abs(theta.val)
		return stat_h
	}

	return DomainError(theta)

}

// These routines force the angle theta to lie in the range (-pi,pi].
// Note that the mathematical value of pi is slightly greater than gsl.Pi, so the machine numbers gsl.Pi and -gsl.Pi are included in the range.
func Angle_restrict_symm_err_e(theta float64, result *Result) err.GSLError {
	/* synthetic extended precision constants */
	const (
		P1    = 4 * 7.8539812564849853515625e-01
		P2    = 4 * 3.7748947079307981766760e-08
		P3    = 4 * 2.6951514290790594840552e-15
		TwoPi = 2 * (P1 + P2 + P3)
	)

	y := gsl.Sign(theta) * 2 * math.Floor(math.Abs(theta)/TwoPi)
	r := ((theta - y*P1) - y*P2) - y*P3

	if r > gsl.Pi {
		r = (((r - 2*P1) - 2*P2) - 2*P3) /* r-TwoPi */
	} else if r < -gsl.Pi {
		r = (((r + 2*P1) + 2*P2) + 2*P3) /* r+TwoPi */
	}

	result.val = r

	if math.Abs(theta) > 0.0625/gsl.Float64Eps {
		result.val = math.NaN()
		result.err = math.NaN()
		return err.ERROR("error", err.ELOSS)
	} else if math.Abs(theta) > 0.0625/gsl.SqrtFloat64Eps {
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val-theta)
		return nil
	} else {
		delta := math.Abs(result.val - theta)
		result.err = 2.0 * gsl.Float64Eps
		if delta < gsl.Pi {
			result.err *= delta
		} else {
			result.err *= gsl.Pi
		}
		return nil
	}

}

func Angle_restrict_pos_err_e(theta float64, result *Result) err.GSLError {
	/* synthetic extended precision constants */
	P1 := 4 * 7.85398125648498535156e-01
	P2 := 4 * 3.77489470793079817668e-08
	P3 := 4 * 2.69515142907905952645e-15
	TwoPi := 2 * (P1 + P2 + P3)

	y := 2 * math.Floor(theta/TwoPi)

	r := ((theta - y*P1) - y*P2) - y*P3

	if r > TwoPi { /* r-TwoPi */
		r = (((r - 2*P1) - 2*P2) - 2*P3)
	} else if r < 0 { /* may happen due to FP rounding */
		r = (((r + 2*P1) + 2*P2) + 2*P3) /* r+TwoPi */
	}

	result.val = r

	if math.Abs(theta) > 0.0625/gsl.Float64Eps {
		result.val = math.NaN()
		result.err = math.Abs(result.val)
		return err.ERROR("error", err.ELOSS)
	} else if math.Abs(theta) > 0.0625/gsl.SqrtFloat64Eps {
		result.err = gsl.Float64Eps * math.Abs(result.val-theta)
		return nil
	} else {
		delta := math.Abs(result.val - theta)
		result.err = 2.0 * gsl.Float64Eps
		if delta < gsl.Pi {
			result.err *= delta
		} else {
			result.err *= gsl.Pi
		}
		return nil
	}
}

func Angle_restrict_symm_e(theta *float64) err.GSLError {
	r := new(Result)
	stat := Angle_restrict_symm_err_e(*theta, r)
	*theta = r.val
	return stat
}

func Angle_restrict_pos_e(theta *float64) err.GSLError {
	r := new(Result)
	stat := Angle_restrict_pos_err_e(*theta, r)
	*theta = r.val
	return stat
}

func Sin_err_e(x, dx float64, result *Result) err.GSLError {
	stat_s := Sin_e(x, result)
	result.err += math.Abs(math.Cos(x) * dx)
	result.err += gsl.Float64Eps * math.Abs(result.val)
	return stat_s
}

func Cos_err_e(x, dx float64, result *Result) err.GSLError {
	stat_c := Cos_e(x, result)
	result.err += math.Abs(math.Sin(x) * dx)
	result.err += gsl.Float64Eps * math.Abs(result.val)
	return stat_c
}

func Sinc_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	if ax < 0.8 {
		/* Do not go to the limit of the fit since
		 * there is a zero there and the Chebyshev
		 * accuracy will go to zero.
		 */
		return sinc_cs.Evaluate(2.0*ax-1.0, result)
	} else if ax < 100.0 {
		/* Small arguments are no problem.
		 * We trust the library sin() to
		 * roughly machine precision.
		 */
		result.val = math.Sin(gsl.Pi*ax) / (gsl.Pi * ax)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		/* Large arguments must be handled separately.
		 */
		r := gsl.Pi * ax
		s := new(Result)
		stat_s := Sin_e(r, s)
		result.val = s.val / r
		result.err = s.err/r + 2.0*gsl.Float64Eps*math.Abs(result.val)
		return stat_s
	}

}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Sin(x float64) float64 {
	result := new(Result)
	status := Sin_e(x, result)
	return EvalResult(result, status)
}

func Cos(x float64) float64 {
	result := new(Result)
	status := Cos_e(x, result)
	return EvalResult(result, status)
}

func Hypot(x, y float64) float64 {
	result := new(Result)
	status := Hypot_e(x, y, result)
	return EvalResult(result, status)
}

func Lnsinh(x float64) float64 {
	result := new(Result)
	status := Lnsinh_e(x, result)
	return EvalResult(result, status)
}

func Lncosh(x float64) float64 {
	result := new(Result)
	status := Lncosh_e(x, result)
	return EvalResult(result, status)
}

func Angle_restrict_symm(theta float64) float64 {
	var result float64 = theta
	status := Angle_restrict_symm_e(&result)
	return EvalFloat64(result, status)
}

func Angle_restrict_pos(theta float64) float64 {
	var result float64 = theta
	status := Angle_restrict_pos_e(&result)
	return EvalFloat64(result, status)
}

func Sinc(x float64) float64 {
	result := new(Result)
	status := Sinc_e(x, result)
	return EvalResult(result, status)
}
