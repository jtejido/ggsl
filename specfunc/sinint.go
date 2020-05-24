/* specfunc/sinint.c
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

/* based on SLATEC r9sifg.f, W. Fullerton */

/*
 series for f1   on the interval  2.00000e-02 to  6.25000e-02
                                        with weighted error   2.82e-17
                                         log weighted error  16.55
                               significant figures required  15.36
                                    decimal places required  17.20
*/
var (
	f1_data = []float64{
		-0.1191081969051363610,
		-0.0247823144996236248,
		0.0011910281453357821,
		-0.0000927027714388562,
		0.0000093373141568271,
		-0.0000011058287820557,
		0.0000001464772071460,
		-0.0000000210694496288,
		0.0000000032293492367,
		-0.0000000005206529618,
		0.0000000000874878885,
		-0.0000000000152176187,
		0.0000000000027257192,
		-0.0000000000005007053,
		0.0000000000000940241,
		-0.0000000000000180014,
		0.0000000000000035063,
		-0.0000000000000006935,
		0.0000000000000001391,
		-0.0000000000000000282,
	}

	f1_cs = &chebyshevSeries{
		f1_data,
		19,
		-1, 1,
		10,
	}
	/*

	   series for f2   on the interval  0.00000e+00 to  2.00000e-02
	                                          with weighted error   4.32e-17
	                                           log weighted error  16.36
	                                 significant figures required  14.75
	                                      decimal places required  17.10
	*/
	f2_data = []float64{
		-0.0348409253897013234,
		-0.0166842205677959686,
		0.0006752901241237738,
		-0.0000535066622544701,
		0.0000062693421779007,
		-0.0000009526638801991,
		0.0000001745629224251,
		-0.0000000368795403065,
		0.0000000087202677705,
		-0.0000000022601970392,
		0.0000000006324624977,
		-0.0000000001888911889,
		0.0000000000596774674,
		-0.0000000000198044313,
		0.0000000000068641396,
		-0.0000000000024731020,
		0.0000000000009226360,
		-0.0000000000003552364,
		0.0000000000001407606,
		-0.0000000000000572623,
		0.0000000000000238654,
		-0.0000000000000101714,
		0.0000000000000044259,
		-0.0000000000000019634,
		0.0000000000000008868,
		-0.0000000000000004074,
		0.0000000000000001901,
		-0.0000000000000000900,
		0.0000000000000000432,
	}

	f2_cs = &chebyshevSeries{
		f2_data,
		28,
		-1, 1,
		14,
	}
	/*

	   series for g1   on the interval  2.00000e-02 to  6.25000e-02
	                                          with weighted error   5.48e-17
	                                           log weighted error  16.26
	                                 significant figures required  15.47
	                                      decimal places required  16.92
	*/
	g1_data = []float64{
		-0.3040578798253495954,
		-0.0566890984597120588,
		0.0039046158173275644,
		-0.0003746075959202261,
		0.0000435431556559844,
		-0.0000057417294453025,
		0.0000008282552104503,
		-0.0000001278245892595,
		0.0000000207978352949,
		-0.0000000035313205922,
		0.0000000006210824236,
		-0.0000000001125215474,
		0.0000000000209088918,
		-0.0000000000039715832,
		0.0000000000007690431,
		-0.0000000000001514697,
		0.0000000000000302892,
		-0.0000000000000061400,
		0.0000000000000012601,
		-0.0000000000000002615,
		0.0000000000000000548,
	}

	g1_cs = &chebyshevSeries{
		g1_data,
		20,
		-1, 1,
		13,
	}
	/*

	   series for g2   on the interval  0.00000e+00 to  2.00000e-02
	                                          with weighted error   5.01e-17
	                                           log weighted error  16.30
	                                 significant figures required  15.12
	                                      decimal places required  17.07
	*/
	g2_data = []float64{
		-0.0967329367532432218,
		-0.0452077907957459871,
		0.0028190005352706523,
		-0.0002899167740759160,
		0.0000407444664601121,
		-0.0000071056382192354,
		0.0000014534723163019,
		-0.0000003364116512503,
		0.0000000859774367886,
		-0.0000000238437656302,
		0.0000000070831906340,
		-0.0000000022318068154,
		0.0000000007401087359,
		-0.0000000002567171162,
		0.0000000000926707021,
		-0.0000000000346693311,
		0.0000000000133950573,
		-0.0000000000053290754,
		0.0000000000021775312,
		-0.0000000000009118621,
		0.0000000000003905864,
		-0.0000000000001708459,
		0.0000000000000762015,
		-0.0000000000000346151,
		0.0000000000000159996,
		-0.0000000000000075213,
		0.0000000000000035970,
		-0.0000000000000017530,
		0.0000000000000008738,
		-0.0000000000000004487,
		0.0000000000000002397,
		-0.0000000000000001347,
		0.0000000000000000801,
		-0.0000000000000000501,
	}

	g2_cs = &chebyshevSeries{
		g2_data,
		33,
		-1, 1,
		20,
	}
)

/* x >= 4.0 */
func fg_asymp(x float64, f, g *Result) {
	/*
	   xbig = sqrt (1.0/r1mach(3))
	   xmaxf = exp (amin1(-alog(r1mach(1)), alog(r1mach(2))) - 0.01)
	   xmaxg = 1.0/sqrt(r1mach(1))
	   xbnd = sqrt(50.0)
	*/
	xbig := 1.0 / gsl.SqrtFloat64Eps
	xmaxf := 1.0 / gsl.MinFloat64
	xmaxg := 1.0 / gsl.SqrtMinFloat64
	xbnd := 7.07106781187
	x2 := x * x

	if x <= xbnd {
		result_c1, result_c2 := new(Result), new(Result)
		f1_cs.Evaluate((1.0/x2-0.04125)/0.02125, result_c1)
		g1_cs.Evaluate((1.0/x2-0.04125)/0.02125, result_c2)
		f.val = (1.0 + result_c1.val) / x
		g.val = (1.0 + result_c2.val) / x2
		f.err = result_c1.err/x + 2.0*gsl.Float64Eps*math.Abs(f.val)
		g.err = result_c2.err/x2 + 2.0*gsl.Float64Eps*math.Abs(g.val)
	} else if x <= xbig {
		result_c1, result_c2 := new(Result), new(Result)
		f2_cs.Evaluate(100.0/x2-1.0, result_c1)
		g2_cs.Evaluate(100.0/x2-1.0, result_c2)
		f.val = (1.0 + result_c1.val) / x
		g.val = (1.0 + result_c2.val) / x2
		f.err = result_c1.err/x + 2.0*gsl.Float64Eps*math.Abs(f.val)
		g.err = result_c2.err/x2 + 2.0*gsl.Float64Eps*math.Abs(g.val)
	} else {
		if x < xmaxf {
			f.val = 1.0 / x
		}

		if x < xmaxg {
			g.val = 1.0 / x2
		}

		f.err = 2.0 * gsl.Float64Eps * math.Abs(f.val)
		g.err = 2.0 * gsl.Float64Eps * math.Abs(g.val)
	}

	return
}

/* based on SLATEC si.f, W. Fullerton

series for si   on the interval  0.00000e+00 to  1.60000e+01
                                       with weighted error   1.22e-17
                                        log weighted error  16.91
                              significant figures required  16.37
                                   decimal places required  17.45
*/
var (
	si_data = []float64{
		-0.1315646598184841929,
		-0.2776578526973601892,
		0.0354414054866659180,
		-0.0025631631447933978,
		0.0001162365390497009,
		-0.0000035904327241606,
		0.0000000802342123706,
		-0.0000000013562997693,
		0.0000000000179440722,
		-0.0000000000001908387,
		0.0000000000000016670,
		-0.0000000000000000122,
	}

	si_cs = &chebyshevSeries{
		si_data,
		11,
		-1, 1,
		9,
	}
	/*
	   series for ci   on the interval  0.00000e+00 to  1.60000e+01
	                                          with weighted error   1.94e-18
	                                           log weighted error  17.71
	                                 significant figures required  17.74
	                                      decimal places required  18.27
	*/
	ci_data = []float64{
		-0.34004281856055363156,
		-1.03302166401177456807,
		0.19388222659917082877,
		-0.01918260436019865894,
		0.00110789252584784967,
		-0.00004157234558247209,
		0.00000109278524300229,
		-0.00000002123285954183,
		0.00000000031733482164,
		-0.00000000000376141548,
		0.00000000000003622653,
		-0.00000000000000028912,
		0.00000000000000000194,
	}

	ci_cs = &chebyshevSeries{
		ci_data,
		12,
		-1, 1,
		9,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Si_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	if ax < gsl.SqrtFloat64Eps {
		result.val = x
		result.err = 0.0
		return nil
	} else if ax <= 4.0 {
		result_c := new(Result)
		si_cs.Evaluate((x*x-8.0)*0.125, result_c)
		result.val = x * (0.75 + result_c.val)
		result.err = ax * result_c.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		/* Note there is no loss of precision
		 * here bcause of the leading constant.
		 */
		f, g := new(Result), new(Result)
		fg_asymp(ax, f, g)
		result.val = 0.5*gsl.Pi - f.val*math.Cos(ax) - g.val*math.Sin(ax)
		result.err = f.err + g.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		if x < 0.0 {
			result.val = -result.val
		}
		return nil
	}
}

func Ci_e(x float64, result *Result) err.GSLError {
	if x <= 0.0 {
		return DomainError(result)
	} else if x <= 4.0 {
		lx := math.Log(x)
		y := (x*x - 8.0) * 0.125
		result_c := new(Result)
		ci_cs.Evaluate(y, result_c)
		result.val = lx - 0.5 + result_c.val
		result.err = 2.0*gsl.Float64Eps*(math.Abs(lx)+0.5) + result_c.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		sin_result, cos_result := new(Result), new(Result)
		stat_sin := Sin_e(x, sin_result)
		stat_cos := Cos_e(x, cos_result)
		f, g := new(Result), new(Result)
		fg_asymp(x, f, g)
		result.val = f.val*sin_result.val - g.val*cos_result.val
		result.err = math.Abs(f.err * sin_result.val)
		result.err += math.Abs(g.err * cos_result.val)
		result.err += math.Abs(f.val * sin_result.err)
		result.err += math.Abs(g.val * cos_result.err)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_sin, stat_cos)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Si(x float64) float64 {
	result := new(Result)
	status := Si_e(x, result)
	return EvalResult(result, status)
}

func Ci(x float64) float64 {
	result := new(Result)
	status := Ci_e(x, result)
	return EvalResult(result, status)
}
