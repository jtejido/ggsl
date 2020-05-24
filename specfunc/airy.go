/* specfunc/airy.c
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
	am21d = []float64{
		0.0065809191761485,
		0.0023675984685722,
		0.0001324741670371,
		0.0000157600904043,
		0.0000027529702663,
		0.0000006102679017,
		0.0000001595088468,
		0.0000000471033947,
		0.0000000152933871,
		0.0000000053590722,
		0.0000000020000910,
		0.0000000007872292,
		0.0000000003243103,
		0.0000000001390106,
		0.0000000000617011,
		0.0000000000282491,
		0.0000000000132979,
		0.0000000000064188,
		0.0000000000031697,
		0.0000000000015981,
		0.0000000000008213,
		0.0000000000004296,
		0.0000000000002284,
		0.0000000000001232,
		0.0000000000000675,
		0.0000000000000374,
		0.0000000000000210,
		0.0000000000000119,
		0.0000000000000068,
		0.0000000000000039,
		0.0000000000000023,
		0.0000000000000013,
		0.0000000000000008,
		0.0000000000000005,
		0.0000000000000003,
		0.0000000000000001,
		0.0000000000000001,
	}

	am21 = &chebyshevSeries{
		am21d,
		36,
		-1, 1,
		20,
	}

	ath1d = []float64{
		-0.07125837815669365,
		-0.00590471979831451,
		-0.00012114544069499,
		-0.00000988608542270,
		-0.00000138084097352,
		-0.00000026142640172,
		-0.00000006050432589,
		-0.00000001618436223,
		-0.00000000483464911,
		-0.00000000157655272,
		-0.00000000055231518,
		-0.00000000020545441,
		-0.00000000008043412,
		-0.00000000003291252,
		-0.00000000001399875,
		-0.00000000000616151,
		-0.00000000000279614,
		-0.00000000000130428,
		-0.00000000000062373,
		-0.00000000000030512,
		-0.00000000000015239,
		-0.00000000000007758,
		-0.00000000000004020,
		-0.00000000000002117,
		-0.00000000000001132,
		-0.00000000000000614,
		-0.00000000000000337,
		-0.00000000000000188,
		-0.00000000000000105,
		-0.00000000000000060,
		-0.00000000000000034,
		-0.00000000000000020,
		-0.00000000000000011,
		-0.00000000000000007,
		-0.00000000000000004,
		-0.00000000000000002,
	}

	ath1 = &chebyshevSeries{
		ath1d,
		35,
		-1, 1,
		15,
	}

	am22d = []float64{
		-0.01562844480625341,
		0.00778336445239681,
		0.00086705777047718,
		0.00015696627315611,
		0.00003563962571432,
		0.00000924598335425,
		0.00000262110161850,
		0.00000079188221651,
		0.00000025104152792,
		0.00000008265223206,
		0.00000002805711662,
		0.00000000976821090,
		0.00000000347407923,
		0.00000000125828132,
		0.00000000046298826,
		0.00000000017272825,
		0.00000000006523192,
		0.00000000002490471,
		0.00000000000960156,
		0.00000000000373448,
		0.00000000000146417,
		0.00000000000057826,
		0.00000000000022991,
		0.00000000000009197,
		0.00000000000003700,
		0.00000000000001496,
		0.00000000000000608,
		0.00000000000000248,
		0.00000000000000101,
		0.00000000000000041,
		0.00000000000000017,
		0.00000000000000007,
		0.00000000000000002,
	}

	am22 = &chebyshevSeries{
		am22d,
		32,
		-1, 1,
		15,
	}

	ath2d = []float64{
		0.00440527345871877,
		-0.03042919452318455,
		-0.00138565328377179,
		-0.00018044439089549,
		-0.00003380847108327,
		-0.00000767818353522,
		-0.00000196783944371,
		-0.00000054837271158,
		-0.00000016254615505,
		-0.00000005053049981,
		-0.00000001631580701,
		-0.00000000543420411,
		-0.00000000185739855,
		-0.00000000064895120,
		-0.00000000023105948,
		-0.00000000008363282,
		-0.00000000003071196,
		-0.00000000001142367,
		-0.00000000000429811,
		-0.00000000000163389,
		-0.00000000000062693,
		-0.00000000000024260,
		-0.00000000000009461,
		-0.00000000000003716,
		-0.00000000000001469,
		-0.00000000000000584,
		-0.00000000000000233,
		-0.00000000000000093,
		-0.00000000000000037,
		-0.00000000000000015,
		-0.00000000000000006,
		-0.00000000000000002,
	}

	ath2 = &chebyshevSeries{
		ath2d,
		31,
		-1, 1,
		16,
	}
)

func airy_mod_phase(x float64, mode gsl.MODE_T, mod, phase *Result) err.GSLError {
	result_m, result_p := new(Result), new(Result)
	var m, p, sqx float64
	if x < -2.0 {
		z := 16.0/(x*x*x) + 1.0
		am21.EvaluateMode(z, mode, result_m)
		ath1.EvaluateMode(z, mode, result_p)
	} else if x <= -1.0 {
		z := (16.0/(x*x*x) + 9.0) / 7.0
		am22.EvaluateMode(z, mode, result_m)
		ath2.EvaluateMode(z, mode, result_p)
	} else {
		mod.val = 0.0
		mod.err = 0.0
		phase.val = 0.0
		phase.err = 0.0
		return err.ERROR("x is greater than 1.0", err.EDOM)
	}

	m = 0.3125 + result_m.val
	p = -0.625 + result_p.val

	sqx = math.Sqrt(-x)

	mod.val = math.Sqrt(m / sqx)
	mod.err = math.Abs(mod.val) * (gsl.Float64Eps + math.Abs(result_m.err/result_m.val))
	phase.val = gsl.PiOver4 - x*sqx*p
	phase.err = math.Abs(phase.val) * (gsl.Float64Eps + math.Abs(result_p.err/result_p.val))

	return nil
}

var (
	aifd = []float64{
		-0.03797135849666999750,
		0.05919188853726363857,
		0.00098629280577279975,
		0.00000684884381907656,
		0.00000002594202596219,
		0.00000000006176612774,
		0.00000000000010092454,
		0.00000000000000012014,
		0.00000000000000000010,
	}

	aibf = &chebyshevSeries{
		aifd,
		8,
		-1, 1,
		8,
	}

	aigd = []float64{
		0.01815236558116127,
		0.02157256316601076,
		0.00025678356987483,
		0.00000142652141197,
		0.00000000457211492,
		0.00000000000952517,
		0.00000000000001392,
		0.00000000000000001,
	}

	aibg = &chebyshevSeries{
		aigd,
		7,
		-1, 1,
		7,
	}

	bifd = []float64{
		-0.01673021647198664948,
		0.10252335834249445610,
		0.00170830925073815165,
		0.00001186254546774468,
		0.00000004493290701779,
		0.00000000010698207143,
		0.00000000000017480643,
		0.00000000000000020810,
		0.00000000000000000018,
	}

	aibif = &chebyshevSeries{
		bifd,
		8,
		-1, 1,
		8,
	}

	bigd = []float64{
		0.02246622324857452,
		0.03736477545301955,
		0.00044476218957212,
		0.00000247080756363,
		0.00000000791913533,
		0.00000000001649807,
		0.00000000000002411,
		0.00000000000000002,
	}

	aibig = &chebyshevSeries{
		bigd,
		7,
		-1, 1,
		7,
	}

	bif2d = []float64{
		0.0998457269381604100,
		0.4786249778630055380,
		0.0251552119604330118,
		0.0005820693885232645,
		0.0000074997659644377,
		0.0000000613460287034,
		0.0000000003462753885,
		0.0000000000014288910,
		0.0000000000000044962,
		0.0000000000000000111,
	}

	aibif2 = &chebyshevSeries{
		bif2d,
		9,
		-1, 1,
		9,
	}

	big2d = []float64{
		0.033305662145514340,
		0.161309215123197068,
		0.0063190073096134286,
		0.0001187904568162517,
		0.0000013045345886200,
		0.0000000093741259955,
		0.0000000000474580188,
		0.0000000000001783107,
		0.0000000000000005167,
		0.0000000000000000011,
	}

	aibig2 = &chebyshevSeries{
		big2d,
		9,
		-1, 1,
		9,
	}

	aipd = []float64{
		-0.0187519297793867540198,
		-0.0091443848250055004725,
		0.0009010457337825074652,
		-0.0001394184127221491507,
		0.0000273815815785209370,
		-0.0000062750421119959424,
		0.0000016064844184831521,
		-0.0000004476392158510354,
		0.0000001334635874651668,
		-0.0000000420735334263215,
		0.0000000139021990246364,
		-0.0000000047831848068048,
		0.0000000017047897907465,
		-0.0000000006268389576018,
		0.0000000002369824276612,
		-0.0000000000918641139267,
		0.0000000000364278543037,
		-0.0000000000147475551725,
		0.0000000000060851006556,
		-0.0000000000025552772234,
		0.0000000000010906187250,
		-0.0000000000004725870319,
		0.0000000000002076969064,
		-0.0000000000000924976214,
		0.0000000000000417096723,
		-0.0000000000000190299093,
		0.0000000000000087790676,
		-0.0000000000000040927557,
		0.0000000000000019271068,
		-0.0000000000000009160199,
		0.0000000000000004393567,
		-0.0000000000000002125503,
		0.0000000000000001036735,
		-0.0000000000000000509642,
		0.0000000000000000252377,
		-0.0000000000000000125793,
	}

	aip = &chebyshevSeries{
		aipd,
		35,
		-1, 1,
		17,
	}

	bipd = []float64{
		-0.08322047477943447,
		0.01146118927371174,
		0.00042896440718911,
		-0.00014906639379950,
		-0.00001307659726787,
		0.00000632759839610,
		-0.00000042226696982,
		-0.00000019147186298,
		0.00000006453106284,
		-0.00000000784485467,
		-0.00000000096077216,
		0.00000000070004713,
		-0.00000000017731789,
		0.00000000002272089,
		0.00000000000165404,
		-0.00000000000185171,
		0.00000000000059576,
		-0.00000000000012194,
		0.00000000000001334,
		0.00000000000000172,
		-0.00000000000000145,
		0.00000000000000049,
		-0.00000000000000011,
		0.00000000000000001,
	}

	bip = &chebyshevSeries{
		bipd,
		23,
		-1, 1,
		14,
	}

	bip2d = []float64{
		-0.113596737585988679,
		0.0041381473947881595,
		0.0001353470622119332,
		0.0000104273166530153,
		0.0000013474954767849,
		0.0000001696537405438,
		-0.0000000100965008656,
		-0.0000000167291194937,
		-0.0000000045815364485,
		0.0000000003736681366,
		0.0000000005766930320,
		0.0000000000621812650,
		-0.0000000000632941202,
		-0.0000000000149150479,
		0.0000000000078896213,
		0.0000000000024960513,
		-0.0000000000012130075,
		-0.0000000000003740493,
		0.0000000000002237727,
		0.0000000000000474902,
		-0.0000000000000452616,
		-0.0000000000000030172,
		0.0000000000000091058,
		-0.0000000000000009814,
		-0.0000000000000016429,
		0.0000000000000005533,
		0.0000000000000002175,
		-0.0000000000000001737,
		-0.0000000000000000010,
	}

	aibip2 = &chebyshevSeries{
		bip2d,
		28,
		-1, 1,
		10,
	}
)

func airy_aie(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	sqx := math.Sqrt(x)
	z := 2.0/(x*sqx) - 1.0
	y := math.Sqrt(sqx)
	result_c := new(Result)
	aip.EvaluateMode(z, mode, result_c)
	result.val = (0.28125 + result_c.val) / y
	result.err = result_c.err/y + gsl.Float64Eps*math.Abs(result.val)
	return nil
}

func airy_bie(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	ATR := 8.7506905708484345
	BTR := -2.0938363213560543

	if x < 4.0 {
		sqx := math.Sqrt(x)
		z := ATR/(x*sqx) + BTR
		y := math.Sqrt(sqx)
		result_c := new(Result)
		bip.EvaluateMode(z, mode, result_c)
		result.val = (0.625 + result_c.val) / y
		result.err = result_c.err/y + gsl.Float64Eps*math.Abs(result.val)
	} else {
		sqx := math.Sqrt(x)
		z := 16.0/(x*sqx) - 1.0
		y := math.Sqrt(sqx)
		result_c := new(Result)
		aibip2.EvaluateMode(z, mode, result_c)
		result.val = (0.625 + result_c.val) / y
		result.err = result_c.err/y + gsl.Float64Eps*math.Abs(result.val)
	}

	return nil

}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

// Computes the Airy function Ai(x).
func Airy_Ai_e(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	if x < -1.0 {
		mod, theta, cos_result := new(Result), new(Result), new(Result)
		stat_mp := airy_mod_phase(x, mode, mod, theta)
		stat_cos := Cos_err_e(theta.val, theta.err, cos_result)
		result.val = mod.val * cos_result.val
		result.err = math.Abs(mod.val*cos_result.err) + math.Abs(cos_result.val*mod.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_mp, stat_cos)
	} else if x <= 1.0 {
		z := x * x * x
		result_c0, result_c1 := new(Result), new(Result)
		aibf.EvaluateMode(z, mode, result_c0)
		aibg.EvaluateMode(z, mode, result_c1)
		result.val = 0.375 + (result_c0.val - x*(0.25+result_c1.val))
		result.err = result_c0.err + math.Abs(x*result_c1.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	}

	x32 := x * math.Sqrt(x)
	s := math.Exp(-2.0 * x32 / 3.0)
	result_aie := new(Result)
	stat_aie := airy_aie(x, mode, result_aie)
	result.val = result_aie.val * s
	result.err = result_aie.err*s + result.val*x32*gsl.Float64Eps
	result.err += gsl.Float64Eps * math.Abs(result.val)
	CheckUnderflow(result)
	return stat_aie
}

// This computes a scaled version of the Airy function S_A(x) Ai(x).
// For x > 0 the scaling factor S_A(x) is exp(+(2/3) x^{3/2}), and is 1 for x < 0.
func Airy_Ai_scaled_e(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	if x < -1.0 {
		mod, theta, cos_result := new(Result), new(Result), new(Result)
		stat_mp := airy_mod_phase(x, mode, mod, theta)
		stat_cos := Cos_err_e(theta.val, theta.err, cos_result)
		result.val = mod.val * cos_result.val
		result.err = math.Abs(mod.val*cos_result.err) + math.Abs(cos_result.val*mod.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_mp, stat_cos)
	} else if x <= 1.0 {
		z := x * x * x
		result_c0, result_c1 := new(Result), new(Result)
		aibf.EvaluateMode(z, mode, result_c0)
		aibg.EvaluateMode(z, mode, result_c1)
		result.val = 0.375 + (result_c0.val - x*(0.25+result_c1.val))
		result.err = result_c0.err + math.Abs(x*result_c1.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)

		if x > 0.0 {
			scale := math.Exp(2.0 / 3.0 * math.Sqrt(z))
			result.val *= scale
			result.err *= scale
		}

		return nil
	} else {
		return airy_aie(x, mode, result)
	}

}

// Computes the Airy function Bi(x).
func Airy_Bi_e(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	if x < -1.0 {
		mod, theta, sin_result := new(Result), new(Result), new(Result)
		stat_mp := airy_mod_phase(x, mode, mod, theta)
		stat_sin := Sin_err_e(theta.val, theta.err, sin_result)
		result.val = mod.val * sin_result.val
		result.err = math.Abs(mod.val*sin_result.err) + math.Abs(sin_result.val*mod.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_mp, stat_sin)
	} else if x < 1.0 {
		z := x * x * x
		result_c0, result_c1 := new(Result), new(Result)
		aibif.EvaluateMode(z, mode, result_c0)
		aibig.EvaluateMode(z, mode, result_c1)
		result.val = 0.625 + result_c0.val + x*(0.4375+result_c1.val)
		result.err = result_c0.err + math.Abs(x*result_c1.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x <= 2.0 {
		z := (2.0*x*x*x - 9.0) / 7.0
		result_c0, result_c1 := new(Result), new(Result)
		aibif2.EvaluateMode(z, mode, result_c0)
		aibig2.EvaluateMode(z, mode, result_c1)
		result.val = 1.125 + result_c0.val + x*(0.625+result_c1.val)
		result.err = result_c0.err + math.Abs(x*result_c1.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		y := 2.0 * x * math.Sqrt(x) / 3.0
		s := math.Exp(y)

		if y > gsl.LnMaxFloat64-1.0 {
			return OverflowError(result)
		} else {
			result_bie := new(Result)
			stat_bie := airy_bie(x, mode, result_bie)
			result.val = result_bie.val * s
			result.err = result_bie.err*s + math.Abs(1.5*y*(gsl.Float64Eps*result.val))
			result.err += gsl.Float64Eps * math.Abs(result.val)
			return stat_bie
		}
	}
}

// This computes a scaled version of the Airy function S_B(x) Bi(x).
// For x > 0 the scaling factor S_B(x) is exp(-(2/3) x^{3/2}), and is 1 for x < 0.
func Airy_Bi_scaled_e(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	if x < -1.0 {
		mod, theta, sin_result := new(Result), new(Result), new(Result)
		stat_mp := airy_mod_phase(x, mode, mod, theta)
		stat_sin := Sin_err_e(theta.val, theta.err, sin_result)
		result.val = mod.val * sin_result.val
		result.err = math.Abs(mod.val*sin_result.err) + math.Abs(sin_result.val*mod.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_mp, stat_sin)
	} else if x < 1.0 {
		z := x * x * x
		result_c0, result_c1 := new(Result), new(Result)
		aibif.EvaluateMode(z, mode, result_c0)
		aibig.EvaluateMode(z, mode, result_c1)
		result.val = 0.625 + result_c0.val + x*(0.4375+result_c1.val)
		result.err = result_c0.err + math.Abs(x*result_c1.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		if x > 0.0 {
			scale := math.Exp(-2.0 / 3.0 * math.Sqrt(z))
			result.val *= scale
			result.err *= scale
		}
		return nil
	} else if x <= 2.0 {
		x3 := x * x * x
		z := (2.0*x3 - 9.0) / 7.0
		s := math.Exp(-2.0 / 3.0 * math.Sqrt(x3))
		result_c0, result_c1 := new(Result), new(Result)
		aibif2.EvaluateMode(z, mode, result_c0)
		aibig2.EvaluateMode(z, mode, result_c1)
		result.val = s * (1.125 + result_c0.val + x*(0.625+result_c1.val))
		result.err = s * (result_c0.err + math.Abs(x*result_c1.err))
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		return airy_bie(x, mode, result)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Airy_Ai(x float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Airy_Ai_e(x, mode, result)
	return EvalResult(result, status)
}

func Airy_Ai_scaled(x float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Airy_Ai_scaled_e(x, mode, result)
	return EvalResult(result, status)
}

func Airy_Bi(x float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Airy_Bi_e(x, mode, result)
	return EvalResult(result, status)
}

func Airy_Bi_scaled(x float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Airy_Bi_scaled_e(x, mode, result)
	return EvalResult(result, status)
}
