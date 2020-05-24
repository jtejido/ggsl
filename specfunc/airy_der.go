/* specfunc/airy_der.c
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
	aif_data = []float64{
		0.10527461226531408809,
		0.01183613628152997844,
		0.00012328104173225664,
		0.00000062261225638140,
		0.00000000185298887844,
		0.00000000000363328873,
		0.00000000000000504622,
		0.00000000000000000522,
	}

	aif = &chebyshevSeries{
		aif_data,
		7,
		-1, 1,
		7,
	}

	aig_data = []float64{
		0.021233878150918666852,
		0.086315930335214406752,
		0.001797594720383231358,
		0.000014265499875550693,
		0.000000059437995283683,
		0.000000000152403366479,
		0.000000000000264587660,
		0.000000000000000331562,
		0.000000000000000000314,
	}

	aig = &chebyshevSeries{
		aig_data,
		8,
		-1, 1,
		8,
	}

	aip2_data = []float64{
		0.0065457691989713757,
		0.0023833724120774592,
		-0.0000430700770220586,
		0.0000015629125858629,
		-0.0000000815417186163,
		0.0000000054103738057,
		-0.0000000004284130883,
		0.0000000000389497963,
		-0.0000000000039623161,
		0.0000000000004428184,
		-0.0000000000000536297,
		0.0000000000000069650,
		-0.0000000000000009620,
		0.0000000000000001403,
		-0.0000000000000000215,
	}

	aip2 = &chebyshevSeries{
		aip2_data,
		14,
		-1, 1,
		9,
	}

	aip1_data = []float64{
		0.0358865097808301538,
		0.0114668575627764899,
		-0.0007592073583861400,
		0.0000869517610893841,
		-0.0000128237294298592,
		0.0000022062695681038,
		-0.0000004222295185921,
		0.0000000874686415726,
		-0.0000000192773588418,
		0.0000000044668460054,
		-0.0000000010790108052,
		0.0000000002700029447,
		-0.0000000000696480108,
		0.0000000000184489907,
		-0.0000000000050027817,
		0.0000000000013852243,
		-0.0000000000003908218,
		0.0000000000001121536,
		-0.0000000000000326862,
		0.0000000000000096619,
		-0.0000000000000028935,
		0.0000000000000008770,
		-0.0000000000000002688,
		0.0000000000000000832,
		-0.0000000000000000260,
	}

	aip1 = &chebyshevSeries{
		aip1_data,
		24,
		-1, 1,
		14,
	}

	bif_data = []float64{
		0.1153536790828570243,
		0.0205007894049192875,
		0.0002135290278902876,
		0.0000010783960614677,
		0.0000000032094708833,
		0.0000000000062930407,
		0.0000000000000087403,
		0.0000000000000000090,
	}

	bif = &chebyshevSeries{
		bif_data,
		7,
		-1, 1,
		7,
	}

	big_data = []float64{
		-0.097196440416443537390,
		0.149503576843167066571,
		0.003113525387121326042,
		0.000024708570579821297,
		0.000000102949627731379,
		0.000000000263970373987,
		0.000000000000458279271,
		0.000000000000000574283,
		0.000000000000000000544,
	}

	big = &chebyshevSeries{
		big_data,
		8,
		-1, 1,
		8,
	}

	bif2_data = []float64{
		0.323493987603522033521,
		0.086297871535563559139,
		0.002994025552655397426,
		0.000051430528364661637,
		0.000000525840250036811,
		0.000000003561751373958,
		0.000000000017146864007,
		0.000000000000061663520,
		0.000000000000000171911,
		0.000000000000000000382,
	}

	bif2 = &chebyshevSeries{
		bif2_data,
		9,
		-1, 1,
		9,
	}

	big2_data = []float64{
		1.6062999463621294578,
		0.7449088819876088652,
		0.0470138738610277380,
		0.0012284422062548239,
		0.0000173222412256624,
		0.0000001521901652368,
		0.0000000009113560249,
		0.0000000000039547918,
		0.0000000000000130017,
		0.0000000000000000335,
	}

	big2 = &chebyshevSeries{
		big2_data,
		9,
		-1, 1,
		9,
	}

	bip2_data = []float64{
		-0.13269705443526630495,
		-0.00568443626045977481,
		-0.00015643601119611610,
		-0.00001136737203679562,
		-0.00000143464350991284,
		-0.00000018098531185164,
		0.00000000926177343611,
		0.00000001710005490721,
		0.00000000476698163504,
		-0.00000000035195022023,
		-0.00000000058890614316,
		-0.00000000006678499608,
		0.00000000006395565102,
		0.00000000001554529427,
		-0.00000000000792397000,
		-0.00000000000258326243,
		0.00000000000121655048,
		0.00000000000038707207,
		-0.00000000000022487045,
		-0.00000000000004953477,
		0.00000000000004563782,
		0.00000000000000332998,
		-0.00000000000000921750,
		0.00000000000000094157,
		0.00000000000000167154,
		-0.00000000000000055134,
		-0.00000000000000022369,
		0.00000000000000017487,
		0.00000000000000000207,
	}

	bip2 = &chebyshevSeries{
		bip2_data,
		28,
		-1, 1,
		14,
	}

	bip1_data = []float64{
		-0.1729187351079553719,
		-0.0149358492984694364,
		-0.0005471104951678566,
		0.0001537966292958408,
		0.0000154353476192179,
		-0.0000065434113851906,
		0.0000003728082407879,
		0.0000002072078388189,
		-0.0000000658173336470,
		0.0000000074926746354,
		0.0000000011101336884,
		-0.0000000007265140553,
		0.0000000001782723560,
		-0.0000000000217346352,
		-0.0000000000020302035,
		0.0000000000019311827,
		-0.0000000000006044953,
		0.0000000000001209450,
		-0.0000000000000125109,
		-0.0000000000000019917,
		0.0000000000000015154,
		-0.0000000000000004977,
		0.0000000000000001155,
		-0.0000000000000000186,
	}

	bip1 = &chebyshevSeries{
		bip1_data,
		23,
		-1, 1,
		13,
	}

	an22_data = []float64{
		0.0537418629629794329,
		-0.0126661435859883193,
		-0.0011924334106593007,
		-0.0002032327627275655,
		-0.0000446468963075164,
		-0.0000113359036053123,
		-0.0000031641352378546,
		-0.0000009446708886149,
		-0.0000002966562236472,
		-0.0000000969118892024,
		-0.0000000326822538653,
		-0.0000000113144618964,
		-0.0000000040042691002,
		-0.0000000014440333684,
		-0.0000000005292853746,
		-0.0000000001967763374,
		-0.0000000000740800096,
		-0.0000000000282016314,
		-0.0000000000108440066,
		-0.0000000000042074801,
		-0.0000000000016459150,
		-0.0000000000006486827,
		-0.0000000000002574095,
		-0.0000000000001027889,
		-0.0000000000000412846,
		-0.0000000000000166711,
		-0.0000000000000067657,
		-0.0000000000000027585,
		-0.0000000000000011296,
		-0.0000000000000004645,
		-0.0000000000000001917,
		-0.0000000000000000794,
		-0.0000000000000000330,
	}

	an22 = &chebyshevSeries{
		an22_data,
		32,
		-1, 1,
		18,
	}

	an21_data = []float64{
		0.0198313155263169394,
		-0.0029376249067087533,
		-0.0001136260695958196,
		-0.0000100554451087156,
		-0.0000013048787116563,
		-0.0000002123881993151,
		-0.0000000402270833384,
		-0.0000000084996745953,
		-0.0000000019514839426,
		-0.0000000004783865344,
		-0.0000000001236733992,
		-0.0000000000334137486,
		-0.0000000000093702824,
		-0.0000000000027130128,
		-0.0000000000008075954,
		-0.0000000000002463214,
		-0.0000000000000767656,
		-0.0000000000000243883,
		-0.0000000000000078831,
		-0.0000000000000025882,
		-0.0000000000000008619,
		-0.0000000000000002908,
		-0.0000000000000000993,
		-0.0000000000000000343,
	}

	an21 = &chebyshevSeries{
		an21_data,
		23,
		-1, 1,
		12,
	}

	an20_data = []float64{
		0.0126732217145738027,
		-0.0005212847072615621,
		-0.0000052672111140370,
		-0.0000001628202185026,
		-0.0000000090991442687,
		-0.0000000007438647126,
		-0.0000000000795494752,
		-0.0000000000104050944,
		-0.0000000000015932426,
		-0.0000000000002770648,
		-0.0000000000000535343,
		-0.0000000000000113062,
		-0.0000000000000025772,
		-0.0000000000000006278,
		-0.0000000000000001621,
		-0.0000000000000000441,
	}

	an20 = &chebyshevSeries{
		an20_data,
		15,
		-1, 1,
		8,
	}

	aph2_data = []float64{
		-0.2057088719781465107,
		0.0422196961357771922,
		0.0020482560511207275,
		0.0002607800735165006,
		0.0000474824268004729,
		0.0000105102756431612,
		0.0000026353534014668,
		0.0000007208824863499,
		0.0000002103236664473,
		0.0000000644975634555,
		0.0000000205802377264,
		0.0000000067836273921,
		0.0000000022974015284,
		0.0000000007961306765,
		0.0000000002813860610,
		0.0000000001011749057,
		0.0000000000369306738,
		0.0000000000136615066,
		0.0000000000051142751,
		0.0000000000019351689,
		0.0000000000007393607,
		0.0000000000002849792,
		0.0000000000001107281,
		0.0000000000000433412,
		0.0000000000000170801,
		0.0000000000000067733,
		0.0000000000000027017,
		0.0000000000000010835,
		0.0000000000000004367,
		0.0000000000000001769,
		0.0000000000000000719,
		0.0000000000000000294,
	}

	aph2 = &chebyshevSeries{
		aph2_data,
		31,
		-1, 1,
		16,
	}

	aph1_data = []float64{
		-0.1024172908077571694,
		0.0071697275146591248,
		0.0001209959363122329,
		0.0000073361512841220,
		0.0000007535382954272,
		0.0000001041478171741,
		0.0000000174358728519,
		0.0000000033399795033,
		0.0000000007073075174,
		0.0000000001619187515,
		0.0000000000394539982,
		0.0000000000101192282,
		0.0000000000027092778,
		0.0000000000007523806,
		0.0000000000002156369,
		0.0000000000000635283,
		0.0000000000000191757,
		0.0000000000000059143,
		0.0000000000000018597,
		0.0000000000000005950,
		0.0000000000000001934,
		0.0000000000000000638,
	}

	aph1 = &chebyshevSeries{
		aph1_data,
		21,
		-1, 1,
		10,
	}

	aph0_data = []float64{
		-0.0855849241130933257,
		0.0011214378867065261,
		0.0000042721029353664,
		0.0000000817607381483,
		0.0000000033907645000,
		0.0000000002253264423,
		0.0000000000206284209,
		0.0000000000023858763,
		0.0000000000003301618,
		0.0000000000000527010,
		0.0000000000000094555,
		0.0000000000000018709,
		0.0000000000000004024,
		0.0000000000000000930,
		0.0000000000000000229,
	}

	aph0 = &chebyshevSeries{
		aph0_data,
		14,
		-1, 1,
		7,
	}
)

func airy_deriv_mod_phase(x float64, mode gsl.MODE_T, ampl, phi *Result) err.GSLError {
	pi34 := 2.356194490192344928847
	result_a, result_p := new(Result), new(Result)

	if x <= -4.0 {
		z := 128.0/(x*x*x) + 1.0
		an20.EvaluateMode(z, mode, result_a)
		aph0.EvaluateMode(z, mode, result_p)
	} else if x <= -2.0 {
		z := (128.0/(x*x*x) + 9.0) / 7.0
		an21.EvaluateMode(z, mode, result_a)
		aph1.EvaluateMode(z, mode, result_p)
	} else if x <= -1.0 {
		z := (16.0/(x*x*x) + 9.0) / 7.0
		an22.EvaluateMode(z, mode, result_a)
		aph2.EvaluateMode(z, mode, result_p)
	} else {
		ampl.val = 0.0
		ampl.err = 0.0
		phi.val = 0.0
		phi.err = 0.0
		return err.ERROR("x is greater than 1.0", err.EDOM)
	}

	a := 0.3125 + result_a.val
	p := -0.625 + result_p.val
	sqx := math.Sqrt(-x)

	ampl.val = math.Sqrt(a * sqx)
	ampl.err = math.Abs(ampl.val) * (gsl.Float64Eps + math.Abs(result_a.err/result_a.val))
	phi.val = pi34 - x*sqx*p
	phi.err = math.Abs(phi.val) * (gsl.Float64Eps + math.Abs(result_p.err/result_p.val))

	return nil
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

func Airy_Ai_deriv_scaled_e(x float64, mode gsl.MODE_T, result *Result) err.GSLError {

	if x < -1.0 {
		a, p := new(Result), new(Result)
		status_ap := airy_deriv_mod_phase(x, mode, a, p)
		c := math.Cos(p.val)
		result.val = a.val * c
		result.err = math.Abs(result.val*p.err) + math.Abs(c*a.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return status_ap
	} else if x <= 1.0 {
		result_c0, result_c1 := new(Result), new(Result)
		x3 := x * x * x
		x2 := x * x
		aif.EvaluateMode(x3, mode, result_c0)
		aig.EvaluateMode(x3, mode, result_c1)

		result.val = (x2*(0.125+result_c0.val) - result_c1.val) - 0.25
		result.err = math.Abs(x2*result_c0.val) + result_c1.err
		result.err += gsl.Float64Eps * math.Abs(result.val)

		if x > gsl.Root3Float64Eps*gsl.Root3Float64Eps {
			/* scale only if x is positive */
			s := math.Exp(2.0 * x * math.Sqrt(x) / 3.0)
			result.val *= s
			result.err *= s
		}

		return nil
	} else if x <= 4.0 {
		result_c0 := new(Result)
		sqrtx := math.Sqrt(x)
		z := (16.0/(x*sqrtx) - 9.0) / 7.0
		s := math.Sqrt(sqrtx)
		aip1.EvaluateMode(z, mode, result_c0)

		result.val = -(0.28125 + result_c0.val) * s
		result.err = result_c0.err * s
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
	result_c0 := new(Result)
	sqrtx := math.Sqrt(x)
	z := 16.0/(x*sqrtx) - 1.0
	s := math.Sqrt(sqrtx)
	aip2.EvaluateMode(z, mode, result_c0)

	result.val = -(0.28125 + result_c0.val) * s
	result.err = result_c0.err * s
	result.err += gsl.Float64Eps * math.Abs(result.val)
	return nil

}

func Airy_Ai_deriv_e(x float64, mode gsl.MODE_T, result *Result) err.GSLError {

	if x < -1.0 {
		a, p := new(Result), new(Result)
		status_ap := airy_deriv_mod_phase(x, mode, a, p)
		c := math.Cos(p.val)
		result.val = a.val * c
		result.err = math.Abs(result.val*p.err) + math.Abs(c*a.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return status_ap
	} else if x < 1.0 {
		x3 := x * x * x
		result_c1, result_c2 := new(Result), new(Result)
		aif.EvaluateMode(x3, mode, result_c1)
		aig.EvaluateMode(x3, mode, result_c2)
		result.val = (x*x*(0.125+result_c1.val) - result_c2.val) - 0.25
		result.err = math.Abs(x*x*result_c1.err) + result_c2.err
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x*x*x < 9.0/4.0*gsl.LnMinFloat64*gsl.LnMinFloat64 {
		result_aps := new(Result)
		arg := -2.0 * x * math.Sqrt(x) / 3.0
		stat_a := Airy_Ai_deriv_scaled_e(x, mode, result_aps)
		stat_e := Exp_mult_err_e(arg, 1.5*math.Abs(arg*gsl.Float64Eps), result_aps.val, result_aps.err, result)

		return err.ErrorSelect(stat_e, stat_a)
	}

	return UnderflowError(result)
}

func Airy_Bi_deriv_scaled_e(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	atr := 8.7506905708484345
	btr := -2.0938363213560543

	if x < -1.0 {
		a, p := new(Result), new(Result)
		status_ap := airy_deriv_mod_phase(x, mode, a, p)
		s := math.Sin(p.val)
		result.val = a.val * s
		result.err = math.Abs(result.val*p.err) + math.Abs(s*a.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return status_ap
	} else if x < 1.0 {
		x3 := x * x * x
		x2 := x * x
		result_c1, result_c2 := new(Result), new(Result)
		bif.EvaluateMode(x3, mode, result_c1)
		big.EvaluateMode(x3, mode, result_c2)
		result.val = x2*(result_c1.val+0.25) + result_c2.val + 0.5
		result.err = x2*result_c1.err + result_c2.err
		result.err += gsl.Float64Eps * math.Abs(result.val)

		if x > gsl.Root3Float64Eps*gsl.Root3Float64Eps {
			s := math.Exp(-2.0 * x * math.Sqrt(x) / 3.0)
			result.val *= s
			result.err *= s
		}

		return nil
	} else if x < 2.0 {
		z := (2.0*x*x*x - 9.0) / 7.0
		s := math.Exp(-2.0 * x * math.Sqrt(x) / 3.0)
		result_c0, result_c1 := new(Result), new(Result)
		bif2.EvaluateMode(z, mode, result_c0)
		big2.EvaluateMode(z, mode, result_c1)

		result.val = s * (x*x*(0.25+result_c0.val) + 0.5 + result_c1.val)
		result.err = s * (x*x*result_c0.err + result_c1.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < 4.0 {
		sqrtx := math.Sqrt(x)
		z := atr/(x*sqrtx) + btr
		s := math.Sqrt(sqrtx)
		result_c0 := new(Result)
		bip1.EvaluateMode(z, mode, result_c0)
		result.val = s * (0.625 + result_c0.val)
		result.err = s * result_c0.err
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	}

	sqrtx := math.Sqrt(x)
	z := 16.0/(x*sqrtx) - 1.0
	s := math.Sqrt(sqrtx)
	result_c0 := new(Result)
	bip2.EvaluateMode(z, mode, result_c0)

	result.val = s * (0.625 + result_c0.val)
	result.err = s * result_c0.err
	result.err += gsl.Float64Eps * math.Abs(result.val)
	return nil

}

func Airy_Bi_deriv_e(x float64, mode gsl.MODE_T, result *Result) err.GSLError {
	if x < -1.0 {
		a, p := new(Result), new(Result)
		status_ap := airy_deriv_mod_phase(x, mode, a, p)
		s := math.Sin(p.val)

		result.val = a.val * s
		result.err = math.Abs(result.val*p.err) + math.Abs(s*a.err)
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return status_ap
	} else if x < 1.0 {
		x3 := x * x * x
		x2 := x * x
		result_c1, result_c2 := new(Result), new(Result)
		bif.EvaluateMode(x3, mode, result_c1)
		big.EvaluateMode(x3, mode, result_c2)

		result.val = x2*(result_c1.val+0.25) + result_c2.val + 0.5
		result.err = x2*result_c1.err + result_c2.err
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < 2.0 {
		z := (2.0*x*x*x - 9.0) / 7.0
		result_c1, result_c2 := new(Result), new(Result)
		bif2.EvaluateMode(z, mode, result_c1)
		big2.EvaluateMode(z, mode, result_c2)

		result.val = x*x*(result_c1.val+0.25) + result_c2.val + 0.5
		result.err = x*x*result_c1.err + result_c2.err
		result.err += gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x < gsl.Root3MaxFloat64*gsl.Root3MaxFloat64 {
		arg := 2.0 * (x * math.Sqrt(x) / 3.0)
		result_bps := new(Result)
		stat_b := Airy_Bi_deriv_scaled_e(x, mode, result_bps)
		stat_e := Exp_mult_err_e(arg, 1.5*math.Abs(arg*gsl.Float64Eps), result_bps.val, result_bps.err, result)

		return err.ErrorSelect(stat_e, stat_b)
	}

	return OverflowError(result)
}
