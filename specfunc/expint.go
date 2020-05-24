/* specfunc/expint.c
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
	ae11_data = []float64{
		0.121503239716065790,
		-0.065088778513550150,
		0.004897651357459670,
		-0.000649237843027216,
		0.000093840434587471,
		0.000000420236380882,
		-0.000008113374735904,
		0.000002804247688663,
		0.000000056487164441,
		-0.000000344809174450,
		0.000000058209273578,
		0.000000038711426349,
		-0.000000012453235014,
		-0.000000005118504888,
		0.000000002148771527,
		0.000000000868459898,
		-0.000000000343650105,
		-0.000000000179796603,
		0.000000000047442060,
		0.000000000040423282,
		-0.000000000003543928,
		-0.000000000008853444,
		-0.000000000000960151,
		0.000000000001692921,
		0.000000000000607990,
		-0.000000000000224338,
		-0.000000000000200327,
		-0.000000000000006246,
		0.000000000000045571,
		0.000000000000016383,
		-0.000000000000005561,
		-0.000000000000006074,
		-0.000000000000000862,
		0.000000000000001223,
		0.000000000000000716,
		-0.000000000000000024,
		-0.000000000000000201,
		-0.000000000000000082,
		0.000000000000000017,
	}

	ae11_cs = &chebyshevSeries{
		ae11_data,
		38,
		-1, 1,
		20,
	}

	ae12_data = []float64{
		0.582417495134726740,
		-0.158348850905782750,
		-0.006764275590323141,
		0.005125843950185725,
		0.000435232492169391,
		-0.000143613366305483,
		-0.000041801320556301,
		-0.000002713395758640,
		0.000001151381913647,
		0.000000420650022012,
		0.000000066581901391,
		0.000000000662143777,
		-0.000000002844104870,
		-0.000000000940724197,
		-0.000000000177476602,
		-0.000000000015830222,
		0.000000000002905732,
		0.000000000001769356,
		0.000000000000492735,
		0.000000000000093709,
		0.000000000000010707,
		-0.000000000000000537,
		-0.000000000000000716,
		-0.000000000000000244,
		-0.000000000000000058,
	}

	ae12_cs = &chebyshevSeries{
		ae12_data,
		24,
		-1, 1,
		15,
	}

	e11_data = []float64{
		-16.11346165557149402600,
		7.79407277874268027690,
		-1.95540581886314195070,
		0.37337293866277945612,
		-0.05692503191092901938,
		0.00721107776966009185,
		-0.00078104901449841593,
		0.00007388093356262168,
		-0.00000620286187580820,
		0.00000046816002303176,
		-0.00000003209288853329,
		0.00000000201519974874,
		-0.00000000011673686816,
		0.00000000000627627066,
		-0.00000000000031481541,
		0.00000000000001479904,
		-0.00000000000000065457,
		0.00000000000000002733,
		-0.00000000000000000108,
	}

	e11_cs = &chebyshevSeries{
		e11_data,
		18,
		-1, 1,
		13,
	}

	e12_data = []float64{
		-0.03739021479220279500,
		0.04272398606220957700,
		-0.13031820798497005440,
		0.01441912402469889073,
		-0.00134617078051068022,
		0.00010731029253063780,
		-0.00000742999951611943,
		0.00000045377325690753,
		-0.00000002476417211390,
		0.00000000122076581374,
		-0.00000000005485141480,
		0.00000000000226362142,
		-0.00000000000008635897,
		0.00000000000000306291,
		-0.00000000000000010148,
		0.00000000000000000315,
	}

	e12_cs = &chebyshevSeries{
		e12_data,
		15,
		-1, 1,
		10,
	}

	ae13_data = []float64{
		-0.605773246640603460,
		-0.112535243483660900,
		0.013432266247902779,
		-0.001926845187381145,
		0.000309118337720603,
		-0.000053564132129618,
		0.000009827812880247,
		-0.000001885368984916,
		0.000000374943193568,
		-0.000000076823455870,
		0.000000016143270567,
		-0.000000003466802211,
		0.000000000758754209,
		-0.000000000168864333,
		0.000000000038145706,
		-0.000000000008733026,
		0.000000000002023672,
		-0.000000000000474132,
		0.000000000000112211,
		-0.000000000000026804,
		0.000000000000006457,
		-0.000000000000001568,
		0.000000000000000383,
		-0.000000000000000094,
		0.000000000000000023,
	}

	ae13_cs = &chebyshevSeries{
		ae13_data,
		24,
		-1, 1,
		15,
	}

	ae14_data = []float64{
		-0.18929180007530170,
		-0.08648117855259871,
		0.00722410154374659,
		-0.00080975594575573,
		0.00010999134432661,
		-0.00001717332998937,
		0.00000298562751447,
		-0.00000056596491457,
		0.00000011526808397,
		-0.00000002495030440,
		0.00000000569232420,
		-0.00000000135995766,
		0.00000000033846628,
		-0.00000000008737853,
		0.00000000002331588,
		-0.00000000000641148,
		0.00000000000181224,
		-0.00000000000052538,
		0.00000000000015592,
		-0.00000000000004729,
		0.00000000000001463,
		-0.00000000000000461,
		0.00000000000000148,
		-0.00000000000000048,
		0.00000000000000016,
		-0.00000000000000005,
	}

	ae14_cs = &chebyshevSeries{
		ae14_data,
		25,
		-1, 1,
		13,
	}
)

/* implementation for E1, allowing for scaling by exp(x) */
func expint_E1_impl(x float64, result *Result, scale int) err.GSLError {
	xmaxt := -gsl.LnMinFloat64      /* XMAXT = -LOG (R1MACH(1)) */
	xmax := xmaxt - math.Log(xmaxt) /* XMAX = XMAXT - LOG(XMAXT) */

	/* CHECK_POINTER(result) */

	if x < -xmax && scale != 1 {
		return OverflowError(result)
	} else if x <= -10.0 {
		s := 1.0 / x
		if scale != 1 {
			s *= math.Exp(-x)
		}

		result_c := new(Result)
		ae11_cs.Evaluate(20.0/x+1.0, result_c)
		result.val = s * (1.0 + result_c.val)
		result.err = s * result_c.err
		result.err += 2.0 * gsl.Float64Eps * (math.Abs(x) + 1.0) * math.Abs(result.val)
		return nil
	} else if x <= -4.0 {
		s := 1.0 / x
		if scale != 1 {
			s *= math.Exp(-x)
		}
		result_c := new(Result)
		ae12_cs.Evaluate((40.0/x+7.0)/3.0, result_c)
		result.val = s * (1.0 + result_c.val)
		result.err = s * result_c.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x <= -1.0 {
		ln_term := -math.Log(math.Abs(x))
		scale_factor := 1.
		if scale == 1 {
			scale_factor = math.Exp(x)
		}

		result_c := new(Result)
		e11_cs.Evaluate((2.0*x+5.0)/3.0, result_c)
		result.val = scale_factor * (ln_term + result_c.val)
		result.err = scale_factor * (result_c.err + gsl.Float64Eps*math.Abs(ln_term))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x == 0.0 {
		return DomainError(result)
	} else if x <= 1.0 {
		ln_term := -math.Log(math.Abs(x))
		scale_factor := 1.
		if scale == 1 {
			scale_factor = math.Exp(x)
		}
		result_c := new(Result)
		e12_cs.Evaluate(x, result_c)
		result.val = scale_factor * (ln_term - 0.6875 + x + result_c.val)
		result.err = scale_factor * (result_c.err + gsl.Float64Eps*math.Abs(ln_term))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x <= 4.0 {
		s := 1.0 / x
		if scale != 1 {
			s *= math.Exp(-x)
		}

		result_c := new(Result)
		ae13_cs.Evaluate((8.0/x-5.0)/3.0, result_c)
		result.val = s * (1.0 + result_c.val)
		result.err = s * result_c.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x <= xmax || scale == 1 {
		s := 1.0 / x
		if scale != 1 {
			s *= math.Exp(-x)
		}
		result_c := new(Result)
		ae14_cs.Evaluate(8.0/x-1.0, result_c)
		result.val = s * (1.0 + result_c.val)
		result.err = s * (gsl.Float64Eps + result_c.err)
		result.err += 2.0 * (x + 1.0) * gsl.Float64Eps * math.Abs(result.val)
		if result.val == 0.0 {
			return UnderflowError(result)
		}

		return nil
	}

	return UnderflowError(result)

}

func expint_E2_impl(x float64, result *Result, scale int) err.GSLError {
	xmaxt := -gsl.LnMinFloat64
	xmax := xmaxt - math.Log(xmaxt)

	if x < -xmax && scale != 1 {
		return OverflowError(result)
	} else if x == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if x < 100.0 {
		ex := 1.0
		if scale != 1 {
			ex = math.Exp(-x)
		}

		result_E1 := new(Result)
		stat_E1 := expint_E1_impl(x, result_E1, scale)
		result.val = ex - x*result_E1.val
		result.err = gsl.Float64Eps*ex + math.Abs(x)*result_E1.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_E1
	} else if x < xmax || scale == 1 {
		s := 1.
		if scale != 1 {
			s = math.Exp(-x)
		}

		var (
			c1   = -2.0
			c2   = 6.0
			c3   = -24.0
			c4   = 120.0
			c5   = -720.0
			c6   = 5040.0
			c7   = -40320.0
			c8   = 362880.0
			c9   = -3628800.0
			c10  = 39916800.0
			c11  = -479001600.0
			c12  = 6227020800.0
			c13  = -87178291200.0
			y    = 1.0 / x
			sum6 = c6 + y*(c7+y*(c8+y*(c9+y*(c10+y*(c11+y*(c12+y*c13))))))
			sum  = y * (c1 + y*(c2+y*(c3+y*(c4+y*(c5+y*sum6)))))
		)

		result.val = s * (1.0 + sum) / x
		result.err = 2.0 * (x + 1.0) * gsl.Float64Eps * result.val
		if result.val == 0.0 {
			return UnderflowError(result)
		}

		return nil
	}

	return UnderflowError(result)

}

func expint_En_impl(n int, x float64, result *Result, scale int) err.GSLError {
	if n < 0 {
		return DomainError(result)
	} else if n == 0 {
		if x == 0 {
			return DomainError(result)
		} else {
			if scale == 1 {
				result.val = 1.
			} else {
				result.val = math.Exp(-x)
			}
			result.val /= x

			result.err = 2 * gsl.Float64Eps * math.Abs(result.val)
			CheckUnderflow(result)
			return nil
		}
	} else if n == 1 {
		return expint_E1_impl(x, result, scale)
	} else if n == 2 {
		return expint_E2_impl(x, result, scale)
	} else {
		if x < 0 {
			return DomainError(result)
		}
		if x == 0 {
			if scale == 1 {
				result.val = math.Exp(x)
			} else {
				result.val = 1
			}

			result.val *= (1 / (float64(n) - 1.0))

			result.err = 2 * gsl.Float64Eps * math.Abs(result.val)
			CheckUnderflow(result)
			return nil
		} else {
			result_g := new(Result)
			prefactor := math.Pow(x, float64(n)-1)
			status := Gamma_inc_e(1-float64(n), x, result_g)
			scale_factor := 1.0
			if scale == 1 {
				scale_factor = math.Exp(x)
			}
			result.val = scale_factor * prefactor * result_g.val
			result.err = 2 * gsl.Float64Eps * math.Abs(result.val)
			result.err += 2 * math.Abs(scale_factor*prefactor*result_g.err)
			if status == nil {
				CheckUnderflow(result)
			}
			return status
		}
	}
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

func Expint_E1_e(x float64, result *Result) err.GSLError {
	return expint_E1_impl(x, result, 0)
}

func Expint_E1_scaled_e(x float64, result *Result) err.GSLError {
	return expint_E1_impl(x, result, 1)
}

func Expint_E2_e(x float64, result *Result) err.GSLError {
	return expint_E2_impl(x, result, 0)
}

func Expint_E2_scaled_e(x float64, result *Result) err.GSLError {
	return expint_E2_impl(x, result, 1)
}

func Expint_En_e(n int, x float64, result *Result) err.GSLError {
	return expint_En_impl(n, x, result, 0)
}

func Expint_En_scaled_e(n int, x float64, result *Result) err.GSLError {
	return expint_En_impl(n, x, result, 1)
}

func Expint_Ei_e(x float64, result *Result) err.GSLError {
	status := Expint_E1_e(-x, result)
	result.val = -result.val
	return status
}

func Expint_Ei_scaled_e(x float64, result *Result) err.GSLError {
	status := Expint_E1_scaled_e(-x, result)
	result.val = -result.val
	return status
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

func Expint_E1(x float64) float64 {
	result := new(Result)
	status := Expint_E1_e(x, result)
	return EvalResult(result, status)
}

func Expint_E1_scaled(x float64) float64 {
	result := new(Result)
	status := Expint_E1_scaled_e(x, result)
	return EvalResult(result, status)
}

func Expint_E2(x float64) float64 {
	result := new(Result)
	status := Expint_E2_e(x, result)
	return EvalResult(result, status)
}

func Expint_E2_scaled(x float64) float64 {
	result := new(Result)
	status := Expint_E2_scaled_e(x, result)
	return EvalResult(result, status)
}

func Expint_En(n int, x float64) float64 {
	result := new(Result)
	status := Expint_En_e(n, x, result)
	return EvalResult(result, status)
}

func Expint_En_scaled(n int, x float64) float64 {
	result := new(Result)
	status := Expint_En_scaled_e(n, x, result)
	return EvalResult(result, status)
}

func Expint_Ei(x float64) float64 {
	result := new(Result)
	status := Expint_Ei_e(x, result)
	return EvalResult(result, status)
}

func Expint_Ei_scaled(x float64) float64 {
	result := new(Result)
	status := Expint_Ei_scaled_e(x, result)
	return EvalResult(result, status)
}
