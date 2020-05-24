/* specfunc/erfc.c
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

const (
	LogRootPi_ = 0.57236494292470008706
)

func erfc8_sum(x float64) float64 {
	/* estimates erfc(x) valid for 8 < x < 100 */
	/* This is based on index 5725 in Hart et al */

	P := []float64{
		2.97886562639399288862,
		7.409740605964741794425,
		6.1602098531096305440906,
		5.019049726784267463450058,
		1.275366644729965952479585264,
		0.5641895835477550741253201704,
	}
	Q := []float64{
		3.3690752069827527677,
		9.608965327192787870698,
		17.08144074746600431571095,
		12.0489519278551290360340491,
		9.396034016235054150430579648,
		2.260528520767326969591866945,
		1.0,
	}

	var num, den float64

	num = P[5]
	for i := 4; i >= 0; i-- {
		num = x*num + P[i]
	}

	den = Q[6]

	for i := 5; i >= 0; i-- {
		den = x*den + Q[i]
	}

	return num / den
}

func erfc8(x float64) float64 {
	e := erfc8_sum(x)
	e *= math.Exp(-x * x)
	return e
}

func log_erfc8(x float64) float64 {
	e := erfc8_sum(x)
	e = math.Log(e) - x*x
	return e
}

/* Abramowitz+Stegun, 7.1.5 */
func erfseries(x float64, result *Result) err.GSLError {
	coef := x
	e := coef
	var del float64
	for k := 1; k < 30; k++ {
		coef *= -x * x / float64(k)
		del = coef / (2.0*float64(k) + 1.0)
		e += del
	}
	result.val = 2.0 / math.SqrtPi * e
	result.err = 2.0 / math.SqrtPi * (math.Abs(del) + gsl.Float64Eps)
	return nil
}

/* Chebyshev fit for erfc((t+1)/2), -1 < t < 1
 */

var (
	erfc_xlt1_data = []float64{
		1.06073416421769980345174155056,
		-0.42582445804381043569204735291,
		0.04955262679620434040357683080,
		0.00449293488768382749558001242,
		-0.00129194104658496953494224761,
		-0.00001836389292149396270416979,
		0.00002211114704099526291538556,
		-5.23337485234257134673693179020e-7,
		-2.78184788833537885382530989578e-7,
		1.41158092748813114560316684249e-8,
		2.72571296330561699984539141865e-9,
		-2.06343904872070629406401492476e-10,
		-2.14273991996785367924201401812e-11,
		2.22990255539358204580285098119e-12,
		1.36250074650698280575807934155e-13,
		-1.95144010922293091898995913038e-14,
		-6.85627169231704599442806370690e-16,
		1.44506492869699938239521607493e-16,
		2.45935306460536488037576200030e-18,
		-9.29599561220523396007359328540e-19,
	}

	erfc_xlt1_cs = &chebyshevSeries{
		erfc_xlt1_data,
		19,
		-1, 1,
		12,
	}

	/* Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1
	 */
	erfc_x15_data = []float64{
		0.44045832024338111077637466616,
		-0.143958836762168335790826895326,
		0.044786499817939267247056666937,
		-0.013343124200271211203618353102,
		0.003824682739750469767692372556,
		-0.001058699227195126547306482530,
		0.000283859419210073742736310108,
		-0.000073906170662206760483959432,
		0.000018725312521489179015872934,
		-4.62530981164919445131297264430e-6,
		1.11558657244432857487884006422e-6,
		-2.63098662650834130067808832725e-7,
		6.07462122724551777372119408710e-8,
		-1.37460865539865444777251011793e-8,
		3.05157051905475145520096717210e-9,
		-6.65174789720310713757307724790e-10,
		1.42483346273207784489792999706e-10,
		-3.00141127395323902092018744545e-11,
		6.22171792645348091472914001250e-12,
		-1.26994639225668496876152836555e-12,
		2.55385883033257575402681845385e-13,
		-5.06258237507038698392265499770e-14,
		9.89705409478327321641264227110e-15,
		-1.90685978789192181051961024995e-15,
		3.50826648032737849245113757340e-16,
	}

	erfc_x15_cs = &chebyshevSeries{
		erfc_x15_data,
		24,
		-1, 1,
		16,
	}

	/* Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1
	 */
	erfc_x510_data = []float64{
		1.11684990123545698684297865808,
		0.003736240359381998520654927536,
		-0.000916623948045470238763619870,
		0.000199094325044940833965078819,
		-0.000040276384918650072591781859,
		7.76515264697061049477127605790e-6,
		-1.44464794206689070402099225301e-6,
		2.61311930343463958393485241947e-7,
		-4.61833026634844152345304095560e-8,
		8.00253111512943601598732144340e-9,
		-1.36291114862793031395712122089e-9,
		2.28570483090160869607683087722e-10,
		-3.78022521563251805044056974560e-11,
		6.17253683874528285729910462130e-12,
		-9.96019290955316888445830597430e-13,
		1.58953143706980770269506726000e-13,
		-2.51045971047162509999527428316e-14,
		3.92607828989125810013581287560e-15,
		-6.07970619384160374392535453420e-16,
		9.12600607264794717315507477670e-17,
	}

	erfc_x510_cs = &chebyshevSeries{
		erfc_x510_data,
		19,
		-1, 1,
		12,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

func Erfc_e(x float64, result *Result) err.GSLError {
	ax := math.Abs(x)
	var e_val, e_err float64

	if ax <= 1.0 {
		t := 2.0*ax - 1.0
		c := new(Result)
		erfc_xlt1_cs.Evaluate(t, c)
		e_val = c.val
		e_err = c.err
	} else if ax <= 5.0 {
		ex2 := math.Exp(-x * x)
		t := 0.5 * (ax - 3.0)
		c := new(Result)
		erfc_x15_cs.Evaluate(t, c)
		e_val = ex2 * c.val
		e_err = ex2 * (c.err + 2.0*math.Abs(x)*gsl.Float64Eps)
	} else if ax < 10.0 {
		exterm := math.Exp(-x*x) / ax
		t := (2.0*ax - 15.0) / 5.0
		c := new(Result)
		erfc_x510_cs.Evaluate(t, c)
		e_val = exterm * c.val
		e_err = exterm * (c.err + 2.0*math.Abs(x)*gsl.Float64Eps + gsl.Float64Eps)
	} else {
		e_val = erfc8(ax)
		e_err = (x*x + 1.0) * gsl.Float64Eps * math.Abs(e_val)
	}

	if x < 0.0 {
		result.val = 2.0 - e_val
		result.err = e_err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	} else {
		result.val = e_val
		result.err = e_err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	}

	return nil
}

func Log_erfc_e(x float64, result *Result) err.GSLError {
	if x*x < 10.0*gsl.Root6Float64Eps {
		var (
			y = x / math.SqrtPi
			/* series for -1/2 Log[Erfc[Sqrt[Pi] y]] */
			c3     = (4.0 - gsl.Pi) / 3.0
			c4     = 2.0 * (1.0 - gsl.Pi/3.0)
			c5     = -0.001829764677455021 /* (96.0 - 40.0*M_PI + 3.0*M_PI*M_PI)/30.0  */
			c6     = 0.02629651521057465   /* 2.0*(120.0 - 60.0*M_PI + 7.0*M_PI*M_PI)/45.0 */
			c7     = -0.01621575378835404
			c8     = 0.00125993961762116
			c9     = 0.00556964649138
			c10    = -0.0045563339802
			c11    = 0.0009461589032
			c12    = 0.0013200243174
			c13    = -0.00142906
			c14    = 0.00048204
			series = c8 + y*(c9+y*(c10+y*(c11+y*(c12+y*(c13+c14*y)))))
		)
		series = y * (1.0 + y*(1.0+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*series)))))))
		result.val = -2.0 * series
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if x > 8.0 {
		/*
		  don't like use of log1p(); added above series stuff for small x instead, should be ok [GJ]
		  else if (fabs(x) < 1.0) {
		    gsl_sf_result result_erf;
		    gsl_sf_erf_e(x, &result_erf);
		    result->val  = log1p(-result_erf.val);
		    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
		    return GSL_SUCCESS;
		  }
		*/

		result.val = log_erfc8(x)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		result_erfc := new(Result)
		Erfc_e(x, result_erfc)
		result.val = math.Log(result_erfc.val)
		result.err = math.Abs(result_erfc.err / result_erfc.val)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Erf_e(x float64, result *Result) err.GSLError {
	if math.Abs(x) < 1.0 {
		return erfseries(x, result)
	} else {
		result_erfc := new(Result)
		Erfc_e(x, result_erfc)
		result.val = 1.0 - result_erfc.val
		result.err = result_erfc.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

func Erf_Z_e(x float64, result *Result) err.GSLError {
	ex2 := math.Exp(-x * x / 2.0)
	result.val = ex2 / (math.Sqrt2 * math.SqrtPi)
	result.err = math.Abs(x*result.val) * gsl.Float64Eps
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	CheckUnderflow(result)
	return nil
}

func Erf_Q_e(x float64, result *Result) err.GSLError {
	result_erfc := new(Result)
	stat := Erfc_e(x/math.Sqrt2, result_erfc)
	result.val = 0.5 * result_erfc.val
	result.err = 0.5 * result_erfc.err
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return stat
}

func Hazard_e(x float64, result *Result) err.GSLError {
	if x < 25.0 {
		result_ln_erfc := new(Result)
		stat_l := Log_erfc_e(x/math.Sqrt2, result_ln_erfc)
		lnc := -0.22579135264472743236 /* ln(sqrt(2/pi)) */
		arg := lnc - 0.5*x*x - result_ln_erfc.val
		stat_e := Exp_e(arg, result)
		result.err += 3.0 * (1.0 + math.Abs(x)) * gsl.Float64Eps * math.Abs(result.val)
		result.err += math.Abs(result_ln_erfc.err * result.val)
		return err.ErrorSelect(stat_l, stat_e)
	} else {
		ix2 := 1.0 / (x * x)
		corrB := 1.0 - 9.0*ix2*(1.0-11.0*ix2)
		corrM := 1.0 - 5.0*ix2*(1.0-7.0*ix2*corrB)
		corrT := 1.0 - ix2*(1.0-3.0*ix2*corrM)
		result.val = x / corrT
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Erfc(x float64) float64 {
	result := new(Result)
	status := Erfc_e(x, result)
	return EvalResult(result, status)
}

func Log_erfc(x float64) float64 {
	result := new(Result)
	status := Log_erfc_e(x, result)
	return EvalResult(result, status)
}

func Erf(x float64) float64 {
	result := new(Result)
	status := Erf_e(x, result)
	return EvalResult(result, status)
}

func Erf_Z(x float64) float64 {
	result := new(Result)
	status := Erf_Z_e(x, result)
	return EvalResult(result, status)
}

func Erf_Q(x float64) float64 {
	result := new(Result)
	status := Erf_Q_e(x, result)
	return EvalResult(result, status)
}

func Hazard(x float64) float64 {
	result := new(Result)
	status := Hazard_e(x, result)
	return EvalResult(result, status)
}
