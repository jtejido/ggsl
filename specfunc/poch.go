/* specfunc/poch.c
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

const (
	GAMMA_XMAX      = 171.0
	FACT_NMAX       = 170
	DOUBLEFACT_NMAX = 297
)

var (
	bern = [21]float64{
		0.0,
		+0.833333333333333333333333333333333e-01,
		-0.138888888888888888888888888888888e-02,
		+0.330687830687830687830687830687830e-04,
		-0.826719576719576719576719576719576e-06,
		+0.208767569878680989792100903212014e-07,
		-0.528419013868749318484768220217955e-09,
		+0.133825365306846788328269809751291e-10,
		-0.338968029632258286683019539124944e-12,
		+0.858606205627784456413590545042562e-14,
		-0.217486869855806187304151642386591e-15,
		+0.550900282836022951520265260890225e-17,
		-0.139544646858125233407076862640635e-18,
		+0.353470703962946747169322997780379e-20,
		-0.895351742703754685040261131811274e-22,
		+0.226795245233768306031095073886816e-23,
		-0.574472439520264523834847971943400e-24,
		+0.145517247561486490186626486727132e-26,
		-0.368599494066531017818178247990866e-28,
		+0.933673425709504467203255515278562e-30,
		-0.236502241570062993455963519636983e-31,
	}
)

func pochrel_smallx(a, x float64, result *Result) err.GSLError {
	SQTBIG := 1.0 / (2.0 * math.Sqrt2 * gsl.Sqrt3 * gsl.SqrtMinFloat64)
	ALNEPS := gsl.LnFloat64Eps - math.Ln2

	if x == 0.0 {
		return Psi_e(a, result)
	} else {
		var bp float64
		var incr int
		if a < -0.5 {
			bp = 1.0 - a - x
		} else {
			bp = a
		}

		if bp < 10.0 {
			incr = 11 - int(bp)
		}

		b := bp + float64(incr)
		var dpoch1 float64
		dexprl := new(Result)
		v := b + 0.5*(x-1.0)
		alnvar := math.Log(v)
		q := x * alnvar

		poly1 := 0.0

		if v < SQTBIG {
			nterms := int(-0.5*ALNEPS/alnvar + 1.0)
			var2 := (1.0 / v) / v
			rho := 0.5 * (x + 1.0)
			term := var2
			var gbern [24]float64

			gbern[1] = 1.0
			gbern[2] = -rho / 12.0
			poly1 = gbern[2] * term

			if nterms > 20 {
				/* NTERMS IS TOO BIG, MAYBE D1MACH(3) IS BAD */
				/* nterms = 20; */
				result.val = 0.0
				result.err = 0.0
				return err.ERROR("error", err.ESANITY)
			}

			for kk := 2; kk <= nterms; kk++ {
				k := float64(kk)
				gbk := 0.0
				for j := 1; j <= kk; j++ {
					gbk += bern[kk-j+1] * gbern[j]
				}
				gbern[kk+1] = -rho * gbk / k

				term *= (2*k - 2 - x) * (2*k - 1 - x) * var2
				poly1 += gbern[kk+1] * term
			}
		}

		stat_dexprl := Expm1_e(q, dexprl)
		if stat_dexprl != nil {
			result.val = 0.0
			result.err = 0.0
			return stat_dexprl
		}
		dexprl.val = dexprl.val / q
		poly1 *= (x - 1.0)
		dpoch1 = dexprl.val*(alnvar+q*poly1) + poly1

		for i := incr - 1; i >= 0; i-- {
			/*
			   WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
			   TO OBTAIN DPOCH1(BP,X).
			*/
			binv := 1.0 / (bp + float64(i))
			dpoch1 = (dpoch1 - binv) / (1.0 + x*binv)
		}

		if bp == a {
			result.val = dpoch1
			result.err = 2.0 * gsl.Float64Eps * (math.Abs(float64(incr)) + 1.0) * math.Abs(result.val)
			return nil
		}

		sinpxx := math.Sin(gsl.Pi*x) / x
		sinpx2 := math.Sin(0.5 * gsl.Pi * x)
		t1 := sinpxx / math.Tan(gsl.Pi*b)
		t2 := 2.0 * sinpx2 * (sinpx2 / x)
		trig := t1 - t2

		result.val = dpoch1*(1.0+x*trig) + trig
		result.err = (math.Abs(dpoch1*x) + 1.0) * gsl.Float64Eps * (math.Abs(t1) + math.Abs(t2))
		result.err += 2.0 * gsl.Float64Eps * (math.Abs(float64(incr)) + 1.0) * math.Abs(result.val)
		return nil
	}

}

func lnpoch_pos(a, x float64, result *Result) err.GSLError {
	absx := math.Abs(x)

	if absx > 0.1*a || absx*math.Log(math.Max(a, 2.0)) > 0.1 {
		if a < GAMMA_XMAX && a+x < GAMMA_XMAX {
			g1, g2 := new(Result), new(Result)
			Gammainv_e(a, g1)
			Gammainv_e(a+x, g2)
			result.val = -math.Log(g2.val / g1.val)
			result.err = g1.err/math.Abs(g1.val) + g2.err/math.Abs(g2.val)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return nil
		} else {
			lg1, lg2 := new(Result), new(Result)
			stat_1 := Lngamma_e(a, lg1)
			stat_2 := Lngamma_e(a+x, lg2)
			result.val = lg2.val - lg1.val
			result.err = lg2.err + lg1.err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return err.ErrorSelect(stat_1, stat_2)
		}
	} else if absx < 0.1*a && a > 15.0 {
		eps := x / a
		den := 1.0 + eps
		d3 := den * den * den
		d5 := d3 * den * den
		d7 := d5 * den * den
		c1 := -eps / den
		c3 := -eps * (3.0 + eps*(3.0+eps)) / d3
		c5 := -eps * (5.0 + eps*(10.0+eps*(10.0+eps*(5.0+eps)))) / d5
		c7 := -eps * (7.0 + eps*(21.0+eps*(35.0+eps*(35.0+eps*(21.0+eps*(7.0+eps)))))) / d7
		p8 := math.Pow(1.0+eps, 8)
		c8 := 1.0/p8 - 1.0
		c9 := 1.0/(p8*(1.0+eps)) - 1.0
		a4 := a * a * a * a
		a6 := a4 * a * a
		ser_1 := c1 + c3/(30.0*a*a) + c5/(105.0*a4) + c7/(140.0*a6)
		ser_2 := c8/(99.0*a6*a*a) - 691.0/360360.0*c9/(a6*a4)
		ser := (ser_1 + ser_2) / (12.0 * a)

		term1 := x * math.Log(a/math.E)
		var term2 float64

		ln_1peps := new(Result)
		Log_1plusx_e(eps, ln_1peps)
		term2 = (x + a - 0.5) * ln_1peps.val

		result.val = term1 + term2 + ser
		result.err = gsl.Float64Eps * math.Abs(term1)
		result.err += math.Abs((x + a - 0.5) * ln_1peps.err)
		result.err += math.Abs(ln_1peps.val) * gsl.Float64Eps * (math.Abs(x) + math.Abs(a) + 0.5)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		poch_rel := new(Result)
		stat_p := pochrel_smallx(a, x, poch_rel)
		eps := x * poch_rel.val
		stat_e := Log_1plusx_e(eps, result)
		result.err = 2.0 * math.Abs(x*poch_rel.err/(1.0+eps))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return err.ErrorSelect(stat_e, stat_p)
	}
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

func Lnpoch_e(a, x float64, result *Result) err.GSLError {
	if a <= 0.0 || a+x <= 0.0 {
		return DomainError(result)
	} else if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		return lnpoch_pos(a, x, result)
	}
}

func Lnpoch_sgn_e(a, x float64, result *Result, sgn *float64) err.GSLError {

	if x == 0.0 {
		*sgn = 1.
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if a > 0.0 && a+x > 0.0 {
		*sgn = 1.
		return lnpoch_pos(a, x, result)
	} else if a <= 0 && a == math.Floor(a) {

		if a+x < 0 && x == math.Floor(x) {
			result_pos := new(Result)
			stat := lnpoch_pos(-a, -x, result_pos)
			f := math.Log(a / (a + x))
			s := -1.
			if math.Mod(x, 2) == 0 {
				s = 1.
			}
			result.val = f - result_pos.val
			result.err = result_pos.err + 2.0*gsl.Float64Eps*f
			*sgn = s
			return stat
		} else if a+x == 0 {
			stat := Lngamma_sgn_e(-a+1, result, sgn)
			s := -1.
			if math.Mod(-a, 2) == 0 {
				s = 1.
			}
			*sgn *= s
			return stat
		} else {
			result.val = math.Inf(-1)
			result.err = 0.0
			*sgn = 1
			return nil
		}
	} else if a < 0.0 && a+x < 0.0 {
		sin_1 := math.Sin(gsl.Pi * (1.0 - a))
		sin_2 := math.Sin(gsl.Pi * (1.0 - a - x))
		if sin_1 == 0.0 || sin_2 == 0.0 {
			*sgn = 0.0
			return DomainError(result)
		} else {
			lnp_pos := new(Result)
			stat_pp := lnpoch_pos(1.0-a, -x, lnp_pos)
			lnterm := math.Log(math.Abs(sin_1 / sin_2))
			result.val = lnterm - lnp_pos.val
			result.err = lnp_pos.err
			result.err += 2.0 * gsl.Float64Eps * (math.Abs(1.0-a) + math.Abs(1.0-a-x)) * math.Abs(lnterm)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			*sgn = gsl.Sign(sin_1 * sin_2)
			return stat_pp
		}
	}
	lg_apn, lg_a := new(Result), new(Result)
	var s_apn, s_a float64
	stat_apn := Lngamma_sgn_e(a+x, lg_apn, &s_apn)
	stat_a := Lngamma_sgn_e(a, lg_a, &s_a)

	if stat_apn == nil && stat_a == nil {
		result.val = lg_apn.val - lg_a.val
		result.err = lg_apn.err + lg_a.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		*sgn = s_a * s_apn
		return nil
	} else {
		if stat_apn.Status() == err.EDOM || stat_a.Status() == err.EDOM {
			*sgn = 0.0
			return DomainError(result)
		} else {
			result.val = 0.0
			result.err = 0.0
			*sgn = 0.0
			return err.Failure()
		}
	}

}

func Poch_e(a, x float64, result *Result) err.GSLError {
	if x == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else {
		lnpoch := new(Result)
		var sgn float64
		stat_lnpoch := Lnpoch_sgn_e(a, x, lnpoch, &sgn)
		if math.IsInf(lnpoch.val, -1) {
			result.val = 0
			result.err = 0
			return stat_lnpoch
		} else {
			stat_exp := Exp_err_e(lnpoch.val, lnpoch.err, result)
			result.val *= sgn
			result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
			return err.ErrorSelect(stat_exp, stat_lnpoch)
		}
	}
}

func Pochrel_e(a, x float64, result *Result) err.GSLError {
	absx := math.Abs(x)
	absa := math.Abs(a)

	if absx > 0.1*absa || absx*math.Log(math.Max(absa, 2.0)) > 0.1 {
		var sgn float64
		lnpoch := new(Result)
		stat_poch := Lnpoch_sgn_e(a, x, lnpoch, &sgn)
		if lnpoch.val > gsl.LnMaxFloat64 {
			return OverflowError(result)
		} else {
			el := math.Exp(lnpoch.val)
			result.val = (sgn*el - 1.0) / x
			result.err = math.Abs(result.val) * (lnpoch.err + 2.0*gsl.Float64Eps)
			result.err += 2.0 * gsl.Float64Eps * (math.Abs(sgn*el) + 1.0) / math.Abs(x)
			return stat_poch
		}

	}

	return pochrel_smallx(a, x, result)

}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

func Lnpoch(a, x float64) float64 {
	result := new(Result)
	status := Lnpoch_e(a, x, result)
	return EvalResult(result, status)
}

func Poch(a, x float64) float64 {
	result := new(Result)
	status := Poch_e(a, x, result)
	return EvalResult(result, status)
}

func Pochrel(a, x float64) float64 {
	result := new(Result)
	status := Pochrel_e(a, x, result)
	return EvalResult(result, status)
}
