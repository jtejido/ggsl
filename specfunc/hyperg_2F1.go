/* specfunc/hyperg_2F1.c
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

/* Assumes c != negative integer.
 */
func hyperg_2F1_series(a, b, c, x float64, result *Result) err.GSLError {
	var (
		sum_pos  = 1.0
		sum_neg  = 0.0
		del_pos  = 1.0
		del_neg  = 0.0
		del      = 1.0
		del_prev float64
		k        = 0.0
		i        = 0
	)

	if math.Abs(c) < gsl.Float64Eps {
		result.val = 0.0 /* FIXME: ?? */
		result.err = 1.0
		return err.ERROR("error", err.EDOM)
	}

	for math.Abs((del_pos+del_neg)/(sum_pos-sum_neg)) > gsl.Float64Eps {
		i++
		if i > 30000 {
			result.val = sum_pos - sum_neg
			result.err = del_pos + del_neg
			result.err += 2.0 * gsl.Float64Eps * (sum_pos + sum_neg)
			result.err += 2.0 * gsl.Float64Eps * (2.0*math.Sqrt(k) + 1.0) * math.Abs(result.val)
			return err.ERROR("error", err.EMAXITER)
		}
		del_prev = del
		del *= (a + k) * (b + k) * x / ((c + k) * (k + 1.0)) /* Gauss series */

		if del > 0.0 {
			del_pos = del
			sum_pos += del
		} else if del == 0.0 {
			/* Exact termination (a or b was a negative integer).
			 */
			del_pos = 0.0
			del_neg = 0.0
			break
		} else {
			del_neg = -del
			sum_neg -= del
		}

		/*
		 * This stopping criteria is taken from the thesis
		 * "Computation of Hypergeometic Functions" by J. Pearson, pg. 31
		 * (http://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf)
		 * and fixes bug #45926
		 */
		if math.Abs(del_prev/(sum_pos-sum_neg)) < gsl.Float64Eps && math.Abs(del/(sum_pos-sum_neg)) < gsl.Float64Eps {
			break
		}

		k += 1.0
	}

	result.val = sum_pos - sum_neg
	result.err = del_pos + del_neg
	result.err += 2.0 * gsl.Float64Eps * (sum_pos + sum_neg)
	result.err += 2.0 * gsl.Float64Eps * (2.0*math.Sqrt(k) + 1.0) * math.Abs(result.val)

	return nil
}

/* a = aR + i aI, b = aR - i aI */
func hyperg_2F1_conj_series(aR, aI, c, x float64, result *Result) err.GSLError {
	if c == 0.0 {
		result.val = 0.0 /* FIXME: should be Inf */
		result.err = 0.0
		return err.ERROR("error", err.EDOM)
	} else {
		var (
			sum_pos = 1.0
			sum_neg = 0.0
			del_pos = 1.0
			del_neg = 0.0
			del     = 1.0
			k       = 0.0
		)
		for math.Abs((del_pos+del_neg)/(sum_pos-sum_neg)) > gsl.Float64Eps {
			del *= ((aR+k)*(aR+k) + aI*aI) / ((k + 1.0) * (c + k)) * x

			if del >= 0.0 {
				del_pos = del
				sum_pos += del
			} else {
				del_neg = -del
				sum_neg -= del
			}

			if k > 30000 {
				result.val = sum_pos - sum_neg
				result.err = del_pos + del_neg
				result.err += 2.0 * gsl.Float64Eps * (sum_pos + sum_neg)
				result.err += 2.0 * gsl.Float64Eps * (2.0*math.Sqrt(k) + 1.0) * math.Abs(result.val)
				return err.ERROR("error", err.EMAXITER)
			}

			k += 1.0
		}

		result.val = sum_pos - sum_neg
		result.err = del_pos + del_neg
		result.err += 2.0 * gsl.Float64Eps * (sum_pos + sum_neg)
		result.err += 2.0 * gsl.Float64Eps * (2.0*math.Sqrt(k) + 1.0) * math.Abs(result.val)

		return nil
	}
}

/* Luke's rational approximation. The most accesible
 * discussion is in [Kolbig, CPC 23, 51 (1981)].
 * The convergence is supposedly guaranteed for x < 0.
 * You have to read Luke's books to see this and other
 * results. Unfortunately, the stability is not so
 * clear to me, although it seems very efficient when
 * it works.
 */
func hyperg_2F1_luke(a, b, c, xin float64, result *Result) err.GSLError {
	var (
		stat_iter err.GSLError
		RECUR_BIG = 1.0e+50
		nmax      = 20000
		n         = 3
		x         = -xin
		x3        = x * x * x
		t0        = a * b / c
		t1        = (a + 1.0) * (b + 1.0) / (2.0 * c)
		t2        = (a + 2.0) * (b + 2.0) / (2.0 * (c + 1.0))
		F         = 1.0
		prec      float64
		Bnm3      = 1.0                                            /* B0 */
		Bnm2      = 1.0 + t1*x                                     /* B1 */
		Bnm1      = 1.0 + t2*x*(1.0+t1/3.0*x)                      /* B2 */
		Anm3      = 1.0                                            /* A0 */
		Anm2      = Bnm2 - t0*x                                    /* A1 */
		Anm1      = Bnm1 - t0*(1.0+t2*x)*x + t0*t1*(c/(c+1.0))*x*x /* A2 */
	)

	for {
		var (
			fn    = float64(n)
			npam1 = fn + a - 1
			npbm1 = fn + b - 1
			npcm1 = fn + c - 1
			npam2 = fn + a - 2
			npbm2 = fn + b - 2
			npcm2 = fn + c - 2
			tnm1  = 2*fn - 1
			tnm3  = 2*fn - 3
			tnm5  = 2*fn - 5
			n2    = fn * fn
			F1    = (3.0*n2 + (a+b-6)*fn + 2 - a*b - 2*(a+b)) / (2 * tnm3 * npcm1)
			F2    = -(3.0*n2 - (a+b+6)*fn + 2 - a*b) * npam1 * npbm1 / (4 * tnm1 * tnm3 * npcm2 * npcm1)
			F3    = (npam2 * npam1 * npbm2 * npbm1 * (fn - a - 2) * (fn - b - 2)) / (8 * tnm3 * tnm3 * tnm5 * (fn + c - 3) * npcm2 * npcm1)
			E     = -npam1 * npbm1 * (fn - c - 1) / (2 * tnm3 * npcm2 * npcm1)
			An    = (1.0+F1*x)*Anm1 + (E+F2*x)*x*Anm2 + F3*x3*Anm3
			Bn    = (1.0+F1*x)*Bnm1 + (E+F2*x)*x*Bnm2 + F3*x3*Bnm3
			r     = An / Bn
		)

		prec = math.Abs((F - r) / F)
		F = r

		if prec < gsl.Float64Eps || n > nmax {
			break
		}

		if math.Abs(An) > RECUR_BIG || math.Abs(Bn) > RECUR_BIG {
			An /= RECUR_BIG
			Bn /= RECUR_BIG
			Anm1 /= RECUR_BIG
			Bnm1 /= RECUR_BIG
			Anm2 /= RECUR_BIG
			Bnm2 /= RECUR_BIG
			Anm3 /= RECUR_BIG
			Bnm3 /= RECUR_BIG
		} else if math.Abs(An) < 1.0/RECUR_BIG || math.Abs(Bn) < 1.0/RECUR_BIG {
			An *= RECUR_BIG
			Bn *= RECUR_BIG
			Anm1 *= RECUR_BIG
			Bnm1 *= RECUR_BIG
			Anm2 *= RECUR_BIG
			Bnm2 *= RECUR_BIG
			Anm3 *= RECUR_BIG
			Bnm3 *= RECUR_BIG
		}

		n++
		Bnm3 = Bnm2
		Bnm2 = Bnm1
		Bnm1 = Bn
		Anm3 = Anm2
		Anm2 = Anm1
		Anm1 = An
	}

	result.val = F
	result.err = 2.0 * math.Abs(prec*F)
	result.err += 2.0 * gsl.Float64Eps * (float64(n) + 1.0) * math.Abs(F)

	/* FIXME: just a hack: there's a lot of shit going on here */
	result.err *= 8.0 * (math.Abs(a) + math.Abs(b) + 1.0)

	if n >= nmax {
		stat_iter = err.MaxIteration()
	}

	return stat_iter
}

/* Luke's rational approximation for the
 * case a = aR + i aI, b = aR - i aI.
 */
func hyperg_2F1_conj_luke(aR, aI, c, xin float64, result *Result) err.GSLError {
	var (
		stat_iter err.GSLError
		RECUR_BIG = 1.0e+50
		nmax      = 10000
		n         = 3
		x         = -xin
		x3        = x * x * x
		atimesb   = aR*aR + aI*aI
		apb       = 2.0 * aR
		t0        = atimesb / c
		t1        = (atimesb + apb + 1.0) / (2.0 * c)
		t2        = (atimesb + 2.0*apb + 4.0) / (2.0 * (c + 1.0))
		F         = 1.0
		prec      float64
		Bnm3      = 1.0                                            /* B0 */
		Bnm2      = 1.0 + t1*x                                     /* B1 */
		Bnm1      = 1.0 + t2*x*(1.0+t1/3.0*x)                      /* B2 */
		Anm3      = 1.0                                            /* A0 */
		Anm2      = Bnm2 - t0*x                                    /* A1 */
		Anm1      = Bnm1 - t0*(1.0+t2*x)*x + t0*t1*(c/(c+1.0))*x*x /* A2 */
	)
	for {
		var (
			fn          = float64(n)
			nm1         = fn - 1
			nm2         = fn - 2
			npam1_npbm1 = atimesb + nm1*apb + nm1*nm1
			npam2_npbm2 = atimesb + nm2*apb + nm2*nm2
			npcm1       = nm1 + c
			npcm2       = nm2 + c
			tnm1        = 2*fn - 1
			tnm3        = 2*fn - 3
			tnm5        = 2*fn - 5
			n2          = fn * fn
			F1          = (3.0*n2 + (apb-6)*fn + 2 - atimesb - 2*apb) / (2 * tnm3 * npcm1)
			F2          = -(3.0*n2 - (apb+6)*fn + 2 - atimesb) * npam1_npbm1 / (4 * tnm1 * tnm3 * npcm2 * npcm1)
			F3          = (npam2_npbm2 * npam1_npbm1 * (nm2*nm2 - nm2*apb + atimesb)) / (8 * tnm3 * tnm3 * tnm5 * (fn + c - 3) * npcm2 * npcm1)
			E           = -npam1_npbm1 * (fn - c - 1) / (2 * tnm3 * npcm2 * npcm1)

			An = (1.0+F1*x)*Anm1 + (E+F2*x)*x*Anm2 + F3*x3*Anm3
			Bn = (1.0+F1*x)*Bnm1 + (E+F2*x)*x*Bnm2 + F3*x3*Bnm3
			r  = An / Bn
		)

		prec = math.Abs(F-r) / math.Abs(F)
		F = r

		if prec < gsl.Float64Eps || n > nmax {
			break
		}

		if math.Abs(An) > RECUR_BIG || math.Abs(Bn) > RECUR_BIG {
			An /= RECUR_BIG
			Bn /= RECUR_BIG
			Anm1 /= RECUR_BIG
			Bnm1 /= RECUR_BIG
			Anm2 /= RECUR_BIG
			Bnm2 /= RECUR_BIG
			Anm3 /= RECUR_BIG
			Bnm3 /= RECUR_BIG
		} else if math.Abs(An) < 1.0/RECUR_BIG || math.Abs(Bn) < 1.0/RECUR_BIG {
			An *= RECUR_BIG
			Bn *= RECUR_BIG
			Anm1 *= RECUR_BIG
			Bnm1 *= RECUR_BIG
			Anm2 *= RECUR_BIG
			Bnm2 *= RECUR_BIG
			Anm3 *= RECUR_BIG
			Bnm3 *= RECUR_BIG
		}

		n++
		Bnm3 = Bnm2
		Bnm2 = Bnm1
		Bnm1 = Bn
		Anm3 = Anm2
		Anm2 = Anm1
		Anm1 = An
	}

	result.val = F
	result.err = 2.0 * math.Abs(prec*F)
	result.err += 2.0 * gsl.Float64Eps * (float64(n) + 1.0) * math.Abs(F)

	/* FIXME: see above */
	result.err *= 8.0 * (math.Abs(aR) + math.Abs(aI) + 1.0)

	if n >= nmax {
		stat_iter = err.MaxIteration()
	}

	return stat_iter
}

/* Do the reflection described in [Moshier, p. 334].
 * Assumes a,b,c != neg integer.
 */
func hyperg_2F1_reflect(a, b, c, x float64, result *Result) err.GSLError {
	var (
		d         = c - a - b
		intd      = int(math.Floor(d + 0.5))
		d_integer = (math.Abs(d-float64(intd)) < locEPS)
	)

	if d_integer {
		ln_omx := math.Log(1.0 - x)
		ad := math.Abs(d)
		var stat_F2 err.GSLError
		sgn_2 := 1.
		F1, F2 := new(Result), new(Result)
		var d1, d2 float64
		lng_c, lng_ad2, lng_bd2 := new(Result), new(Result), new(Result)
		var stat_ad2, stat_bd2 err.GSLError

		if d >= 0.0 {
			d1 = d
			d2 = 0.0
		} else {
			d1 = 0.0
			d2 = d
		}

		stat_ad2 = Lngamma_e(a+d2, lng_ad2)
		stat_bd2 = Lngamma_e(b+d2, lng_bd2)
		Lngamma_e(c, lng_c) // stat_c

		/* Evaluate F1.
		 */
		if ad < gsl.Float64Eps {
			/* d = 0 */
			F1.val = 0.0
			F1.err = 0.0
		} else {
			lng_ad, lng_ad1, lng_bd1 := new(Result), new(Result), new(Result)
			stat_ad := Lngamma_e(ad, lng_ad)
			stat_ad1 := Lngamma_e(a+d1, lng_ad1)
			stat_bd1 := Lngamma_e(b+d1, lng_bd1)

			if stat_ad1 == nil && stat_bd1 == nil && stat_ad == nil {
				/* Gamma functions in the denominator are ok.
				 * Proceed with evaluation.
				 */
				var (
					i           int
					sum1        = 1.0
					term        = 1.0
					ln_pre1_val = lng_ad.val + lng_c.val + d2*ln_omx - lng_ad1.val - lng_bd1.val
					ln_pre1_err = lng_ad.err + lng_c.err + lng_ad1.err + lng_bd1.err + gsl.Float64Eps*math.Abs(ln_pre1_val)
					stat_e      err.GSLError
				)

				/* Do F1 sum.
				 */
				for i = 1; float64(i) < ad; i++ {
					j := float64(i) - 1
					term *= (a + d2 + j) * (b + d2 + j) / (1.0 + d2 + j) / float64(i) * (1.0 - x)
					sum1 += term
				}

				stat_e = Exp_mult_err_e(ln_pre1_val, ln_pre1_err, sum1, gsl.Float64Eps*math.Abs(sum1), F1)

				if stat_e != nil {
					if stat_e.Status() == err.EOVRFLW {
						return OverflowError(result)
					}
				}

			} else {
				/* Gamma functions in the denominator were not ok.
				 * So the F1 term is zero.
				 */
				F1.val = 0.0
				F1.err = 0.0
			}
		} /* end F1 evaluation */

		/* Evaluate F2.
		 */
		if stat_ad2 == nil && stat_bd2 == nil {
			/* Gamma functions in the denominator are ok.
			 * Proceed with evaluation.
			 */
			var (
				maxiter                     = 2000
				psi_1                       = -gsl.Euler
				psi_1pd, psi_apd1, psi_bpd1 Result
				stat_1pd                    = Psi_e(1.0+ad, &psi_1pd)
				stat_apd1                   = Psi_e(a+d1, &psi_apd1)
				stat_bpd1                   = Psi_e(b+d1, &psi_bpd1)
				stat_dall                   = err.ErrorSelect(stat_1pd, stat_apd1, stat_bpd1)
				psi_val                     = psi_1 + psi_1pd.val - psi_apd1.val - psi_bpd1.val - ln_omx
				psi_err                     = psi_1pd.err + psi_apd1.err + psi_bpd1.err + gsl.Float64Eps*math.Abs(psi_val)
				fact                        = 1.0
				sum2_val                    = psi_val
				sum2_err                    = psi_err
				ln_pre2_val                 = lng_c.val + d1*ln_omx - lng_ad2.val - lng_bd2.val
				ln_pre2_err                 = lng_c.err + lng_ad2.err + lng_bd2.err + gsl.Float64Eps*math.Abs(ln_pre2_val)
				stat_e                      err.GSLError
				j                           int
			)

			/* Do F2 sum.
			 */
			for j = 1; j < maxiter; j++ {
				fj := float64(j)
				/* values for psi functions use recurrence; Abramowitz+Stegun 6.3.5 */
				term1 := 1.0/fj + 1.0/(ad+fj)
				term2 := 1.0/(a+d1+fj-1.0) + 1.0/(b+d1+fj-1.0)
				delta := 0.0
				psi_val += term1 - term2
				psi_err += gsl.Float64Eps * (math.Abs(term1) + math.Abs(term2))
				fact *= (a + d1 + fj - 1.0) * (b + d1 + fj - 1.0) / ((ad + fj) * fj) * (1.0 - x)
				delta = fact * psi_val
				sum2_val += delta
				sum2_err += math.Abs(fact*psi_err) + gsl.Float64Eps*math.Abs(delta)
				if math.Abs(delta) < gsl.Float64Eps*math.Abs(sum2_val) {
					break
				}
			}

			if j == maxiter {
				stat_F2 = err.MaxIteration()
			}

			if sum2_val == 0.0 {
				F2.val = 0.0
				F2.err = 0.0
			} else {
				stat_e = Exp_mult_err_e(ln_pre2_val, ln_pre2_err, sum2_val, sum2_err, F2)
				if stat_e != nil {
					if stat_e.Status() == err.EOVRFLW {
						result.val = 0.0
						result.err = 0.0
						return err.ERROR("error", err.EOVRFLW)
					}
				}

			}
			stat_F2 = err.ErrorSelect(stat_F2, stat_dall)
		} else {
			/* Gamma functions in the denominator not ok.
			 * So the F2 term is zero.
			 */
			F2.val = 0.0
			F2.err = 0.0
		} /* end F2 evaluation */

		if gsl.IsOdd(intd) {
			sgn_2 = -1.0
		}

		result.val = F1.val + sgn_2*F2.val
		result.err = F1.err + F2.err
		result.err += 2.0 * gsl.Float64Eps * (math.Abs(F1.val) + math.Abs(F2.val))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_F2
	} else {
		/* d not an integer */

		pre1, pre2 := new(Result), new(Result)
		var sgn1, sgn2 float64
		F1, F2 := new(Result), new(Result)
		var status_F1, status_F2 err.GSLError

		/* These gamma functions appear in the denominator, so we
		 * catch their harmless domain errors and set the terms to zero.
		 */
		ln_g1ca, ln_g1cb, ln_g2a, ln_g2b := new(Result), new(Result), new(Result), new(Result)
		var sgn_g1ca, sgn_g1cb, sgn_g2a, sgn_g2b float64
		stat_1ca := Lngamma_sgn_e(c-a, ln_g1ca, &sgn_g1ca)
		stat_1cb := Lngamma_sgn_e(c-b, ln_g1cb, &sgn_g1cb)
		stat_2a := Lngamma_sgn_e(a, ln_g2a, &sgn_g2a)
		stat_2b := Lngamma_sgn_e(b, ln_g2b, &sgn_g2b)
		ok1 := (stat_1ca == nil && stat_1cb == nil)
		ok2 := (stat_2a == nil && stat_2b == nil)

		ln_gc, ln_gd, ln_gmd := new(Result), new(Result), new(Result)
		var sgn_gc, sgn_gd, sgn_gmd float64
		Lngamma_sgn_e(c, ln_gc, &sgn_gc)
		Lngamma_sgn_e(d, ln_gd, &sgn_gd)
		Lngamma_sgn_e(-d, ln_gmd, &sgn_gmd)

		sgn1 = sgn_gc * sgn_gd * sgn_g1ca * sgn_g1cb
		sgn2 = sgn_gc * sgn_gmd * sgn_g2a * sgn_g2b

		if ok1 && ok2 {
			ln_pre1_val := ln_gc.val + ln_gd.val - ln_g1ca.val - ln_g1cb.val
			ln_pre2_val := ln_gc.val + ln_gmd.val - ln_g2a.val - ln_g2b.val + d*math.Log(1.0-x)
			ln_pre1_err := ln_gc.err + ln_gd.err + ln_g1ca.err + ln_g1cb.err
			ln_pre2_err := ln_gc.err + ln_gmd.err + ln_g2a.err + ln_g2b.err
			if ln_pre1_val < gsl.LnMaxFloat64 && ln_pre2_val < gsl.LnMaxFloat64 {
				Exp_err_e(ln_pre1_val, ln_pre1_err, pre1)
				Exp_err_e(ln_pre2_val, ln_pre2_err, pre2)
				pre1.val *= sgn1
				pre2.val *= sgn2
			} else {
				return OverflowError(result)
			}
		} else if ok1 && !ok2 {
			ln_pre1_val := ln_gc.val + ln_gd.val - ln_g1ca.val - ln_g1cb.val
			ln_pre1_err := ln_gc.err + ln_gd.err + ln_g1ca.err + ln_g1cb.err
			if ln_pre1_val < gsl.LnMaxFloat64 {
				Exp_err_e(ln_pre1_val, ln_pre1_err, pre1)
				pre1.val *= sgn1
				pre2.val = 0.0
				pre2.err = 0.0
			} else {
				return OverflowError(result)
			}
		} else if !ok1 && ok2 {
			ln_pre2_val := ln_gc.val + ln_gmd.val - ln_g2a.val - ln_g2b.val + d*math.Log(1.0-x)
			ln_pre2_err := ln_gc.err + ln_gmd.err + ln_g2a.err + ln_g2b.err
			if ln_pre2_val < gsl.LnMaxFloat64 {
				pre1.val = 0.0
				pre1.err = 0.0
				Exp_err_e(ln_pre2_val, ln_pre2_err, pre2)
				pre2.val *= sgn2
			} else {
				return OverflowError(result)
			}
		} else {
			pre1.val = 0.0
			pre2.val = 0.0
			return UnderflowError(result)
		}

		status_F1 = hyperg_2F1_series(a, b, 1.0-d, 1.0-x, F1)
		status_F2 = hyperg_2F1_series(c-a, c-b, 1.0+d, 1.0-x, F2)

		result.val = pre1.val*F1.val + pre2.val*F2.val
		result.err = math.Abs(pre1.val*F1.err) + math.Abs(pre2.val*F2.err)
		result.err += math.Abs(pre1.err*F1.val) + math.Abs(pre2.err*F2.val)
		result.err += 2.0 * gsl.Float64Eps * (math.Abs(pre1.val*F1.val) + math.Abs(pre2.val*F2.val))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)

		if status_F1 != nil {
			return status_F1
		}

		if status_F2 != nil {
			return status_F2
		}

		return nil
	}
}

func pow_omx(x, p float64, result *Result) err.GSLError {
	var ln_omx, ln_result float64
	if math.Abs(x) < gsl.Root5Float64Eps {
		ln_omx = -x * (1.0 + x*(1.0/2.0+x*(1.0/3.0+x/4.0+x*x/5.0)))
	} else {
		ln_omx = math.Log(1.0 - x)
	}
	ln_result = p * ln_omx
	return Exp_err_e(ln_result, gsl.Float64Eps*math.Abs(ln_result), result)
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Hyperg_2F1_e(a, b, c, x float64, result *Result) err.GSLError {
	var (
		d             = c - a - b
		rinta         = math.Floor(a + 0.5)
		rintb         = math.Floor(b + 0.5)
		rintc         = math.Floor(c + 0.5)
		a_neg_integer = (a < 0.0 && math.Abs(a-rinta) < locEPS)
		b_neg_integer = (b < 0.0 && math.Abs(b-rintb) < locEPS)
		c_neg_integer = (c < 0.0 && math.Abs(c-rintc) < locEPS)
	)

	result.val = 0.0
	result.err = 0.0

	/* Handle x == 1.0 RJM */

	if math.Abs(x-1.0) < locEPS && (c-a-b) > 0 && c != 0 && !c_neg_integer {
		lngamc, lngamcab, lngamca, lngamcb := new(Result), new(Result), new(Result), new(Result)
		var lngamc_sgn, lngamca_sgn, lngamcb_sgn float64
		var status err.GSLError
		stat1 := Lngamma_sgn_e(c, lngamc, &lngamc_sgn)
		stat2 := Lngamma_e(c-a-b, lngamcab)
		stat3 := Lngamma_sgn_e(c-a, lngamca, &lngamca_sgn)
		stat4 := Lngamma_sgn_e(c-b, lngamcb, &lngamcb_sgn)

		if stat1 != nil || stat2 != nil || stat3 != nil || stat4 != nil {
			return DomainError(result)
		}

		status = Exp_err_e(lngamc.val+lngamcab.val-lngamca.val-lngamcb.val, lngamc.err+lngamcab.err+lngamca.err+lngamcb.err, result)

		result.val *= lngamc_sgn / (lngamca_sgn * lngamcb_sgn)
		return status
	}

	if x < -1.0 || 1.0 <= x {
		return DomainError(result)
	}

	if c_neg_integer {
		/* If c is a negative integer, then either a or b must be a
		   negative integer of smaller magnitude than c to ensure
		   cancellation of the series. */
		if !(a_neg_integer && a > c+0.1) && !(b_neg_integer && b > c+0.1) {
			return DomainError(result)
		}
	}

	if math.Abs(c-b) < locEPS || math.Abs(c-a) < locEPS {
		return pow_omx(x, d, result) /* (1-x)^(c-a-b) */
	}

	if a >= 0.0 && b >= 0.0 && c >= 0.0 && x >= 0.0 && x < 0.995 {
		/* Series has all positive definite
		 * terms and x is not close to 1.
		 */
		return hyperg_2F1_series(a, b, c, x, result)
	}

	if math.Abs(a) < 10.0 && math.Abs(b) < 10.0 {
		/* a and b are not too large, so we attempt
		 * variations on the series summation.
		 */
		if a_neg_integer {
			return hyperg_2F1_series(rinta, b, c, x, result)
		}
		if b_neg_integer {
			return hyperg_2F1_series(a, rintb, c, x, result)
		}

		if x < -0.25 {
			return hyperg_2F1_luke(a, b, c, x, result)
		} else if x < 0.5 {
			return hyperg_2F1_series(a, b, c, x, result)
		} else {
			if math.Abs(c) > 10.0 {
				return hyperg_2F1_series(a, b, c, x, result)
			} else {
				return hyperg_2F1_reflect(a, b, c, x, result)
			}
		}
	} else {
		/* Either a or b or both large.
		 * Introduce some new variables ap,bp so that bp is
		 * the larger in magnitude.
		 */
		var ap, bp float64
		if math.Abs(a) > math.Abs(b) {
			bp = a
			ap = b
		} else {
			bp = b
			ap = a
		}

		if x < 0.0 {
			/* What the hell, maybe Luke will converge.
			 */
			return hyperg_2F1_luke(a, b, c, x, result)
		}

		if math.Max(math.Abs(ap), 1.0)*math.Abs(bp)*math.Abs(x) < 2.0*math.Abs(c) {
			/* If c is large enough or x is small enough,
			 * we can attempt the series anyway.
			 */
			return hyperg_2F1_series(a, b, c, x, result)
		}

		if math.Abs(bp*bp*x*x) < 0.001*math.Abs(bp) && math.Abs(ap) < 10.0 {
			/* The famous but nearly worthless "large b" asymptotic.
			 */
			stat := Hyperg_1F1_e(ap, c, bp*x, result)
			result.err = 0.001 * math.Abs(result.val)
			return stat
		}

		/* We give up. */
		result.val = 0.0
		result.err = 0.0
		return err.ERROR("error", err.EUNIMPL)
	}
}

func Hyperg_2F1_conj_e(aR, aI, c, x float64, result *Result) err.GSLError {
	var (
		ax            = math.Abs(x)
		rintc         = math.Floor(c + 0.5)
		c_neg_integer = (c < 0.0 && math.Abs(c-rintc) < locEPS)
	)

	result.val = 0.0
	result.err = 0.0

	if ax >= 1.0 || c_neg_integer || c == 0.0 {
		return DomainError(result)
	}

	if (ax < 0.25 && math.Abs(aR) < 20.0 && math.Abs(aI) < 20.0) || (c > 0.0 && x > 0.0) {
		return hyperg_2F1_conj_series(aR, aI, c, x, result)
	} else if math.Abs(aR) < 10.0 && math.Abs(aI) < 10.0 {
		if x < -0.25 {
			return hyperg_2F1_conj_luke(aR, aI, c, x, result)
		} else {
			return hyperg_2F1_conj_series(aR, aI, c, x, result)
		}
	} else {
		if x < 0.0 {
			/* What the hell, maybe Luke will converge.
			 */
			return hyperg_2F1_conj_luke(aR, aI, c, x, result)
		}

		/* Give up. */
		result.val = 0.0
		result.err = 0.0
		return err.ERROR("error", err.EUNIMPL)
	}
}

func Hyperg_2F1_renorm_e(a, b, c, x float64, result *Result) err.GSLError {
	var (
		rinta         = math.Floor(a + 0.5)
		rintb         = math.Floor(b + 0.5)
		rintc         = math.Floor(c + 0.5)
		a_neg_integer = (a < 0.0 && math.Abs(a-rinta) < locEPS)
		b_neg_integer = (b < 0.0 && math.Abs(b-rintb) < locEPS)
		c_neg_integer = (c < 0.0 && math.Abs(c-rintc) < locEPS)
	)

	if c_neg_integer {
		if (a_neg_integer && a > c+0.1) || (b_neg_integer && b > c+0.1) {
			/* 2F1 terminates early */
			result.val = 0.0
			result.err = 0.0
			return nil
		} else {
			/* 2F1 does not terminate early enough, so something survives */
			/* [Abramowitz+Stegun, 15.1.2] */
			g1, g2, g3, g4, g5 := new(Result), new(Result), new(Result), new(Result), new(Result)
			var s1, s2, s3, s4, s5 float64
			var stat err.GSLError
			stat = Lngamma_sgn_e(a-c+1, g1, &s1)
			stat = Lngamma_sgn_e(b-c+1, g2, &s2)
			stat = Lngamma_sgn_e(a, g3, &s3)
			stat = Lngamma_sgn_e(b, g4, &s4)
			stat = Lngamma_sgn_e(-c+2, g5, &s5)
			if stat != nil {
				return DomainError(result)
			} else {
				F := new(Result)
				stat_F := Hyperg_2F1_e(a-c+1, b-c+1, -c+2, x, F)
				ln_pre_val := g1.val + g2.val - g3.val - g4.val - g5.val
				ln_pre_err := g1.err + g2.err + g3.err + g4.err + g5.err
				sg := s1 * s2 * s3 * s4 * s5
				stat_e := Exp_mult_err_e(ln_pre_val, ln_pre_err, sg*F.val, F.err, result)
				return err.ErrorSelect(stat_e, stat_F)
			}
		}
	} else {
		/* generic c */
		F, lng := new(Result), new(Result)
		var sgn float64
		stat_g := Lngamma_sgn_e(c, lng, &sgn)
		stat_F := Hyperg_2F1_e(a, b, c, x, F)
		stat_e := Exp_mult_err_e(-lng.val, lng.err, sgn*F.val, F.err, result)
		return err.ErrorSelect(stat_e, stat_F, stat_g)
	}
}

func Hyperg_2F1_conj_renorm_e(aR, aI, c, x float64, result *Result) err.GSLError {
	var (
		rintc         = math.Floor(c + 0.5)
		rinta         = math.Floor(aR + 0.5)
		a_neg_integer = (aR < 0.0 && math.Abs(aR-rinta) < locEPS && aI == 0.0)
		c_neg_integer = (c < 0.0 && math.Abs(c-rintc) < locEPS)
	)

	if c_neg_integer {
		if a_neg_integer && aR > c+0.1 {
			/* 2F1 terminates early */
			result.val = 0.0
			result.err = 0.0
			return nil
		} else {
			/* 2F1 does not terminate early enough, so something survives */
			/* [Abramowitz+Stegun, 15.1.2] */
			g1, g2, g3, a1, a2 := new(Result), new(Result), new(Result), new(Result), new(Result)
			var stat err.GSLError
			stat = Lngamma_complex_e(aR-c+1, aI, g1, a1)
			stat = Lngamma_complex_e(aR, aI, g2, a2)
			stat = Lngamma_e(-c+2.0, g3)
			if stat != nil {
				return DomainError(result)
			} else {
				F := new(Result)
				stat_F := Hyperg_2F1_conj_e(aR-c+1, aI, -c+2, x, F)
				ln_pre_val := 2.0*(g1.val-g2.val) - g3.val
				ln_pre_err := 2.0*(g1.err+g2.err) + g3.err
				stat_e := Exp_mult_err_e(ln_pre_val, ln_pre_err, F.val, F.err, result)
				return err.ErrorSelect(stat_e, stat_F)
			}
		}
	} else {
		/* generic c */
		F, lng := new(Result), new(Result)
		var sgn float64
		stat_g := Lngamma_sgn_e(c, lng, &sgn)
		stat_F := Hyperg_2F1_conj_e(aR, aI, c, x, F)
		stat_e := Exp_mult_err_e(-lng.val, lng.err, sgn*F.val, F.err, result)
		return err.ErrorSelect(stat_e, stat_F, stat_g)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Hyperg_2F1(a, b, c, x float64) float64 {
	result := new(Result)
	status := Hyperg_2F1_e(a, b, c, x, result)
	return EvalResult(result, status)
}

func Hyperg_2F1_conj(aR, aI, c, x float64) float64 {
	result := new(Result)
	status := Hyperg_2F1_conj_e(aR, aI, c, x, result)
	return EvalResult(result, status)
}

func Hyperg_2F1_renorm(a, b, c, x float64) float64 {
	result := new(Result)
	status := Hyperg_2F1_renorm_e(a, b, c, x, result)
	return EvalResult(result, status)
}

func Hyperg_2F1_conj_renorm(aR, aI, c, x float64) float64 {
	result := new(Result)
	status := Hyperg_2F1_conj_renorm_e(aR, aI, c, x, result)
	return EvalResult(result, status)
}
