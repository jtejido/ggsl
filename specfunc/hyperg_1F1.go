/* specfunc/hyperg_1F1.c
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
	SMALL_ANGLE        bool
	_1F1_INT_THRESHOLD = 100.0 * gsl.Float64Eps
)

/* Asymptotic result for 1F1(a, b, x)  x . -Infinity.
 * Assumes b-a != neg integer and b != neg integer.
 */
func hyperg_1F1_asymp_negx(a, b, x float64, result *Result) err.GSLError {
	lg_b, lg_bma := new(Result), new(Result)
	var sgn_b, sgn_bma float64

	stat_b := Lngamma_sgn_e(b, lg_b, &sgn_b)
	stat_bma := Lngamma_sgn_e(b-a, lg_bma, &sgn_bma)

	if stat_b == nil && stat_bma == nil {
		F := new(Result)
		stat_F := Hyperg_2F0_series_e(a, 1.0+a-b, -1.0/x, -1, F)
		if F.val != 0 {
			ln_term_val := a * math.Log(-x)
			ln_term_err := 2.0 * gsl.Float64Eps * (math.Abs(a) + math.Abs(ln_term_val))
			ln_pre_val := lg_b.val - lg_bma.val - ln_term_val
			ln_pre_err := lg_b.err + lg_bma.err + ln_term_err
			stat_e := Exp_mult_err_e(ln_pre_val, ln_pre_err, sgn_bma*sgn_b*F.val, F.err, result)
			return err.ErrorSelect(stat_e, stat_F)
		} else {
			result.val = 0.0
			result.err = 0.0
			return stat_F
		}
	} else {
		return DomainError(result)
	}
}

/* Asymptotic result for 1F1(a, b, x)  x . +Infinity
 * Assumes b != neg integer and a != neg integer
 */
func hyperg_1F1_asymp_posx(a, b, x float64, result *Result) err.GSLError {
	lg_b, lg_a := new(Result), new(Result)
	var sgn_b, sgn_a float64

	stat_b := Lngamma_sgn_e(b, lg_b, &sgn_b)
	stat_a := Lngamma_sgn_e(a, lg_a, &sgn_a)

	if stat_a == nil && stat_b == nil {
		F := new(Result)
		stat_F := Hyperg_2F0_series_e(b-a, 1.0-a, 1.0/x, -1, F)
		if stat_F == nil && F.val != 0 {
			lnx := math.Log(x)
			ln_term_val := (a - b) * lnx
			ln_term_err := 2.0*gsl.Float64Eps*(math.Abs(a)+math.Abs(b))*math.Abs(lnx) + 2.0*gsl.Float64Eps*math.Abs(a-b)
			ln_pre_val := lg_b.val - lg_a.val + ln_term_val + x
			ln_pre_err := lg_b.err + lg_a.err + ln_term_err + 2.0*gsl.Float64Eps*math.Abs(x)
			stat_e := Exp_mult_err_e(ln_pre_val, ln_pre_err,
				sgn_a*sgn_b*F.val, F.err,
				result)
			return err.ErrorSelect(stat_e, stat_F)
		} else {
			result.val = 0.0
			result.err = 0.0
			return stat_F
		}
	} else {
		return DomainError(result)
	}
}

/* Asymptotic result from Slater 4.3.7
 *
 * To get the general series, write M(a,b,x) as
 *
 *  M(a,b,x)=sum ((a)_n/(b)_n) (x^n / n!)
 *
 * and expand (b)_n in inverse powers of b as follows
 *
 * -math.Log(1/(b)_n) = sum_(k=0)^(n-1) math.Log(b+k)
 *             = n math.Log(b) + sum_(k=0)^(n-1) math.Log(1+k/b)
 *
 * Do a taylor expansion of the log in 1/b and sum the resulting terms
 * using the standard algebraic formulas for finite sums of powers of
 * k.  This should then give
 *
 * M(a,b,x) = sum_(n=0)^(inf) (a_n/n!) (x/b)^n * (1 - n(n-1)/(2b)
 *                          + (n-1)n(n+1)(3n-2)/(24b^2) + ...
 *
 * which can be summed explicitly. The trick for summing it is to take
 * derivatives of sum_(i=0)^(inf) a_n*y^n/n! = (1-y)^(-a);
 *
 * [BJG 16/01/2007]
 */

func hyperg_1F1_largebx(a, b, x float64, result *Result) err.GSLError {
	y := x / b
	f := math.Exp(-a * math.Log1p(-y))
	t1 := -((a * (a + 1.0)) / (2 * b)) * math.Pow((y/(1.0-y)), 2.0)
	t2 := (1 / (24 * b * b)) * ((a * (a + 1) * y * y) / math.Pow(1-y, 4)) * (12 + 8*(2*a+1)*y + (3*a*a-a-2)*y*y)
	t3 := (-1 / (48 * b * b * b * math.Pow(1-y, 6))) * a * ((a + 1) * ((y*((a+1)*(a*(y*(y*((y*(a-2)+16)*(a-1))+72))+96)) + 24) * math.Pow(y, 2)))
	result.val = f * (1 + t1 + t2 + t3)
	result.err = 2*math.Abs(f*t3) + 2*gsl.Float64Eps*math.Abs(result.val)
	return nil
}

/* Asymptotic result for x < 2b-4a, 2b-4a large.
 * [Abramowitz+Stegun, 13.5.21]
 *
 * assumes 0 <= x/(2b-4a) <= 1
 */
func hyperg_1F1_large2bm4a(a, b, x float64, result *Result) err.GSLError {
	eta := 2.0*b - 4.0*a
	cos2th := x / eta
	sin2th := 1.0 - cos2th
	th := math.Acos(math.Sqrt(cos2th))
	pre_h := 0.25 * gsl.Pi * gsl.Pi * eta * eta * cos2th * sin2th
	lg_b := new(Result)
	stat_lg := Lngamma_e(b, lg_b)
	t1 := 0.5 * (1.0 - b) * math.Log(0.25*x*eta)
	t2 := 0.25 * math.Log(pre_h)
	lnpre_val := lg_b.val + 0.5*x + t1 - t2
	lnpre_err := lg_b.err + 2.0*gsl.Float64Eps*(math.Abs(0.5*x)+math.Abs(t1)+math.Abs(t2))
	var s1, s2 float64
	if SMALL_ANGLE {
		eps := math.Asin(math.Sqrt(cos2th)) /* theta = pi/2 - eps */
		if math.Mod(a, 1.0) != 0.0 {
			s1 = math.Sin(a * gsl.Pi)
		}
		var eta_reduc float64
		if math.Mod(eta+1, 4.0) != 0.0 {
			eta_reduc = math.Mod(eta+1, 8.0)
		}
		phi1 := 0.25 * eta_reduc * gsl.Pi
		phi2 := 0.25 * eta * (2*eps + math.Sin(2.0*eps))
		s2 = math.Sin(phi1 - phi2)
	} else {
		s1 = math.Sin(a * gsl.Pi)
		s2 = math.Sin(0.25*eta*(2.0*th-math.Sin(2.0*th)) + 0.25*gsl.Pi)
	}
	ser_val := s1 + s2
	ser_err := 2.0 * gsl.Float64Eps * (math.Abs(s1) + math.Abs(s2))
	stat_e := Exp_mult_err_e(lnpre_val, lnpre_err, ser_val, ser_err, result)
	return err.ErrorSelect(stat_e, stat_lg)
}

/* Luke's rational approximation.
 * See [Luke, Algorithms for the Computation of Mathematical Functions, p.182]
 *
 * Like the case of the 2F1 rational approximations, these are
 * probably guaranteed to converge for x < 0, barring gross
 * numerical instability in the pre-asymptotic regime.
 */
func hyperg_1F1_luke(a, c, xin float64, result *Result) err.GSLError {
	RECUR_BIG := 1.0e+50
	nmax := 5000
	n := 3
	x := -xin
	x3 := x * x * x
	t0 := a / c
	t1 := (a + 1.0) / (2.0 * c)
	t2 := (a + 2.0) / (2.0 * (c + 1.0))
	F := 1.0
	var prec float64

	Bnm3 := 1.0                       /* B0 */
	Bnm2 := 1.0 + t1*x                /* B1 */
	Bnm1 := 1.0 + t2*x*(1.0+t1/3.0*x) /* B2 */

	Anm3 := 1.0                                            /* A0 */
	Anm2 := Bnm2 - t0*x                                    /* A1 */
	Anm1 := Bnm1 - t0*(1.0+t2*x)*x + t0*t1*(c/(c+1.0))*x*x /* A2 */

	for {
		fn := float64(n)
		npam1 := fn + a - 1
		npcm1 := fn + c - 1
		npam2 := fn + a - 2
		npcm2 := fn + c - 2
		tnm1 := 2*fn - 1
		tnm3 := 2*fn - 3
		tnm5 := 2*fn - 5
		F1 := (fn - a - 2) / (2 * tnm3 * npcm1)
		F2 := (fn + a) * npam1 / (4 * tnm1 * tnm3 * npcm2 * npcm1)
		F3 := -npam2 * npam1 * (fn - a - 2) / (8 * tnm3 * tnm3 * tnm5 * (fn + c - 3) * npcm2 * npcm1)
		E := -npam1 * (fn - c - 1) / (2 * tnm3 * npcm2 * npcm1)

		An := (1.0+F1*x)*Anm1 + (E+F2*x)*x*Anm2 + F3*x3*Anm3
		Bn := (1.0+F1*x)*Bnm1 + (E+F2*x)*x*Bnm2 + F3*x3*Bnm3
		r := An / Bn

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
	result.err = 2.0 * math.Abs(F*prec)
	result.err += 2.0 * gsl.Float64Eps * (float64(n) - 1.0) * math.Abs(F)

	return nil
}

/* Series for 1F1(1,b,x)
 * b > 0
 */
func hyperg_1F1_1_series(b, x float64, result *Result) err.GSLError {
	sum_val := 1.0
	sum_err := 0.0
	term := 1.0
	n := 1.0
	for math.Abs(term/sum_val) > 0.25*gsl.Float64Eps {
		term *= x / (b + n - 1)
		sum_val += term
		sum_err += 8.0*gsl.Float64Eps*math.Abs(term) + gsl.Float64Eps*math.Abs(sum_val)
		n += 1.0
	}
	result.val = sum_val
	result.err = sum_err
	result.err += 2.0 * math.Abs(term)
	return nil
}

/* 1F1(1,b,x)
 * b >= 1, b integer
 */
func hyperg_1F1_1_int(b int, x float64, result *Result) err.GSLError {
	if b < 1 {
		return DomainError(result)
	} else if b == 1 {
		return Exp_e(x, result)
	} else if b == 2 {
		return Exprel_e(x, result)
	} else if b == 3 {
		return Exprel_2_e(x, result)
	}

	return Exprel_n_e(b-1, x, result)

}

/* 1F1(1,b,x)
 * b >=1, b real
 *
 * checked OK: [GJ] Thu Oct  1 16:46:35 MDT 1998
 */
func hyperg_1F1_1(b, x float64, result *Result) err.GSLError {
	ax := math.Abs(x)
	ib := math.Floor(b + 0.1)

	if b < 1.0 {
		return DomainError(result)
	} else if b == 1.0 {
		return Exp_e(x, result)
	} else if b >= 1.4*ax {
		return hyperg_1F1_1_series(b, x, result)
	} else if math.Abs(b-ib) < _1F1_INT_THRESHOLD && ib < math.MaxInt32 {
		return hyperg_1F1_1_int(int(ib), x, result)
	} else if x > 0.0 {
		if x > 100.0 && b < 0.75*x {
			return hyperg_1F1_asymp_posx(1.0, b, x, result)
		} else if b < 1.0e+05 {
			/* Recurse backward on b, from a
			 * chosen offset point. For x > 0,
			 * which holds here, this should
			 * be a stable direction.
			 */
			off := math.Ceil(1.4*x-b) + 1.0
			bp := b + off
			M := new(Result)
			stat_s := hyperg_1F1_1_series(bp, x, M)
			err_rat := M.err / math.Abs(M.val)
			for bp > b+0.1 {
				/* M(1,b-1) = x/(b-1) M(1,b) + 1 */
				bp -= 1.0
				M.val = 1.0 + x/bp*M.val
			}
			result.val = M.val
			result.err = err_rat * math.Abs(M.val)
			result.err += 2.0 * gsl.Float64Eps * (math.Abs(off) + 1.0) * math.Abs(M.val)
			return stat_s
		} else if math.Abs(x) < math.Abs(b) && math.Abs(x) < math.Sqrt(math.Abs(b))*math.Abs(b-x) {
			return hyperg_1F1_largebx(1.0, b, x, result)
		} else if math.Abs(x) > math.Abs(b) {
			return hyperg_1F1_1_series(b, x, result)
		} else {
			return hyperg_1F1_large2bm4a(1.0, b, x, result)
		}
	} else {
		/* x <= 0 and b not large compared to |x|
		 */
		if ax < 10.0 && b < 10.0 {
			return hyperg_1F1_1_series(b, x, result)
		} else if ax >= 100.0 && math.Max(math.Abs(2.0-b), 1.0) < 0.99*ax {
			return hyperg_1F1_asymp_negx(1.0, b, x, result)
		} else {
			return hyperg_1F1_luke(1.0, b, x, result)
		}
	}
}

/* 1F1(a,b,x)/Gamma(b) for b.0
 * [limit of Abramowitz+Stegun 13.3.7]
 */
func hyperg_1F1_renorm_b0(a, x float64, result *Result) err.GSLError {
	eta := a * x
	if eta > 0.0 {
		root_eta := math.Sqrt(eta)
		I1_scaled := new(Result)
		stat_I := Bessel_I1_scaled_e(2.0*root_eta, I1_scaled)
		if I1_scaled.val <= 0.0 {
			result.val = 0.0
			result.err = 0.0
			return err.ErrorSelect(stat_I, err.Domain())
		} else {
			/* Note that 13.3.7 contains higher terms which are zeroth order
			   in b.  These make a non-negligible contribution to the sum.
			   With the first correction term, the I1 above is replaced by
			   I1 + (2/3)*a*(x/(4a))**(3/2)*I2(2*root_eta).  We will add
			   this as part of the result and error estimate. */

			corr1 := (2.0 / 3.0) * a * math.Pow(x/(4.0*a), 1.5) * Bessel_In_scaled(2, 2.0*root_eta)

			lnr_val := 0.5*x + 0.5*math.Log(eta) + math.Abs(2.0*root_eta) + math.Log(I1_scaled.val+corr1)
			lnr_err := gsl.Float64Eps*(1.5*math.Abs(x)+1.0) + math.Abs((I1_scaled.err+corr1)/I1_scaled.val)
			return Exp_err_e(lnr_val, lnr_err, result)
		}
	} else if eta == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		/* eta < 0 */
		root_eta := math.Sqrt(-eta)
		J1 := new(Result)
		stat_J := Bessel_J1_e(2.0*root_eta, J1)
		if J1.val <= 0.0 {
			result.val = 0.0
			result.err = 0.0
			return err.ErrorSelect(stat_J, err.Domain())
		} else {
			t1 := 0.5 * x
			t2 := 0.5 * math.Log(-eta)
			t3 := math.Abs(x)
			t4 := math.Log(J1.val)
			lnr_val := t1 + t2 + t3 + t4
			lnr_err := gsl.Float64Eps*(1.5*math.Abs(x)+1.0) + math.Abs(J1.err/J1.val)
			ex := new(Result)
			stat_e := Exp_err_e(lnr_val, lnr_err, ex)
			result.val = -ex.val
			result.err = ex.err
			return stat_e
		}
	}

}

/* 1F1'(a,b,x)/1F1(a,b,x)
 * Uses Gautschi's series transformation of the
 * continued fraction. This is apparently the best
 * method for getting this ratio in the stable region.
 * The convergence is monotone and supergeometric
 * when b > x.
 * Assumes a >= -1.
 */
func hyperg_1F1_CF1_p_ser(a, b, x float64, result *float64) err.GSLError {
	if a == 0.0 {
		*result = 0.0
		return nil
	} else {
		maxiter := 5000
		sum := 1.0
		pk := 1.0
		rhok := 0.0
		var k int
		for k = 1; k < maxiter; k++ {
			fk := float64(k)
			ak := (a + fk) * x / ((b - x + fk - 1.0) * (b - x + fk))
			rhok = -ak * (1.0 + rhok) / (1.0 + ak*(1.0+rhok))
			pk *= rhok
			sum += pk
			if math.Abs(pk/sum) < 2.0*gsl.Float64Eps {
				break
			}
		}
		*result = a / (b - x) * sum
		if k == maxiter {
			return err.ERROR("error", err.EMAXITER)
		}

		return nil
	}
}

/* 1F1(a,b,x)
 * |a| <= 1, b > 0
 */
func hyperg_1F1_small_a_bgt0(a, b, x float64, result *Result) err.GSLError {
	bma := b - a
	oma := 1.0 - a
	ap1mb := 1.0 + a - b
	abs_bma := math.Abs(bma)
	abs_oma := math.Abs(oma)
	abs_ap1mb := math.Abs(ap1mb)

	ax := math.Abs(x)

	if a == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if a == 1.0 && b >= 1.0 {
		return hyperg_1F1_1(b, x, result)
	} else if a == -1.0 {
		result.val = 1.0 + a/b*x
		result.err = gsl.Float64Eps * (1.0 + math.Abs(a/b*x))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if b >= 1.4*ax {
		return Hyperg_1F1_series_e(a, b, x, result)
	} else if x > 0.0 {
		if x > 100.0 && abs_bma*abs_oma < 0.5*x {
			return hyperg_1F1_asymp_posx(a, b, x, result)
		} else if b < 5.0e+06 {
			/* Recurse backward on b from
			 * a suitably high point.
			 */
			b_del := math.Ceil(1.4*x-b) + 1.0
			bp := b + b_del
			r_Mbp1, r_Mb := new(Result), new(Result)
			var Mbp1, Mb, Mbm1 float64
			stat_0 := Hyperg_1F1_series_e(a, bp+1.0, x, r_Mbp1)
			stat_1 := Hyperg_1F1_series_e(a, bp, x, r_Mb)
			err_rat := math.Abs(r_Mbp1.err/r_Mbp1.val) + math.Abs(r_Mb.err/r_Mb.val)
			Mbp1 = r_Mbp1.val
			Mb = r_Mb.val
			for bp > b+0.1 {
				/* Do backward recursion. */
				Mbm1 = ((x+bp-1.0)*Mb - x*(bp-a)/bp*Mbp1) / (bp - 1.0)
				bp -= 1.0
				Mbp1 = Mb
				Mb = Mbm1
			}
			result.val = Mb
			result.err = err_rat * (math.Abs(b_del) + 1.0) * math.Abs(Mb)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(Mb)
			return err.ErrorSelect(stat_0, stat_1)
		} else if math.Abs(x) < math.Abs(b) && math.Abs(a*x) < math.Sqrt(math.Abs(b))*math.Abs(b-x) {
			return hyperg_1F1_largebx(a, b, x, result)
		} else {
			return hyperg_1F1_large2bm4a(a, b, x, result)
		}
	} else {
		/* x < 0 and b not large compared to |x|
		 */
		if ax < 10.0 && b < 10.0 {
			return Hyperg_1F1_series_e(a, b, x, result)
		} else if ax >= 100.0 && math.Max(abs_ap1mb, 1.0) < 0.99*ax {
			return hyperg_1F1_asymp_negx(a, b, x, result)
		} else {
			return hyperg_1F1_luke(a, b, x, result)
		}
	}
}

/* 1F1(b+eps,b,x)
 * |eps|<=1, b > 0
 */
func hyperg_1F1_beps_bgt0(eps, b, x float64, result *Result) err.GSLError {
	if b > math.Abs(x) && math.Abs(eps) < gsl.SqrtFloat64Eps {
		/* If b-a is very small and x/b is not too large we can
		 * use this explicit approximation.
		 *
		 * 1F1(b+eps,b,x) = exp(ax/b) (1 - eps x^2 (v2 + v3 x + ...) + ...)
		 *
		 *   v2 = a/(2b^2(b+1))
		 *   v3 = a(b-2a)/(3b^3(b+1)(b+2))
		 *   ...
		 *
		 * See [Luke, Mathematical Functions and Their Approximations, p.292]
		 *
		 * This cannot be used for b near a negative integer or zero.
		 * Also, if x/b is large the deviation from exp(x) behaviour grows.
		 */
		a := b + eps
		exab := new(Result)
		stat_e := Exp_e(a*x/b, exab)
		v2 := a / (2.0 * b * b * (b + 1.0))
		v3 := a * (b - 2.0*a) / (3.0 * b * b * b * (b + 1.0) * (b + 2.0))
		v := v2 + v3*x
		f := (1.0 - eps*x*x*v)
		result.val = exab.val * f
		result.err = exab.err * math.Abs(f)
		result.err += math.Abs(exab.val) * gsl.Float64Eps * (1.0 + math.Abs(eps*x*x*v))
		result.err += 4.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_e
	} else {
		/* Otherwise use a Kummer transformation to reduce
		 * it to the small a case.
		 */
		Kummer_1F1 := new(Result)
		stat_K := hyperg_1F1_small_a_bgt0(-eps, b, -x, Kummer_1F1)
		if Kummer_1F1.val != 0.0 {
			stat_e := Exp_mult_err_e(x, 2.0*gsl.Float64Eps*math.Abs(x), Kummer_1F1.val, Kummer_1F1.err, result)
			return err.ErrorSelect(stat_e, stat_K)
		} else {
			result.val = 0.0
			result.err = 0.0
			return stat_K
		}
	}
}

/* 1F1(a,2a,x) = Gamma(a + 1/2) E(x) (|x|/4)^(-a+1/2) scaled_I(a-1/2,|x|/2)
 *
 * E(x) = exp(x) x > 0
 *      = 1      x < 0
 *
 * a >= 1/2
 */
func hyperg_1F1_beq2a_pos(a, x float64, result *Result) err.GSLError {
	if x == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	}

	I := new(Result)
	stat_I := Bessel_Inu_scaled_e(a-0.5, 0.5*math.Abs(x), I)
	lg := new(Result)
	stat_g := Lngamma_e(a+0.5, lg)
	ln_term := (0.5 - a) * math.Log(0.25*math.Abs(x))
	lnpre_val := lg.val + math.Max(x, 0.0) + ln_term
	lnpre_err := lg.err + gsl.Float64Eps*(math.Abs(ln_term)+math.Abs(x))
	stat_e := Exp_mult_err_e(lnpre_val, lnpre_err,
		I.val, I.err,
		result)
	return err.ErrorSelect(stat_e, stat_g, stat_I)

}

/* Handle the case of a and b both positive integers.
 * Assumes a > 0 and b > 0.
 */
func hyperg_1F1_ab_posint(a, b int, x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	if a == b {
		return Exp_e(x, result) /* 1F1(a,a,x) */
	} else if a == 1 {
		return Exprel_n_e(b-1, x, result) /* 1F1(1,b,x) */
	} else if b == a+1 {
		K := new(Result)
		stat_K := Exprel_n_e(a, -x, K) /* 1F1(1,1+a,-x) */
		stat_e := Exp_mult_err_e(x, 2.0*gsl.Float64Eps*math.Abs(x),
			K.val, K.err,
			result)
		return err.ErrorSelect(stat_e, stat_K)
	} else if a == b+1 {
		ex := new(Result)
		stat_e := Exp_e(x, ex)
		result.val = ex.val * (1.0 + x/float64(b))
		result.err = ex.err * (1.0 + x/float64(b))
		result.err += ex.val * gsl.Float64Eps * (1.0 + math.Abs(x/float64(b)))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_e
	} else if a == b+2 {
		ex := new(Result)
		stat_e := Exp_e(x, ex)
		poly := (1.0 + x/float64(b)*(2.0+x/(float64(b)+1.0)))
		result.val = ex.val * poly
		result.err = ex.err * math.Abs(poly)
		result.err += ex.val * gsl.Float64Eps * (1.0 + math.Abs(x/float64(b))*(2.0+math.Abs(x/(float64(b)+1.0))))
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_e
	} else if b == 2*a {
		return hyperg_1F1_beq2a_pos(float64(a), x, result) /* 1F1(a,2a,x) */
	} else if (b < 10 && a < 10 && ax < 5.0) || (float64(b) > float64(a)*ax) || (b > a && ax < 5.0) {
		return Hyperg_1F1_series_e(float64(a), float64(b), x, result)
	} else if b > a && float64(b) >= 2*float64(a)+x {
		/* Use the Gautschi CF series, then
		 * recurse backward to a=0 for normalization.
		 * This will work for either sign of x.
		 */
		var rap float64
		stat_CF1 := hyperg_1F1_CF1_p_ser(float64(a), float64(b), x, &rap)
		ra := 1.0 + x/float64(a)*rap
		Ma := gsl.SqrtMinFloat64
		Map1 := ra * Ma
		Mnp1 := Map1
		Mn := Ma

		for n := a; n > 0; n-- {
			fn := float64(n)
			fb := float64(b)
			Mnm1 := (fn*Mnp1 - (2*fn-fb+x)*Mn) / (fb - fn)
			Mnp1 = Mn
			Mn = Mnm1
		}
		result.val = Ma / Mn
		result.err = 2.0 * gsl.Float64Eps * (math.Abs(float64(a)) + 1.0) * math.Abs(Ma/Mn)
		return stat_CF1
	} else if b > a && float64(b) < 2*float64(a)+x && float64(b) > x {
		/* Use the Gautschi series representation of
		 * the continued fraction. Then recurse forward
		 * to the a=b line for normalization. This will
		 * work for either sign of x, although we do need
		 * to check for b > x, for when x is positive.
		 */
		var rap float64
		stat_CF1 := hyperg_1F1_CF1_p_ser(float64(a), float64(b), x, &rap)
		ra := 1.0 + x/float64(a)*rap
		ex := new(Result)

		Ma := gsl.SqrtMinFloat64
		Map1 := ra * Ma
		Mnm1 := Ma
		Mn := Map1
		for n := a + 1; n < b; n++ {
			Mnp1 := (float64(b-n)*Mnm1 + (2*float64(n)-float64(b)+x)*Mn) / float64(n)
			Mnm1 = Mn
			Mn = Mnp1
		}

		stat_ex := Exp_e(x, ex) /* 1F1(b,b,x) */
		result.val = ex.val * Ma / Mn
		result.err = ex.err * math.Abs(Ma/Mn)
		result.err += 4.0 * gsl.Float64Eps * (math.Abs(float64(b)-float64(a)) + 1.0) * math.Abs(result.val)
		return err.ErrorSelect(stat_ex, stat_CF1)
	} else if x >= 0.0 {
		if b < a {
			/* The point b,b is below the b=2a+x line.
			 * Forward recursion on a from b,b+1 is possible.
			 * Note that a > b + 1 as well, since we already tried a = b + 1.
			 */
			if x+math.Log(math.Abs(x/float64(b))) < gsl.LnMaxFloat64-2.0 {
				ex := math.Exp(x)
				Mnm1 := ex                      /* 1F1(b,b,x)   */
				Mn := ex * (1.0 + x/float64(b)) /* 1F1(b+1,b,x) */
				for n := b + 1; n < a; n++ {
					fn := float64(n)
					fb := float64(b)
					Mnp1 := ((fb-fn)*Mnm1 + (2*fn-fb+x)*Mn) / fn
					Mnm1 = Mn
					Mn = Mnp1
				}
				result.val = Mn
				result.err = (x + 1.0) * gsl.Float64Eps * math.Abs(Mn)
				result.err *= math.Abs(float64(a-b)) + 1.0
				return nil
			} else {
				return OverflowError(result)
			}
		} else {
			/* b > a
			 * b < 2a + x
			 * b <= x (otherwise we would have finished above)
			 *
			 * Gautschi anomalous convergence region. However, we can
			 * recurse forward all the way from a=0,1 because we are
			 * always underneath the b=2a+x line.
			 */
			r_Mn := new(Result)
			Mnm1 := 1.0 /* 1F1(0,b,x) */
			Exprel_n_e(b-1, x, r_Mn)
			Mn := r_Mn.val
			for n := 1; n < a; n++ {
				fn := float64(n)
				fb := float64(b)
				Mnp1 := ((fb-fn)*Mnm1 + (2*fn-fb+x)*Mn) / fn
				Mnm1 = Mn
				Mn = Mnp1
			}
			result.val = Mn
			result.err = math.Abs(Mn) * (1.0 + math.Abs(float64(a))) * math.Abs(r_Mn.err/r_Mn.val)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(Mn)
			return nil
		}
	} else {
		/* x < 0
		 * b < a (otherwise we would have tripped one of the above)
		 */

		if float64(a) <= 0.5*(float64(b)-x) || float64(a) >= -x {
			/* Gautschi continued fraction is in the anomalous region,
			 * so we must find another way. We recurse down in b,
			 * from the a=b line.
			 */
			ex := math.Exp(x)
			Manp1 := ex
			Man := ex * (1.0 + x/(float64(a)-1.0))

			for n := a - 1; n > b; n-- {
				fn := float64(n)
				fa := float64(a)
				Manm1 := (-fn*(1-fn-x)*Man - x*(fn-fa)*Manp1) / (fn * (fn - 1.0))
				Manp1 = Man
				Man = Manm1
			}
			result.val = Man
			result.err = (math.Abs(x) + 1.0) * gsl.Float64Eps * math.Abs(Man)
			result.err *= math.Abs(float64(b)-float64(a)) + 1.0
			return nil
		} else {
			/* Pick a0 such that b ~= 2a0 + x, then
			 * recurse down in b from a0,a0 to determine
			 * the values near the line b=2a+x. Then recurse
			 * forward on a from a0.
			 */
			a0 := int(math.Ceil(0.5 * (float64(b) - x)))
			var (
				Mnm1 float64
				Mn   float64
			)

			ex := math.Exp(x)
			Ma0np1 := ex
			Ma0n := ex * (1.0 + x/(float64(a0)-1.0))
			fa0 := float64(a0)
			for n := a0 - 1; n > b; n-- {
				fn := float64(n)
				Ma0nm1 := (-fn*(1-fn-x)*Ma0n - x*(fn-fa0)*Ma0np1) / (fn * (fn - 1.0))
				Ma0np1 = Ma0n
				Ma0n = Ma0nm1
			}
			Ma0bp1 := Ma0np1
			Ma0b := Ma0n
			fb := float64(b)
			Ma0p1b := (fb*(fa0+x)*Ma0b + x*(fa0-fb)*Ma0bp1) / (fa0 * fb)

			/* Initialise the recurrence correctly BJG */

			if a0 >= a {
				Mn = Ma0b
			} else if a0+1 >= a {
				Mn = Ma0p1b
			} else {
				Mnm1 = Ma0b
				Mn = Ma0p1b

				for n := a0 + 1; n < a; n++ {
					fb := float64(b)
					fn := float64(n)
					Mnp1 := ((fb-fn)*Mnm1 + (2*fn-fb+x)*Mn) / fn
					Mnm1 = Mn
					Mn = Mnp1
				}
			}

			result.val = Mn
			result.err = (math.Abs(x) + 1.0) * gsl.Float64Eps * math.Abs(Mn)
			result.err *= math.Abs((float64(b) - float64(a))) + 1.0
			return nil
		}
	}
}

/* Evaluate a <= 0, a integer, cases directly. (Polynomial; Horner)
 * When the terms are all positive, this
 * must work. We will assume this here.
 */
func hyperg_1F1_a_negint_poly(a int, b, x float64, result *Result) err.GSLError {
	if a == 0 {
		result.val = 1.0
		result.err = 1.0
		return nil
	} else {
		N := -a
		poly := 1.0

		for k := N - 1; k >= 0; k-- {
			t := float64(a+k) / (b + float64(k)) * (x / (float64(k) + 1))
			r := t + 1.0/poly
			if r > 0.9*gsl.MaxFloat64/poly {
				return OverflowError(result)
			} else {
				poly *= r /* P_n = 1 + t_n P_{n-1} */
			}
		}
		result.val = poly
		result.err = 2.0 * (math.Sqrt(float64(N)) + 1.0) * gsl.Float64Eps * math.Abs(poly)
		return nil
	}
}

/* Evaluate negative integer a case by relation
 * to Laguerre polynomials. This is more general than
 * the direct polynomial evaluation, but is safe
 * for all values of x.
 *
 * 1F1(-n,b,x) = n!/(b)_n Laguerre[n,b-1,x]
 *             = n B(b,n) Laguerre[n,b-1,x]
 *
 * assumes b is not a negative integer
 */
func hyperg_1F1_a_negint_lag(a int, b, x float64, result *Result) err.GSLError {
	n := -a

	lag := new(Result)
	stat_l := Laguerre_n_e(n, b-1.0, x, lag)
	if b < 0.0 {
		lnfact, lng1, lng2 := new(Result), new(Result), new(Result)
		var s1, s2 float64
		stat_f := Lnfact_e(uint(n), lnfact)
		stat_g1 := Lngamma_sgn_e(b+float64(n), lng1, &s1)
		stat_g2 := Lngamma_sgn_e(b, lng2, &s2)
		lnpre_val := lnfact.val - (lng1.val - lng2.val)
		lnpre_err := lnfact.err + lng1.err + lng2.err + 2.0*gsl.Float64Eps*math.Abs(lnpre_val)
		stat_e := Exp_mult_err_e(lnpre_val, lnpre_err, s1*s2*lag.val, lag.err, result)
		return err.ErrorSelect(stat_e, stat_l, stat_g1, stat_g2, stat_f)
	} else {
		lnbeta := new(Result)
		Lnbeta_e(b, float64(n), lnbeta)
		if math.Abs(lnbeta.val) < 0.1 {
			/* As we have noted, when B(x,y) is near 1,
			 * evaluating math.Log(B(x,y)) is not accurate.
			 * Instead we evaluate B(x,y) directly.
			 */
			ln_term_val := math.Log(1.25 * float64(n))
			ln_term_err := 2.0 * gsl.Float64Eps * ln_term_val
			beta := new(Result)
			stat_b := Beta_e(b, float64(n), beta)
			stat_e := Exp_mult_err_e(ln_term_val, ln_term_err, lag.val, lag.err, result)
			result.val *= beta.val / 1.25
			result.err *= beta.val / 1.25
			return err.ErrorSelect(stat_e, stat_l, stat_b)
		} else {
			/* B(x,y) was not near 1, so it is safe to use
			 * the logarithmic values.
			 */
			ln_n := math.Log(float64(n))
			ln_term_val := lnbeta.val + ln_n
			ln_term_err := lnbeta.err + 2.0*gsl.Float64Eps*math.Abs(ln_n)
			stat_e := Exp_mult_err_e(ln_term_val, ln_term_err, lag.val, lag.err, result)
			return err.ErrorSelect(stat_e, stat_l)
		}
	}
}

/* Assumes a <= -1,  b <= -1, and b <= a.
 */
func hyperg_1F1_ab_negint(a, b int, x float64, result *Result) err.GSLError {
	if x == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if x > 0.0 {
		return hyperg_1F1_a_negint_poly(a, float64(b), x, result)
	} else {
		/* Apply a Kummer transformation to make x > 0 so
		 * we can evaluate the polynomial safely. Of course,
		 * this assumes b <= a, which must be true for
		 * a<0 and b<0, since otherwise the thing is undefined.
		 */
		K := new(Result)
		stat_K := hyperg_1F1_a_negint_poly(b-a, float64(b), -x, K)
		stat_e := Exp_mult_err_e(x, 2.0*gsl.Float64Eps*math.Abs(x), K.val, K.err, result)
		return err.ErrorSelect(stat_e, stat_K)
	}
}

/* [Abramowitz+Stegun, 13.1.3]
 *
 * M(a,b,x) = Gamma(1+a-b)/Gamma(2-b) x^(1-b) *
 *            { Gamma(b)/Gamma(a) M(1+a-b,2-b,x) - (b-1) U(1+a-b,2-b,x) }
 *
 * b not an integer >= 2
 * a-b not a negative integer
 */
func hyperg_1F1_U(a, b, x float64, result *Result) err.GSLError {
	bp := 2.0 - b
	ap := a - b + 1.0

	lg_ap, lg_bp := new(Result), new(Result)
	var sg_ap float64

	stat_lg0 := Lngamma_sgn_e(ap, lg_ap, &sg_ap)
	stat_lg1 := Lngamma_e(bp, lg_bp)
	stat_lg2 := err.ErrorSelect(stat_lg0, stat_lg1)
	t1 := (bp - 1.0) * math.Log(x)
	lnpre_val := lg_ap.val - lg_bp.val + t1
	lnpre_err := lg_ap.err + lg_bp.err + 2.0*gsl.Float64Eps*math.Abs(t1)

	lg_2mbp, lg_1papmbp := new(Result), new(Result)
	var sg_2mbp, sg_1papmbp float64

	stat_lg3 := Lngamma_sgn_e(2.0-bp, lg_2mbp, &sg_2mbp)
	stat_lg4 := Lngamma_sgn_e(1.0+ap-bp, lg_1papmbp, &sg_1papmbp)
	stat_lg5 := err.ErrorSelect(stat_lg3, stat_lg4)
	lnc1_val := lg_2mbp.val - lg_1papmbp.val
	lnc1_err := lg_2mbp.err + lg_1papmbp.err + gsl.Float64Eps*(math.Abs(lg_2mbp.val)+math.Abs(lg_1papmbp.val))

	M := new(Result)
	U := new(Result_e10)
	stat_F := Hyperg_1F1_e(ap, bp, x, M)
	stat_U := Hyperg_U_e10_e(ap, bp, x, U)
	stat_FU := err.ErrorSelect(stat_F, stat_U)
	term_M := new(Result_e10)
	stat_e0 := Exp_mult_err_e10_e(lnc1_val, lnc1_err, sg_2mbp*sg_1papmbp*M.val, M.err, term_M)

	ombp := 1.0 - bp
	Uee_val := float64(U.e10) * gsl.Ln10
	Uee_err := 2.0 * gsl.Float64Eps * math.Abs(Uee_val)
	Mee_val := float64(term_M.e10) * gsl.Ln10
	Mee_err := 2.0 * gsl.Float64Eps * math.Abs(Mee_val)
	var stat_e1 err.GSLError

	/* Do a little dance with the exponential prefactors
	 * to avoid overflows in intermediate results.
	 */
	if Uee_val > Mee_val {
		factorM_val := math.Exp(Mee_val - Uee_val)
		factorM_err := 2.0 * gsl.Float64Eps * (math.Abs(Mee_val-Uee_val) + 1.0) * factorM_val
		inner_val := term_M.val*factorM_val - ombp*U.val
		inner_err := term_M.err*factorM_val + math.Abs(ombp)*U.err + math.Abs(term_M.val)*factorM_err + gsl.Float64Eps*(math.Abs(term_M.val*factorM_val)+math.Abs(ombp*U.val))
		stat_e1 = Exp_mult_err_e(lnpre_val+Uee_val, lnpre_err+Uee_err, sg_ap*inner_val, inner_err, result)
	} else {
		factorU_val := math.Exp(Uee_val - Mee_val)
		factorU_err := 2.0 * gsl.Float64Eps * (math.Abs(Mee_val-Uee_val) + 1.0) * factorU_val
		inner_val := term_M.val - ombp*factorU_val*U.val
		inner_err := term_M.err + math.Abs(ombp*factorU_val*U.err) + math.Abs(ombp*factorU_err*U.val) + gsl.Float64Eps*(math.Abs(term_M.val)+math.Abs(ombp*factorU_val*U.val))

		stat_e1 = Exp_mult_err_e(lnpre_val+Mee_val, lnpre_err+Mee_err, sg_ap*inner_val, inner_err, result)

	}

	return err.ErrorSelect(stat_e1, stat_e0, stat_FU, stat_lg5, stat_lg2)
}

/* Handle case of generic positive a, b.
 * Assumes b-a is not a negative integer.
 */
func hyperg_1F1_ab_pos(a, b, x float64, result *Result) err.GSLError {
	ax := math.Abs(x)

	if (b < 10.0 && a < 10.0 && ax < 5.0) || (b > a*ax) || (b > a && ax < 5.0) {
		return Hyperg_1F1_series_e(a, b, x, result)
	} else if x < -100.0 && math.Max(math.Abs(a), 1.0)*math.Max(math.Abs(1.0+a-b), 1.0) < 0.7*math.Abs(x) {
		/* Large negative x asymptotic.
		 */
		return hyperg_1F1_asymp_negx(a, b, x, result)
	} else if x > 100.0 && math.Max(math.Abs(b-a), 1.0)*math.Max(math.Abs(1.0-a), 1.0) < 0.7*math.Abs(x) {
		/* Large positive x asymptotic.
		 */
		return hyperg_1F1_asymp_posx(a, b, x, result)
	} else if math.Abs(b-a) <= 1.0 {
		/* Directly handle b near a.
		 */
		return hyperg_1F1_beps_bgt0(a-b, b, x, result) /* a = b + eps */
	} else if b > a && b >= 2*a+x {
		/* Use the Gautschi CF series, then
		 * recurse backward to a near 0 for normalization.
		 * This will work for either sign of x.
		 */
		var rap float64
		var n float64
		stat_CF1 := hyperg_1F1_CF1_p_ser(a, b, x, &rap)
		ra := 1.0 + x/a*rap

		Ma := gsl.SqrtMinFloat64
		Map1 := ra * Ma
		Mnp1 := Map1
		Mn := Ma
		Mn_true := new(Result)
		for n = a; n > 0.5; n -= 1.0 {
			Mnm1 := (n*Mnp1 - (2.0*n-b+x)*Mn) / (b - n)
			Mnp1 = Mn
			Mn = Mnm1
		}

		stat_Mt := hyperg_1F1_small_a_bgt0(n, b, x, Mn_true)

		result.val = (Ma / Mn) * Mn_true.val
		result.err = math.Abs(Ma/Mn) * Mn_true.err
		result.err += 2.0 * gsl.Float64Eps * (math.Abs(a) + 1.0) * math.Abs(result.val)
		return err.ErrorSelect(stat_Mt, stat_CF1)
	} else if b > a && b < 2*a+x && b > x {
		/* Use the Gautschi series representation of
		 * the continued fraction. Then recurse forward
		 * to near the a=b line for normalization. This will
		 * work for either sign of x, although we do need
		 * to check for b > x, which is relevant when x is positive.
		 */
		Mn_true := new(Result)
		var rap float64
		var n float64
		stat_CF1 := hyperg_1F1_CF1_p_ser(a, b, x, &rap)
		ra := 1.0 + x/a*rap
		Ma := gsl.SqrtMinFloat64
		Mnm1 := Ma
		Mn := ra * Mnm1

		for n = a + 1.0; n < b-0.5; n += 1.0 {
			Mnp1 := ((b-n)*Mnm1 + (2*n-b+x)*Mn) / n
			Mnm1 = Mn
			Mn = Mnp1
		}
		stat_Mt := hyperg_1F1_beps_bgt0(n-b, b, x, Mn_true)
		result.val = Ma / Mn * Mn_true.val
		result.err = math.Abs(Ma/Mn) * Mn_true.err
		result.err += 2.0 * gsl.Float64Eps * (math.Abs(b-a) + 1.0) * math.Abs(result.val)
		return err.ErrorSelect(stat_Mt, stat_CF1)
	} else if x >= 0.0 {

		if b < a {
			/* Forward recursion on a from a=b+eps-1,b+eps.
			 */
			N := math.Floor(a - b)
			eps := a - b - N
			r_M0, r_M1 := new(Result), new(Result)
			stat_0 := hyperg_1F1_beps_bgt0(eps-1.0, b, x, r_M0)
			stat_1 := hyperg_1F1_beps_bgt0(eps, b, x, r_M1)
			M0 := r_M0.val
			M1 := r_M1.val
			Mam1 := M0
			Ma := M1

			start_pair := math.Abs(M0) + math.Abs(M1)
			minim_pair := gsl.MaxFloat64

			rat_0 := math.Abs(r_M0.err / r_M0.val)
			rat_1 := math.Abs(r_M1.err / r_M1.val)
			for ap := b + eps; ap < a-0.1; ap += 1.0 {
				Map1 := ((b-ap)*Mam1 + (2.0*ap-b+x)*Ma) / ap
				Mam1 = Ma
				Ma = Map1
				minim_pair = math.Min(math.Abs(Mam1)+math.Abs(Ma), minim_pair)
			}
			pair_ratio := start_pair / minim_pair
			result.val = Ma
			result.err = 2.0 * (rat_0 + rat_1 + gsl.Float64Eps) * (math.Abs(b-a) + 1.0) * math.Abs(Ma)
			result.err += 2.0 * (rat_0 + rat_1) * pair_ratio * pair_ratio * math.Abs(Ma)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(Ma)
			return err.ErrorSelect(stat_0, stat_1)
		} else {
			/* b > a
			 * b < 2a + x
			 * b <= x
			 *
			 * Recurse forward on a from a=eps,eps+1.
			 */
			eps := a - math.Floor(a)
			r_Mnm1, r_Mn := new(Result), new(Result)
			stat_0 := hyperg_1F1_small_a_bgt0(eps, b, x, r_Mnm1)
			stat_1 := hyperg_1F1_small_a_bgt0(eps+1.0, b, x, r_Mn)
			Mnm1 := r_Mnm1.val
			Mn := r_Mn.val

			start_pair := math.Abs(Mn) + math.Abs(Mnm1)
			minim_pair := gsl.MaxFloat64
			rat_0 := math.Abs(r_Mnm1.err / r_Mnm1.val)
			rat_1 := math.Abs(r_Mn.err / r_Mn.val)
			for n := eps + 1.0; n < a-0.1; n++ {
				Mnp1 := ((b-n)*Mnm1 + (2*n-b+x)*Mn) / n
				Mnm1 = Mn
				Mn = Mnp1
				minim_pair = math.Min(math.Abs(Mn)+math.Abs(Mnm1), minim_pair)
			}
			pair_ratio := start_pair / minim_pair
			result.val = Mn
			result.err = 2.0 * (rat_0 + rat_1 + gsl.Float64Eps) * (math.Abs(a) + 1.0) * math.Abs(Mn)
			result.err += 2.0 * (rat_0 + rat_1) * pair_ratio * pair_ratio * math.Abs(Mn)
			result.err += 2.0 * gsl.Float64Eps * math.Abs(Mn)
			return err.ErrorSelect(stat_0, stat_1)
		}
	} else {
		/* x < 0
		 * b < a
		 */

		if a <= 0.5*(b-x) || a >= -x {
			/* Recurse down in b, from near the a=b line, b=a+eps,a+eps-1.
			 */
			N := math.Floor(a - b)
			eps := 1.0 + N - a + b
			r_Manp1, r_Man := new(Result), new(Result)
			stat_0 := hyperg_1F1_beps_bgt0(-eps, a+eps, x, r_Manp1)
			stat_1 := hyperg_1F1_beps_bgt0(1.0-eps, a+eps-1.0, x, r_Man)
			Manp1 := r_Manp1.val
			Man := r_Man.val
			start_pair := math.Abs(Manp1) + math.Abs(Man)
			minim_pair := math.MaxFloat64
			rat_0 := math.Abs(r_Manp1.err / r_Manp1.val)
			rat_1 := math.Abs(r_Man.err / r_Man.val)
			for n := a + eps - 1.0; n > b+0.1; n -= 1.0 {
				Manm1 := (-n*(1-n-x)*Man - x*(n-a)*Manp1) / (n * (n - 1.0))
				Manp1 = Man
				Man = Manm1
				minim_pair = math.Min(math.Abs(Manp1)+math.Abs(Man), minim_pair)
			}

			/* FIXME: this is a nasty little hack; there is some
			   (transient?) instability in this recurrence for some
			   values. I can tell when it happens, which is when
			   this pair_ratio is large. But I do not know how to
			   measure the error in terms of it. I guessed quadratic
			   below, but it is probably worse than that.
			*/
			pair_ratio := start_pair / minim_pair
			result.val = Man
			result.err = 2.0 * (rat_0 + rat_1 + gsl.Float64Eps) * (math.Abs(b-a) + 1.0) * math.Abs(Man)
			result.err *= pair_ratio*pair_ratio + 1.0
			return err.ErrorSelect(stat_0, stat_1)
		} else {
			/* Pick a0 such that b ~= 2a0 + x, then
			 * recurse down in b from a0,a0 to determine
			 * the values near the line b=2a+x. Then recurse
			 * forward on a from a0.
			 */
			var (
				epsa    float64 = a - math.Floor(a)
				a0      float64 = math.Floor(0.5*(b-x)) + epsa
				N       float64 = math.Floor(a0 - b)
				epsb    float64 = 1.0 + N - a0 + b
				Ma0b    float64
				Ma0bp1  float64
				Ma0p1b  float64
				Mnm1    float64
				Mn      float64
				Mnp1    float64
				n       float64
				err_rat float64
			)

			r_Ma0np1, r_Ma0n := new(Result), new(Result)
			stat_0 := hyperg_1F1_beps_bgt0(-epsb, a0+epsb, x, r_Ma0np1)
			stat_1 := hyperg_1F1_beps_bgt0(1.0-epsb, a0+epsb-1.0, x, r_Ma0n)
			Ma0np1 := r_Ma0np1.val
			Ma0n := r_Ma0n.val
			err_rat = math.Abs(r_Ma0np1.err/r_Ma0np1.val) + math.Abs(r_Ma0n.err/r_Ma0n.val)

			for n = a0 + epsb - 1.0; n > b+0.1; n -= 1.0 {
				Ma0nm1 := (-n*(1-n-x)*Ma0n - x*(n-a0)*Ma0np1) / (n * (n - 1.0))
				Ma0np1 = Ma0n
				Ma0n = Ma0nm1
			}
			Ma0bp1 = Ma0np1
			Ma0b = Ma0n
			Ma0p1b = (b*(a0+x)*Ma0b + x*(a0-b)*Ma0bp1) / (a0 * b) /* right-down hook */
			stat_a0 := err.ErrorSelect(stat_0, stat_1)

			/* Initialise the recurrence correctly BJG */

			if a0 >= a-0.1 {
				Mn = Ma0b
			} else if a0+1 >= a-0.1 {
				Mn = Ma0p1b
			} else {
				Mnm1 = Ma0b
				Mn = Ma0p1b

				for n = a0 + 1.0; n < a-0.1; n += 1.0 {
					Mnp1 = ((b-n)*Mnm1 + (2*n-b+x)*Mn) / n
					Mnm1 = Mn
					Mn = Mnp1
				}
			}

			result.val = Mn
			result.err = (err_rat + gsl.Float64Eps) * (math.Abs(b-a) + 1.0) * math.Abs(Mn)
			return stat_a0
		}
	}
}

/* Assumes b != integer
 * Assumes a != integer when x > 0
 * Assumes b-a != neg integer when x < 0
 */
func hyperg_1F1_ab_neg(a, b, x float64, result *Result) err.GSLError {
	bma := b - a
	abs_x := math.Abs(x)
	abs_a := math.Abs(a)
	abs_b := math.Abs(b)
	size_a := math.Max(abs_a, 1.0)
	size_b := math.Max(abs_b, 1.0)
	bma_integer := (bma-math.Floor(bma+0.5) < _1F1_INT_THRESHOLD)

	if (abs_a < 10.0 && abs_b < 10.0 && abs_x < 5.0) || (b > 0.8*math.Max(math.Abs(a), 1.0)*math.Abs(x)) {
		return Hyperg_1F1_series_e(a, b, x, result)
	} else if x > 0.0 && size_b > size_a && size_a*math.Log(math.E*x/size_b) < gsl.LnFloat64Eps+7.0 {
		/* Series terms are positive definite up until
		 * there is a sign change. But by then the
		 * terms are small due to the last condition.
		 */
		return Hyperg_1F1_series_e(a, b, x, result)
	} else if (abs_x < 5.0 && math.Abs(bma) < 10.0 && abs_b < 10.0) || (b > 0.8*math.Max(math.Abs(bma), 1.0)*abs_x) {
		/* Use Kummer transformation to render series safe.
		 */
		Kummer_1F1 := new(Result)
		stat_K := Hyperg_1F1_series_e(bma, b, -x, Kummer_1F1)
		stat_e := Exp_mult_err_e(x, gsl.Float64Eps*math.Abs(x), Kummer_1F1.val, Kummer_1F1.err, result)
		return err.ErrorSelect(stat_e, stat_K)
	} else if x < -30.0 && math.Max(math.Abs(a), 1.0)*math.Max(math.Abs(1.0+a-b), 1.0) < 0.99*math.Abs(x) {
		/* Large negative x asymptotic.
		 * Note that we do not check if b-a is a negative integer.
		 */
		return hyperg_1F1_asymp_negx(a, b, x, result)
	} else if x > 100.0 && math.Max(math.Abs(bma), 1.0)*math.Max(math.Abs(1.0-a), 1.0) < 0.99*math.Abs(x) {
		/* Large positive x asymptotic.
		 * Note that we do not check if a is a negative integer.
		 */
		return hyperg_1F1_asymp_posx(a, b, x, result)
	} else if x > 0.0 && !(bma_integer && bma > 0.0) {
		return hyperg_1F1_U(a, b, x, result)
	} else {
		/* FIXME:  if all else fails, try the series... BJG */
		if x < 0.0 {
			/* Apply Kummer Transformation */
			status := Hyperg_1F1_series_e(b-a, b, -x, result)
			K_factor := math.Exp(x)
			result.val *= K_factor
			result.err *= K_factor
			return status
		} else {
			return Hyperg_1F1_series_e(a, b, x, result)
		}

		/* Sadness... */
		/* result.val = 0.0; */
		/* result.err = 0.0; */
		/* GSL_ERROR ("error", GSL_EUNIMPL); */
	}
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Hyperg_1F1_int_e(a, b int, x float64, result *Result) err.GSLError {
	if x == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if a == b {
		return Exp_e(x, result)
	} else if b == 0 {
		return DomainError(result)
	} else if a == 0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if b < 0 && (a < b || a > 0) {
		/* Standard domain error due to singularity. */
		return DomainError(result)
	} else if x > 100.0 && math.Max(1.0, math.Abs(float64(b-a)))*math.Max(1.0, math.Abs(float64(1-a))) < 0.5*x {
		/* x . +Inf asymptotic */
		return hyperg_1F1_asymp_posx(float64(a), float64(b), x, result)
	} else if x < -100.0 && math.Max(1.0, math.Abs(float64(a)))*math.Max(1.0, math.Abs(float64(1+a-b))) < 0.5*math.Abs(x) {
		/* x . -Inf asymptotic */
		return hyperg_1F1_asymp_negx(float64(a), float64(b), x, result)
	} else if a < 0 && b < 0 {
		return hyperg_1F1_ab_negint(a, b, x, result)
	} else if a < 0 && b > 0 {
		/* Use Kummer to reduce it to the positive integer case.
		 * Note that b > a, strictly, since we already trapped b = a.
		 */
		Kummer_1F1 := new(Result)
		stat_K := hyperg_1F1_ab_posint(b-a, b, -x, Kummer_1F1)
		stat_e := Exp_mult_err_e(x, gsl.Float64Eps*math.Abs(x), Kummer_1F1.val, Kummer_1F1.err, result)
		return err.ErrorSelect(stat_e, stat_K)
	} else {
		/* a > 0 and b > 0 */
		return hyperg_1F1_ab_posint(a, b, x, result)
	}
}

func Hyperg_1F1_e(a, b, x float64, result *Result) err.GSLError {
	bma := b - a
	rinta := math.Floor(a + 0.5)
	rintb := math.Floor(b + 0.5)
	rintbma := math.Floor(bma + 0.5)
	a_integer := (math.Abs(a-rinta) < _1F1_INT_THRESHOLD && rinta > math.MinInt32 && rinta < math.MaxInt32)
	b_integer := (math.Abs(b-rintb) < _1F1_INT_THRESHOLD && rintb > math.MinInt32 && rintb < math.MaxInt32)
	bma_integer := (math.Abs(bma-rintbma) < _1F1_INT_THRESHOLD && rintbma > math.MinInt32 && rintbma < math.MaxInt32)
	b_neg_integer := (b < -0.1 && b_integer)
	a_neg_integer := (a < -0.1 && a_integer)
	bma_neg_integer := (bma < -0.1 && bma_integer)

	if x == 0.0 {
		/* Testing for this before testing a and b
		 * is somewhat arbitrary. The result is that
		 * we have 1F1(a,0,0) = 1.
		 */
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if b == 0.0 {
		return DomainError(result)
	} else if a == 0.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if a == b {
		/* case: a==b; exp(x)
		 * It's good to test exact equality now.
		 * We also test approximate equality later.
		 */
		return Exp_e(x, result)
	} else if math.Abs(b) < _1F1_INT_THRESHOLD && math.Abs(a) < _1F1_INT_THRESHOLD {
		/* a and b near zero: 1 + a/b (exp(x)-1)
		 */

		/* Note that neither a nor b is zero, since
		 * we eliminated that with the above tests.
		 */

		exm1 := new(Result)
		stat_e := Expm1_e(x, exm1)

		sa := -1.
		if a > 0.0 {
			sa = 1.
		}

		sb := -1.
		if b > 0.0 {
			sb = 1.
		}
		lnab := math.Log(math.Abs(a / b)) /* safe */
		hx := new(Result)
		stat_hx := Exp_mult_err_e(lnab, gsl.Float64Eps*math.Abs(lnab), sa*sb*exm1.val, exm1.err, hx)
		result.val = hx.val /* FIXME: excessive paranoia ? what is DBL_MAX+1 ?*/
		if hx.val != gsl.MaxFloat64 {
			result.val++
		}
		result.err = hx.err
		return err.ErrorSelect(stat_hx, stat_e)
	} else if math.Abs(b) < _1F1_INT_THRESHOLD && math.Abs(x*a) < 1 {
		/* b near zero and a not near zero
		 */

		m_arg := 1.0 / (0.5 * b)
		F_renorm := new(Result)
		stat_F := hyperg_1F1_renorm_b0(a, x, F_renorm)
		stat_m := Multiply_err_e(m_arg, 2.0*gsl.Float64Eps*m_arg, 0.5*F_renorm.val, 0.5*F_renorm.err, result)
		return err.ErrorSelect(stat_m, stat_F)
	} else if a_integer && b_integer {
		/* Check for reduction to the integer case.
		 * Relies on the arbitrary "near an integer" test.
		 */
		return Hyperg_1F1_int_e(int(rinta), int(rintb), x, result)
	} else if b_neg_integer && !(a_neg_integer && a > b) {
		/* Standard domain error due to
		 * uncancelled singularity.
		 */

		return DomainError(result)
	} else if a_neg_integer {

		return hyperg_1F1_a_negint_lag(int(rinta), b, x, result)
	} else if b > 0.0 {
		if -1.0 <= a && a <= 1.0 {
			/* Handle small a explicitly.
			 */

			return hyperg_1F1_small_a_bgt0(a, b, x, result)
		} else if bma_neg_integer {
			/* Catch this now, to avoid problems in the
			 * generic evaluation code.
			 */

			Kummer_1F1 := new(Result)
			stat_K := hyperg_1F1_a_negint_lag(int(rintbma), b, -x, Kummer_1F1)
			stat_e := Exp_mult_err_e(x, gsl.Float64Eps*math.Abs(x), Kummer_1F1.val, Kummer_1F1.err, result)
			return err.ErrorSelect(stat_e, stat_K)
		} else if a < 0.0 && math.Abs(x) < 2*gsl.LnMaxFloat64 {
			/* Use Kummer to reduce it to the generic positive case.
			 * Note that b > a, strictly, since we already trapped b = a.
			 * Also b-(b-a)=a, and a is not a negative integer here,
			 * so the generic evaluation is safe.
			 */

			Kummer_1F1 := new(Result)
			stat_K := hyperg_1F1_ab_pos(b-a, b, -x, Kummer_1F1)
			stat_e := Exp_mult_err_e(x, gsl.Float64Eps*math.Abs(x), Kummer_1F1.val, Kummer_1F1.err, result)
			return err.ErrorSelect(stat_e, stat_K)
		} else if a > 0 {
			/* a > 0.0 */

			return hyperg_1F1_ab_pos(a, b, x, result)
		} else {

			return Hyperg_1F1_series_e(a, b, x, result)
		}
	} else {
		/* b < 0.0 */

		if bma_neg_integer && x < 0.0 {
			/* Handle this now to prevent problems
			 * in the generic evaluation.
			 */
			K := new(Result)
			var stat_K, stat_e err.GSLError
			if a < 0.0 {
				/* Kummer transformed version of safe polynomial.
				 * The condition a < 0 is equivalent to b < b-a,
				 * which is the condition required for the series
				 * to be positive definite here.
				 */
				stat_K = hyperg_1F1_a_negint_poly(int(rintbma), b, -x, K)
			} else {
				/* Generic eval for negative integer a. */
				stat_K = hyperg_1F1_a_negint_lag(int(rintbma), b, -x, K)
			}
			stat_e = Exp_mult_err_e(x, gsl.Float64Eps*math.Abs(x), K.val, K.err, result)
			return err.ErrorSelect(stat_e, stat_K)
		} else if a > 0.0 {
			/* Use Kummer to reduce it to the generic negative case.
			 */

			K := new(Result)
			stat_K := hyperg_1F1_ab_neg(b-a, b, -x, K)
			stat_e := Exp_mult_err_e(x, gsl.Float64Eps*math.Abs(x), K.val, K.err, result)
			return err.ErrorSelect(stat_e, stat_K)
		} else {
			s := hyperg_1F1_ab_neg(a, b, x, result)
			return s
		}
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Hyperg_1F1_int(m, n int, x float64) float64 {
	result := new(Result)
	status := Hyperg_1F1_int_e(m, n, x, result)
	return EvalResult(result, status)
}

func Hyperg_1F1(a, b, x float64) float64 {
	result := new(Result)
	status := Hyperg_1F1_e(a, b, x, result)
	return EvalResult(result, status)
}
