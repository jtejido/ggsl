/* specfunc/coupling.c
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

func locMax3(a, b, c int) int {
	d := gsl.Max(a, b)
	return int(gsl.Max(d, c))
}

func locMin3(a, b, c int) int {
	d := gsl.Min(a, b)
	return int(gsl.Min(d, c))
}

func locMin5(a, b, c, d, e int) int {
	f := int(gsl.Min(a, b))
	g := int(gsl.Min(c, d))
	h := int(gsl.Min(f, g))
	return int(gsl.Min(e, h))
}

/* See: [Thompson, Atlas for Computing Mathematical Functions] */
func delta(ta, tb, tc int, d *Result) err.GSLError {
	var f1, f2, f3, f4 Result
	if stat := Fact_e(uint(ta+tb-tc)/2, &f1); stat != nil {
		return OverflowError(d)
	}
	if stat := Fact_e(uint(ta+tc-tb)/2, &f2); stat != nil {
		return OverflowError(d)
	}
	if stat := Fact_e(uint(tb+tc-ta)/2, &f3); stat != nil {
		return OverflowError(d)
	}
	if stat := Fact_e(uint(ta+tb+tc)/2+1, &f4); stat != nil {
		return OverflowError(d)
	}

	d.val = f1.val * f2.val * f3.val / f4.val
	d.err = 4.0 * gsl.Float64Eps * math.Abs(d.val)
	return nil
}

func triangle_selection_fails(two_ja, two_jb, two_jc int) bool {
	/*
	 * enough to check the triangle condition for one spin vs. the other two
	 */
	return ((two_jb < gsl.AbsInt(two_ja-two_jc)) || (two_jb > two_ja+two_jc) || gsl.IsOdd(two_ja+two_jb+two_jc))
}

func m_selection_fails(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc int) bool {
	return (gsl.AbsInt(two_ma) > two_ja || gsl.AbsInt(two_mb) > two_jb || gsl.AbsInt(two_mc) > two_jc || gsl.IsOdd(two_ja+two_ma) || gsl.IsOdd(two_jb+two_mb) || gsl.IsOdd(two_jc+two_mc) || (two_ma+two_mb+two_mc) != 0)
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Coupling_3j_e(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc int, result *Result) err.GSLError {

	if two_ja < 0 || two_jb < 0 || two_jc < 0 {
		return DomainError(result)
	} else if triangle_selection_fails(two_ja, two_jb, two_jc) || m_selection_fails(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc) {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if two_ma == 0 && two_mb == 0 && two_mc == 0 && ((two_ja+two_jb+two_jc)%4 == 2) {
		/* Special case for (ja jb jc; 0 0 0) = 0 when ja+jb+jc=odd */
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		var (
			jca  = (-two_ja + two_jb + two_jc) / 2
			jcb  = (two_ja - two_jb + two_jc) / 2
			jcc  = (two_ja + two_jb - two_jc) / 2
			jmma = (two_ja - two_ma) / 2
			jmmb = (two_jb - two_mb) / 2
			jmmc = (two_jc - two_mc) / 2
			jpma = (two_ja + two_ma) / 2
			jpmb = (two_jb + two_mb) / 2
			jpmc = (two_jc + two_mc) / 2
			jsum = (two_ja + two_jb + two_jc) / 2
			kmin = locMax3(0, jpmb-jmmc, jmma-jpmc)
			kmax = locMin3(jcc, jmma, jpmb)
			k    int
			sign = 1
		)

		if gsl.IsOdd(kmin - jpma + jmmb) {
			sign = -1
		}

		var sum_pos, sum_neg, sum_err float64
		var bc1, bc2, bc3, bcn1, bcn2, bcd1, bcd2, bcd3, bcd4, term, lnorm Result
		var s int

		if status := Lnchoose_e(uint(two_ja), uint(jcc), &bcn1); status != nil {
			s += status.Status()
		}
		if status := Lnchoose_e(uint(two_jb), uint(jcc), &bcn2); status != nil {
			s += status.Status()
		}
		if status := Lnchoose_e(uint(jsum+1), uint(jcc), &bcd1); status != nil {
			s += status.Status()
		}
		if status := Lnchoose_e(uint(two_ja), uint(jmma), &bcd2); status != nil {
			s += status.Status()
		}
		if status := Lnchoose_e(uint(two_jb), uint(jmmb), &bcd3); status != nil {
			s += status.Status()
		}
		if status := Lnchoose_e(uint(two_jc), uint(jpmc), &bcd4); status != nil {

			s += status.Status()
		}

		lnorm.val = 0.5 * (bcn1.val + bcn2.val - bcd1.val - bcd2.val - bcd3.val - bcd4.val - math.Log(float64(two_jc)+1.0))
		lnorm.err = 0.5 * (bcn1.err + bcn2.err + bcd1.err + bcd2.err + bcd3.err + bcd4.err + gsl.Float64Eps*math.Log(float64(two_jc)+1.0))

		for k = kmin; k <= kmax; k++ {
			if status := Lnchoose_e(uint(jcc), uint(k), &bc1); status != nil {
				s += status.Status()
			}
			if status := Lnchoose_e(uint(jcb), uint(jmma-k), &bc2); status != nil {
				s += status.Status()
			}
			if status := Lnchoose_e(uint(jca), uint(jpmb-k), &bc3); status != nil {
				s += status.Status()
			}
			if status := Exp_err_e(bc1.val+bc2.val+bc3.val+lnorm.val, bc1.err+bc2.err+bc3.err+lnorm.err, &term); status != nil {
				s += status.Status()
			}

			if s != 0 {
				return OverflowError(result)
			}

			if sign < 0 {
				sum_neg += term.val
			} else {
				sum_pos += term.val
			}

			sum_err += term.err

			sign = -sign
		}

		result.val = sum_pos - sum_neg
		result.err = sum_err
		result.err += 2.0 * gsl.Float64Eps * (sum_pos + sum_neg)
		result.err += 2.0 * gsl.Float64Eps * (float64(kmax) - float64(kmin)) * math.Abs(result.val)

		return nil
	}
}

func Coupling_6j_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf int, result *Result) err.GSLError {

	if two_ja < 0 || two_jb < 0 || two_jc < 0 || two_jd < 0 || two_je < 0 || two_jf < 0 {
		return DomainError(result)
	} else if triangle_selection_fails(two_ja, two_jb, two_jc) || triangle_selection_fails(two_ja, two_je, two_jf) || triangle_selection_fails(two_jb, two_jd, two_jf) || triangle_selection_fails(two_je, two_jd, two_jc) {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		var n1, d1, d2, d3, d4, d5, d6 Result
		var norm float64
		var tk, tkmin, tkmax int
		var phase, sum_pos, sum_neg, sumsq_err float64
		var s int

		if status := delta(two_ja, two_jb, two_jc, &d1); status != nil {
			s++
		}
		if status := delta(two_ja, two_je, two_jf, &d2); status != nil {
			s++
		}
		if status := delta(two_jb, two_jd, two_jf, &d3); status != nil {
			s++
		}
		if status := delta(two_je, two_jd, two_jc, &d4); status != nil {
			s++
		}

		if s != 0 {
			return OverflowError(result)
		}

		norm = math.Sqrt(d1.val) * math.Sqrt(d2.val) * math.Sqrt(d3.val) * math.Sqrt(d4.val)

		tkmin = locMax3(0,
			two_ja+two_jd-two_jc-two_jf,
			two_jb+two_je-two_jc-two_jf)

		tkmax = locMin5(two_ja+two_jb+two_je+two_jd+2,
			two_ja+two_jb-two_jc,
			two_je+two_jd-two_jc,
			two_ja+two_je-two_jf,
			two_jb+two_jd-two_jf)

		phase = 1.0
		if gsl.IsOdd((two_ja + two_jb + two_je + two_jd + tkmin) / 2) {
			phase = -1.0
		}

		for tk = tkmin; tk <= tkmax; tk += 2 {
			var term, term_err float64
			var den_1, den_2, d1_a, d1_b Result
			var ss int

			if status := Fact_e(uint((two_ja+two_jb+two_je+two_jd-tk)/2+1), &n1); status != nil {
				ss++
			}
			if status := Fact_e(uint(tk/2), &d1_a); status != nil {
				ss++
			}
			if status := Fact_e(uint((two_jc+two_jf-two_ja-two_jd+tk)/2), &d1_b); status != nil {
				ss++
			}
			if status := Fact_e(uint((two_jc+two_jf-two_jb-two_je+tk)/2), &d2); status != nil {
				ss++
			}
			if status := Fact_e(uint((two_ja+two_jb-two_jc-tk)/2), &d3); status != nil {
				ss++
			}
			if status := Fact_e(uint((two_je+two_jd-two_jc-tk)/2), &d4); status != nil {
				ss++
			}
			if status := Fact_e(uint((two_ja+two_je-two_jf-tk)/2), &d5); status != nil {
				ss++
			}
			if status := Fact_e(uint((two_jb+two_jd-two_jf-tk)/2), &d6); status != nil {
				ss++
			}

			if ss != 0 {
				return OverflowError(result)
			}

			d1.val = d1_a.val * d1_b.val
			d1.err = d1_a.err*math.Abs(d1_b.val) + math.Abs(d1_a.val)*d1_b.err

			den_1.val = d1.val * d2.val * d3.val
			den_1.err = d1.err * math.Abs(d2.val*d3.val)
			den_1.err += d2.err * math.Abs(d1.val*d3.val)
			den_1.err += d3.err * math.Abs(d1.val*d2.val)

			den_2.val = d4.val * d5.val * d6.val
			den_2.err = d4.err * math.Abs(d5.val*d6.val)
			den_2.err += d5.err * math.Abs(d4.val*d6.val)
			den_2.err += d6.err * math.Abs(d4.val*d5.val)

			term = phase * n1.val / den_1.val / den_2.val
			phase = -phase
			term_err = n1.err / math.Abs(den_1.val) / math.Abs(den_2.val)
			term_err += math.Abs(term/den_1.val) * den_1.err
			term_err += math.Abs(term/den_2.val) * den_2.err

			if term >= 0.0 {
				sum_pos += norm * term
			} else {
				sum_neg -= norm * term
			}

			sumsq_err += norm * norm * term_err * term_err
		}

		result.val = sum_pos - sum_neg
		result.err = 2.0 * gsl.Float64Eps * (sum_pos + sum_neg)
		result.err += math.Sqrt(sumsq_err / (0.5*float64(tkmax-tkmin) + 1.0))
		result.err += 2.0 * gsl.Float64Eps * float64(tkmax-tkmin+2) * math.Abs(result.val)

		return nil
	}
}

func Coupling_RacahW_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf int, result *Result) err.GSLError {
	status := Coupling_6j_e(two_ja, two_jb, two_je, two_jd, two_jc, two_jf, result)
	phase_sum := (two_ja + two_jb + two_jc + two_jd) / 2

	if gsl.IsOdd(phase_sum) {
		result.val *= -1
	}
	return status
}

func Coupling_9j_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji int, result *Result) err.GSLError {

	if two_ja < 0 || two_jb < 0 || two_jc < 0 || two_jd < 0 || two_je < 0 || two_jf < 0 || two_jg < 0 || two_jh < 0 || two_ji < 0 {
		return DomainError(result)
	} else if triangle_selection_fails(two_ja, two_jb, two_jc) || triangle_selection_fails(two_jd, two_je, two_jf) || triangle_selection_fails(two_jg, two_jh, two_ji) || triangle_selection_fails(two_ja, two_jd, two_jg) || triangle_selection_fails(two_jb, two_je, two_jh) || triangle_selection_fails(two_jc, two_jf, two_ji) {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		var tk int
		tkmin := locMax3(gsl.AbsInt(two_ja-two_ji), gsl.AbsInt(two_jh-two_jd), gsl.AbsInt(two_jb-two_jf))
		tkmax := locMin3(two_ja+two_ji, two_jh+two_jd, two_jb+two_jf)
		sum_pos := 0.0
		sum_neg := 0.0
		sumsq_err := 0.0
		phase := 1.0
		for tk = tkmin; tk <= tkmax; tk += 2 {
			var s1, s2, s3 Result
			var term, term_err float64
			s := 0

			if status := Coupling_6j_e(two_ja, two_ji, tk, two_jh, two_jd, two_jg, &s1); status != nil {
				s++
			}
			if status := Coupling_6j_e(two_jb, two_jf, tk, two_jd, two_jh, two_je, &s2); status != nil {
				s++
			}
			if status := Coupling_6j_e(two_ja, two_ji, tk, two_jf, two_jb, two_jc, &s3); status != nil {
				s++
			}

			if s != 0 {
				return OverflowError(result)
			}
			term = s1.val * s2.val * s3.val
			term_err = s1.err * math.Abs(s2.val*s3.val)
			term_err += s2.err * math.Abs(s1.val*s3.val)
			term_err += s3.err * math.Abs(s1.val*s2.val)

			if term >= 0.0 {
				sum_pos += float64(tk+1) * term
			} else {
				sum_neg -= float64(tk+1) * term
			}

			sumsq_err += (float64(tk+1) * term_err) * (float64(tk+1) * term_err)
		}

		if gsl.IsOdd(tkmin) {
			phase = -1.
		}

		result.val = phase * (sum_pos - sum_neg)
		result.err = 2.0 * gsl.Float64Eps * (sum_pos + sum_neg)
		result.err += math.Sqrt(sumsq_err / (0.5*float64(tkmax-tkmin) + 1.0))
		result.err += 2.0 * gsl.Float64Eps * float64(tkmax-tkmin+2) * math.Abs(result.val)

		return nil
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc int) float64 {
	result := new(Result)
	status := Coupling_3j_e(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc, result)
	return EvalResult(result, status)
}

func Coupling_6j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf int) float64 {
	result := new(Result)
	status := Coupling_6j_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, result)
	return EvalResult(result, status)
}

func Coupling_RacahW(two_ja, two_jb, two_jc, two_jd, two_je, two_jf int) float64 {
	result := new(Result)
	status := Coupling_RacahW_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, result)
	return EvalResult(result, status)
}

func Coupling_9j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji int) float64 {
	result := new(Result)
	status := Coupling_9j_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji, result)
	return EvalResult(result, status)

}
