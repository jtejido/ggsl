/* specfunc/legendre_poly.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002 Gerard Jungman
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

/* Calculate P_m^m(x) from the analytic result:
 *   P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2) , m > 0
 *            = 1 , m = 0
 */
func legendre_Pmm(m int, x float64) float64 {
	if m == 0 {
		return 1.0
	} else {
		p_mm := 1.0
		root_factor := math.Sqrt(1.0-x) * math.Sqrt(1.0+x)
		fact_coeff := 1.0
		var i int
		for i = 1; i <= m; i++ {
			p_mm *= -fact_coeff * root_factor
			fact_coeff += 2.0
		}
		return p_mm
	}
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Legendre_P1_e(x float64, result *Result) err.GSLError {
	result.val = x
	result.err = 0.0
	return nil
}

func Legendre_P2_e(x float64, result *Result) err.GSLError {
	result.val = 0.5 * (3.0*x*x - 1.0)
	result.err = gsl.Float64Eps * (math.Abs(3.0*x*x) + 1.0)
	return nil
}

func Legendre_P3_e(x float64, result *Result) err.GSLError {
	result.val = 0.5 * x * (5.0*x*x - 3.0)
	result.err = gsl.Float64Eps * (math.Abs(result.val) + 0.5*math.Abs(x)*(math.Abs(5.0*x*x)+3.0))
	return nil

}

func Legendre_Pl_e(l int, x float64, result *Result) err.GSLError {

	if l < 0 || x < -1.0 || x > 1.0 {
		return DomainError(result)
	} else if l == 0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if l == 1 {
		result.val = x
		result.err = 0.0
		return nil
	} else if l == 2 {
		result.val = 0.5 * (3.0*x*x - 1.0)
		result.err = gsl.Float64Eps * (math.Abs(3.0*x*x) + 1.0)
		/*result.err = 3.0 * gsl.Float64Eps * math.Abs(result.val);
		  removed this old bogus estimate [GJ]
		*/
		return nil
	} else if x == 1.0 {
		result.val = 1.0
		result.err = 0.0
		return nil
	} else if x == -1.0 {
		result.val = 1.0
		if gsl.IsOdd(l) {
			result.val = -1.0
		}
		result.err = 0.0
		return nil
	} else if l < 100000 {
		/* upward recurrence: l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2} */
		var (
			p_ellm2 = 1.0 /* P_0(x) */
			p_ellm1 = x   /* P_1(x) */
			p_ell   = p_ellm1
			e_ellm2 = gsl.Float64Eps
			e_ellm1 = math.Abs(x) * gsl.Float64Eps
			e_ell   = e_ellm1
			ell     int
		)

		for ell = 2; ell <= l; ell++ {
			fell := float64(ell)
			p_ell = (x*(2*fell-1)*p_ellm1 - (fell-1)*p_ellm2) / fell
			p_ellm2 = p_ellm1
			p_ellm1 = p_ell

			e_ell = 0.5 * (math.Abs(x)*(2*fell-1.0)*e_ellm1 + (fell-1.0)*e_ellm2) / fell
			e_ellm2 = e_ellm1
			e_ellm1 = e_ell
		}

		result.val = p_ell
		result.err = e_ell + float64(l)*math.Abs(p_ell)*gsl.Float64Eps
		return nil
	} else {
		/* Asymptotic expansion.
		 * [Olver, p. 473]
		 */
		var (
			u            = float64(l) + 0.5
			th           = math.Acos(x)
			J0, Jm1      Result
			stat_J0      = Bessel_J0_e(u*th, &J0)
			stat_Jm1     = Bessel_Jn_e(-1, u*th, &Jm1)
			pre, B00, c1 float64
		)
		/* B00 = 1/8 (1 - th cot(th) / th^2
		 * pre = math.Sqrt(th/sin(th))
		 */
		if th < gsl.Root4Float64Eps {
			B00 = (1.0 + th*th/15.0) / 24.0
			pre = 1.0 + th*th/12.0
		} else {
			sin_th := math.Sqrt(1.0 - x*x)
			cot_th := x / sin_th
			B00 = 1.0 / 8.0 * (1.0 - th*cot_th) / (th * th)
			pre = math.Sqrt(th / sin_th)
		}

		c1 = th / u * B00

		result.val = pre * (J0.val + c1*Jm1.val)
		result.err = pre * (J0.err + math.Abs(c1)*Jm1.err)
		result.err += gsl.SqrtFloat64Eps * math.Abs(result.val)

		return err.ErrorSelect(stat_J0, stat_Jm1)
	}
}

func Legendre_Pl_array(lmax int, x float64, result_array []float64) err.GSLError {

	if lmax < 0 || x < -1.0 || x > 1.0 {
		return err.ERROR("domain error", err.EDOM)
	} else if lmax == 0 {
		result_array[0] = 1.0
		return nil
	} else if lmax == 1 {
		result_array[0] = 1.0
		result_array[1] = x
		return nil
	} else {
		/* upward recurrence: l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2} */

		p_ellm2 := 1.0 /* P_0(x) */
		p_ellm1 := x   /* P_1(x) */
		p_ell := p_ellm1
		var ell int

		result_array[0] = 1.0
		result_array[1] = x

		for ell = 2; ell <= lmax; ell++ {
			fell := float64(ell)
			p_ell = (x*(2*fell-1)*p_ellm1 - (fell-1)*p_ellm2) / fell
			p_ellm2 = p_ellm1
			p_ellm1 = p_ell
			result_array[ell] = p_ell
		}

		return nil
	}
}

func Legendre_Pl_deriv_array(lmax int, x float64, result_array, result_deriv_array []float64) err.GSLError {
	stat_array := Legendre_Pl_array(lmax, x, result_array)
	flmax := float64(lmax)
	if lmax >= 0 {
		result_deriv_array[0] = 0.0
	}
	if lmax >= 1 {
		result_deriv_array[1] = 1.0
	}

	if stat_array == nil {
		var ell int

		if math.Abs(x-1.0)*(flmax+1.0)*(flmax+1.0) < gsl.SqrtFloat64Eps {
			/* x is near 1 */
			for ell = 2; ell <= lmax; ell++ {
				fell := float64(ell)
				pre := 0.5 * fell * (fell + 1.0)
				result_deriv_array[ell] = pre * (1.0 - 0.25*(1.0-x)*(fell+2.0)*(fell-1.0))
			}
		} else if math.Abs(x+1.0)*(flmax+1.0)*(flmax+1.0) < gsl.SqrtFloat64Eps {
			/* x is near -1 */
			for ell = 2; ell <= lmax; ell++ {
				fell := float64(ell)
				sgn := -1.0 /* derivative is odd in x for even ell */
				if gsl.IsOdd(ell) {
					sgn = 1.0
				}
				pre := sgn * 0.5 * fell * (fell + 1.0)
				result_deriv_array[ell] = pre * (1.0 - 0.25*(1.0+x)*(fell+2.0)*(fell-1.0))
			}
		} else {
			diff_a := 1.0 + x
			diff_b := 1.0 - x
			for ell = 2; ell <= lmax; ell++ {
				result_deriv_array[ell] = -float64(ell) * (x*result_array[ell] - result_array[ell-1]) / (diff_a * diff_b)
			}
		}

		return nil
	} else {
		return stat_array
	}
}

func Legendre_Plm_e(l, m int, x float64, result *Result) err.GSLError {
	/* If l is large and m is large, then we have to worry
	 * about overflow. Calculate an approximate exponent which
	 * measures the normalization of this thing.
	 */
	dif := float64(l - m)
	sum := float64(l + m)
	t_d := 0.5 * dif * (math.Log(dif) - 1.0)
	if dif == 0.0 {
		t_d = 0.0
	}
	t_s := 0.5 * sum * (math.Log(sum) - 1.0)
	if dif == 0.0 {
		t_s = 0.0
	}
	exp_check := 0.5*math.Log(2.0*float64(l)+1.0) + t_d - t_s

	if m < 0 || l < m || x < -1.0 || x > 1.0 {
		return DomainError(result)
	} else if exp_check < gsl.LnMinFloat64+10.0 {
		/* Bail out. */
		return OverflowError(result)
	} else {
		/* Account for the error due to the
		 * representation of 1-x.
		 */
		err_amp := 1.0 / (gsl.Float64Eps + math.Abs(1.0-math.Abs(x)))

		/* P_m^m(x) and P_{m+1}^m(x) */
		p_mm := legendre_Pmm(m, x)
		p_mmp1 := x * (2*float64(m) + 1) * p_mm

		if l == m {
			result.val = p_mm
			result.err = err_amp * 2.0 * gsl.Float64Eps * math.Abs(p_mm)
			return nil
		} else if l == m+1 {
			result.val = p_mmp1
			result.err = err_amp * 2.0 * gsl.Float64Eps * math.Abs(p_mmp1)
			return nil
		} else {
			/* upward recurrence: (l-m) P(l,m) = (2l-1) z P(l-1,m) - (l+m-1) P(l-2,m)
			 * start at P(m,m), P(m+1,m)
			 */

			p_ellm2 := p_mm
			p_ellm1 := p_mmp1
			p_ell := 0.0
			var ell int

			for ell = m + 2; ell <= l; ell++ {
				p_ell = (x*(2*float64(ell)-1)*p_ellm1 - float64(ell+m-1)*p_ellm2) / float64(ell-m)
				p_ellm2 = p_ellm1
				p_ellm1 = p_ell
			}

			result.val = p_ell
			result.err = err_amp * (0.5*float64(l-m) + 1.0) * gsl.Float64Eps * math.Abs(p_ell)

			return nil
		}
	}
}

func Legendre_sphPlm_e(l, m int, x float64, result *Result) err.GSLError {

	if m < 0 || l < m || x < -1.0 || x > 1.0 {
		return DomainError(result)
	} else if m == 0 {
		var P Result
		stat_P := Legendre_Pl_e(l, x, &P)
		pre := math.Sqrt((2.0*float64(l) + 1.0) / (4.0 * gsl.Pi))
		result.val = pre * P.val
		result.err = pre * P.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return stat_P
	} else if x == 1.0 || x == -1.0 {
		/* m > 0 here */
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		/* m > 0 and |x| < 1 here */

		/* Starting value for recursion.
		 * Y_m^m(x) = math.Sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) ) (-1)^m (1-x^2)^(m/2) / pi^(1/4)
		 */
		var (
			fm                 = float64(m)
			lncirc             Result
			lnpoch             Result
			lnpre_val          float64
			lnpre_err          float64
			ex_pre             Result
			sr                 float64
			sgn                = 1.0
			y_mmp1_factor      = x * math.Sqrt(2.0*fm+3.0)
			y_mm, y_mm_err     float64
			y_mmp1, y_mmp1_err float64
		)
		if gsl.IsOdd(m) {
			sgn = -1.0
		}
		Log_1plusx_e(-x*x, &lncirc)
		Lnpoch_e(fm, 0.5, &lnpoch) /* Gamma(m+1/2)/Gamma(m) */
		lnpre_val = -0.25*gsl.LnPi + 0.5*(lnpoch.val+fm*lncirc.val)
		lnpre_err = 0.25*gsl.LnPi*gsl.Float64Eps + 0.5*(lnpoch.err+math.Abs(fm)*lncirc.err)
		/* Compute math.Exp(ln_pre) with error term, avoiding call to gsl_sf_exp_err BJG */
		ex_pre.val = math.Exp(lnpre_val)
		ex_pre.err = 2.0 * (math.Sinh(lnpre_err) + gsl.Float64Eps) * ex_pre.val
		sr = math.Sqrt((2.0 + 1.0/fm) / (4.0 * gsl.Pi))
		y_mm = sgn * sr * ex_pre.val
		y_mm_err = 2.0*gsl.Float64Eps*math.Abs(y_mm) + sr*ex_pre.err
		y_mm_err *= 1.0 + 1.0/(gsl.Float64Eps+math.Abs(1.0-x))
		y_mmp1 = y_mmp1_factor * y_mm
		y_mmp1_err = math.Abs(y_mmp1_factor) * y_mm_err

		if l == m {
			result.val = y_mm
			result.err = y_mm_err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(y_mm)
			return nil
		} else if l == m+1 {
			result.val = y_mmp1
			result.err = y_mmp1_err
			result.err += 2.0 * gsl.Float64Eps * math.Abs(y_mmp1)
			return nil
		} else {
			y_ell := 0.0
			y_ell_err := 0.0
			var ell int

			/* Compute Y_l^m, l > m+1, upward recursion on l. */
			for ell = m + 2; ell <= l; ell++ {
				rat1 := float64(ell-m) / float64(ell+m)
				rat2 := float64(ell-m-1) / float64(ell+m-1)
				factor1 := math.Sqrt(rat1 * (2.0*float64(ell) + 1.0) * (2.0*float64(ell) - 1.0))
				factor2 := math.Sqrt(rat1 * rat2 * (2.0*float64(ell) + 1.0) / (2.0*float64(ell) - 3.0))
				y_ell = (x*y_mmp1*factor1 - float64(ell+m-1)*y_mm*factor2) / float64(ell-m)
				y_mm = y_mmp1
				y_mmp1 = y_ell

				y_ell_err = 0.5 * (math.Abs(x*factor1)*y_mmp1_err + math.Abs(float64(ell+m-1)*factor2)*y_mm_err) / math.Abs(float64(ell-m))
				y_mm_err = y_mmp1_err
				y_mmp1_err = y_ell_err
			}

			result.val = y_ell
			result.err = y_ell_err + (0.5*float64(l-m)+1.0)*gsl.Float64Eps*math.Abs(y_ell)

			return nil
		}
	}
}

func Legendre_Plm_array(lmax, m int, x float64, result_array []float64) err.GSLError {
	/* If l is large and m is large, then we have to worry
	 * about overflow. Calculate an approximate exponent which
	 * measures the normalization of this thing.
	 */

	dif := float64(lmax - m)
	sum := float64(lmax + m)
	t_d := 0.0
	if dif == 0.0 {
		t_d = 0.5 * dif * (math.Log(dif) - 1.0)
	}
	t_s := 0.0
	if dif == 0.0 {
		t_d = 0.5 * sum * (math.Log(sum) - 1.0)
	}
	exp_check := 0.5*math.Log(2.0*float64(lmax)+1.0) + t_d - t_s

	if m < 0 || lmax < m || x < -1.0 || x > 1.0 {
		return err.ERROR("domain error", err.EDOM)
	} else if m > 0 && (x == 1.0 || x == -1.0) {
		var ell int
		for ell = m; ell <= lmax; ell++ {
			result_array[ell-m] = 0.0
		}
		return nil
	} else if exp_check < gsl.LnMinFloat64+10.0 {
		/* Bail out. */
		return err.ERROR("overflow", err.EOVRFLW)
	} else {
		p_mm := legendre_Pmm(m, x)
		p_mmp1 := x * (2.0*float64(m) + 1.0) * p_mm

		if lmax == m {
			result_array[0] = p_mm
			return nil
		} else if lmax == m+1 {
			result_array[0] = p_mm
			result_array[1] = p_mmp1
			return nil
		} else {
			p_ellm2 := p_mm
			p_ellm1 := p_mmp1
			p_ell := 0.0
			var ell int

			result_array[0] = p_mm
			result_array[1] = p_mmp1

			for ell = m + 2; ell <= lmax; ell++ {
				p_ell = (x*(2.0*float64(ell)-1.0)*p_ellm1 - float64(ell+m-1)*p_ellm2) / float64(ell-m)
				p_ellm2 = p_ellm1
				p_ellm1 = p_ell
				result_array[ell-m] = p_ell
			}

			return nil
		}
	}
}

func Legendre_Plm_deriv_array(lmax, m int, x float64, result_array, result_deriv_array []float64) err.GSLError {
	if m < 0 || m > lmax {
		return err.ERROR("m < 0 or m > lmax", err.EDOM)
	} else if m == 0 {
		/* It is better to do m=0 this way, so we can more easily
		 * trap the divergent case which can occur when m == 1.
		 */
		return Legendre_Pl_deriv_array(lmax, x, result_array, result_deriv_array)
	} else {
		stat_array := Legendre_Plm_array(lmax, m, x, result_array)

		if stat_array == nil {
			var ell int

			if m == 1 && (1.0-math.Abs(x) < gsl.Float64Eps) {
				/* This divergence is real and comes from the cusp-like
				 * behaviour for m = 1. For example, P[1,1] = - Sqrt[1-x^2].
				 */
				return err.ERROR("divergence near |x| = 1.0 since m = 1", err.EOVRFLW)
			} else if m == 2 && (1.0-math.Abs(x) < gsl.Float64Eps) {
				/* m = 2 gives a finite nonzero result for |x| near 1 */
				if math.Abs(x-1.0) < gsl.Float64Eps {
					for ell = m; ell <= lmax; ell++ {
						fell := float64(ell)
						result_deriv_array[ell-m] = -0.25 * x * (fell - 1.0) * fell * (fell + 1.0) * (fell + 2.0)
					}
				} else if math.Abs(x+1.0) < gsl.Float64Eps {
					for ell = m; ell <= lmax; ell++ {
						fell := float64(ell)
						sgn := -1.0
						if gsl.IsOdd(ell) {
							sgn = 1.0
						}
						result_deriv_array[ell-m] = -0.25 * sgn * x * (fell - 1.0) * fell * (fell + 1.0) * (fell + 2.0)
					}
				}
				return nil
			} else {
				/* m > 2 is easier to deal with since the endpoints always vanish */
				if 1.0-math.Abs(x) < gsl.Float64Eps {
					for ell = m; ell <= lmax; ell++ {
						result_deriv_array[ell-m] = 0.0
					}
					return nil
				} else {
					diff_a := 1.0 + x
					diff_b := 1.0 - x
					result_deriv_array[0] = -float64(m) * x / (diff_a * diff_b) * result_array[0]
					if lmax-m >= 1 {
						result_deriv_array[1] = float64(2*m+1) * (x*result_deriv_array[0] + result_array[0])
					}
					for ell = m + 2; ell <= lmax; ell++ {
						result_deriv_array[ell-m] = -(float64(ell)*x*result_array[ell-m] - float64(ell+m)*result_array[ell-1-m]) / (diff_a * diff_b)
					}
					return nil
				}
			}
		} else {
			return stat_array
		}
	}
}

func Legendre_sphPlm_array(lmax, m int, x float64, result_array []float64) err.GSLError {
	if m < 0 || lmax < m || x < -1.0 || x > 1.0 {
		return err.ERROR("error", err.EDOM)
	} else if m > 0 && (x == 1.0 || x == -1.0) {
		var ell int
		for ell = m; ell <= lmax; ell++ {
			result_array[ell-m] = 0.0
		}
		return nil
	} else {
		var y_mm, y_mmp1 float64

		if m == 0 {
			y_mm = 0.5 / gsl.SqrtPi /* Y00 = 1/math.Sqrt(4pi) */
			y_mmp1 = x * gsl.Sqrt3 * y_mm
		} else {
			/* |x| < 1 here */

			var lncirc, lnpoch Result
			var lnpre float64
			sgn := 1.0
			if gsl.IsOdd(m) {
				sgn = -1.0
			}
			Log_1plusx_e(-x*x, &lncirc)
			Lnpoch_e(float64(m), 0.5, &lnpoch) /* Gamma(m+1/2)/Gamma(m) */
			lnpre = -0.25*gsl.LnPi + 0.5*(lnpoch.val+float64(m)*lncirc.val)
			y_mm = math.Sqrt((2.0+1.0/float64(m))/(4.0*gsl.Pi)) * sgn * math.Exp(lnpre)
			y_mmp1 = x * math.Sqrt(2.0*float64(m)+3.0) * y_mm
		}

		if lmax == m {
			result_array[0] = y_mm
			return nil
		} else if lmax == m+1 {
			result_array[0] = y_mm
			result_array[1] = y_mmp1
			return nil
		} else {
			var y_ell float64
			var ell int

			result_array[0] = y_mm
			result_array[1] = y_mmp1

			/* Compute Y_l^m, l > m+1, upward recursion on l. */
			for ell = m + 2; ell <= lmax; ell++ {
				fell := float64(ell)
				fm := float64(m)
				rat1 := float64(ell-m) / float64(ell+m)
				rat2 := (fell - fm - 1.0) / (fell + fm - 1.0)
				factor1 := math.Sqrt(rat1 * (2*fell + 1) * (2*fell - 1))
				factor2 := math.Sqrt(rat1 * rat2 * (2*fell + 1) / (2*fell - 3))
				y_ell = (x*y_mmp1*factor1 - (fell+fm-1)*y_mm*factor2) / (fell - fm)
				y_mm = y_mmp1
				y_mmp1 = y_ell
				result_array[ell-m] = y_ell
			}
		}

		return nil
	}
}

func Legendre_sphPlm_deriv_array(lmax, m int, x float64, result_array, result_deriv_array []float64) err.GSLError {
	if m < 0 || lmax < m || x < -1.0 || x > 1.0 {
		return err.ERROR("domain", err.EDOM)
	} else if m == 0 {
		/* m = 0 is easy to trap */
		stat_array := Legendre_Pl_deriv_array(lmax, x, result_array, result_deriv_array)
		var ell int
		for ell = 0; ell <= lmax; ell++ {
			prefactor := math.Sqrt((2.0*float64(ell) + 1.0) / (4.0 * gsl.Pi))
			result_array[ell] *= prefactor
			result_deriv_array[ell] *= prefactor
		}
		return stat_array
	} else if m == 1 {
		/* Trapping m = 1 is necessary because of the possible divergence.
		 * Recall that this divergence is handled properly in ..._Plm_deriv_array(),
		 * and the scaling factor is not large for small m, so we just scale.
		 */
		stat_array := Legendre_Plm_deriv_array(lmax, m, x, result_array, result_deriv_array)
		var ell int
		for ell = 1; ell <= lmax; ell++ {
			prefactor := math.Sqrt((2.0*float64(ell) + 1.0) / (float64(ell) + 1.0) / (4.0 * gsl.Pi * float64(ell)))
			result_array[ell-1] *= prefactor
			result_deriv_array[ell-1] *= prefactor
		}
		return stat_array
	} else {
		/* as for the derivative of P_lm, everything is regular for m >= 2 */

		stat_array := Legendre_sphPlm_array(lmax, m, x, result_array)

		if stat_array == nil {
			var ell int

			if 1.0-math.Abs(x) < gsl.Float64Eps {
				for ell = m; ell <= lmax; ell++ {
					result_deriv_array[ell-m] = 0.0
				}
				return nil
			} else {
				diff_a := 1.0 + x
				diff_b := 1.0 - x
				result_deriv_array[0] = -float64(m) * x / (diff_a * diff_b) * result_array[0]
				if lmax-m >= 1 {
					result_deriv_array[1] = math.Sqrt(2.0*float64(m)+3.0) * (x*result_deriv_array[0] + result_array[0])
				}
				for ell = m + 2; ell <= lmax; ell++ {
					fell := float64(ell)
					c1 := math.Sqrt(((2.0*fell + 1.0) / (2.0*fell - 1.0)) * (float64(ell-m) / float64(ell+m)))
					result_deriv_array[ell-m] = -(fell*x*result_array[ell-m] - c1*float64(ell+m)*result_array[ell-1-m]) / (diff_a * diff_b)
				}
				return nil
			}
		} else {
			return stat_array
		}
	}
}

func Legendre_array_size(lmax, m int) int {
	return lmax - m + 1
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Legendre_P1(x float64) float64 {
	result := new(Result)
	status := Legendre_P1_e(x, result)
	return EvalResult(result, status)
}

func Legendre_P2(x float64) float64 {
	result := new(Result)
	status := Legendre_P2_e(x, result)
	return EvalResult(result, status)
}

func Legendre_P3(x float64) float64 {
	result := new(Result)
	status := Legendre_P3_e(x, result)
	return EvalResult(result, status)
}

func Legendre_Pl(l int, x float64) float64 {
	result := new(Result)
	status := Legendre_Pl_e(l, x, result)
	return EvalResult(result, status)
}

func Legendre_Plm(l, m int, x float64) float64 {
	result := new(Result)
	status := Legendre_Plm_e(l, m, x, result)
	return EvalResult(result, status)
}

func Legendre_sphPlm(l, m int, x float64) float64 {
	result := new(Result)
	status := Legendre_sphPlm_e(l, m, x, result)
	return EvalResult(result, status)
}
