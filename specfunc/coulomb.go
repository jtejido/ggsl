/* specfunc/coulomb.c
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

/* the L=0 normalization constant
 * [Abramowitz+Stegun 14.1.8]
 */
func c0sq(eta float64) float64 {
	twopieta := 2.0 * gsl.Pi * eta

	if math.Abs(eta) < gsl.Float64Eps {
		return 1.0
	} else if twopieta > gsl.LnMaxFloat64 {
		return 0.0
	} else {
		var scale Result
		Expm1_e(twopieta, &scale)
		return twopieta / scale.val
	}
}

/* the full definition of C_L(eta) for any valid L and eta
 * [Abramowitz and Stegun 14.1.7]
 * This depends on the complex gamma function. For large
 * arguments the phase of the complex gamma function is not
 * very accurately determined. However the modulus is, and that
 * is all that we need to calculate C_L.
 *
 * This is not valid for L <= -3/2  or  L = -1.
 */
func cLeta(L, eta float64, result *Result) err.GSLError {
	var ln1 Result /* log of numerator Gamma function */
	var ln2 Result /* log of denominator Gamma function */
	sgn := 1.0
	var arg_val, arg_err float64

	if math.Abs(eta/(L+1.0)) < gsl.Float64Eps {
		Lngamma_e(L+1.0, &ln1)
	} else {
		var p1 Result                            /* phase of numerator Gamma -- not used */
		Lngamma_complex_e(L+1.0, eta, &ln1, &p1) /* should be ok */
	}

	Lngamma_e(2.0*(L+1.0), &ln2)
	if L < -1.0 {
		sgn = -sgn
	}

	arg_val = L*gsl.Ln2 - 0.5*eta*gsl.Pi + ln1.val - ln2.val
	arg_err = ln1.err + ln2.err
	arg_err += gsl.Float64Eps * (math.Abs(L*gsl.Ln2) + math.Abs(0.5*eta*gsl.Pi))
	return Exp_err_e(arg_val, arg_err, result)
}

func Coulomb_CL_e(lam, eta float64, result *Result) err.GSLError {
	if lam <= -1.0 {
		return DomainError(result)
	} else if math.Abs(lam) < gsl.Float64Eps {
		/* saves a calculation of complex_lngamma(), otherwise not necessary */
		result.val = math.Sqrt(c0sq(eta))
		result.err = 2.0 * gsl.Float64Eps * result.val
		return nil
	} else {
		return cLeta(lam, eta, result)
	}
}

/* cl[0] .. cl[kmax] = C_{lam_min}(eta) .. C_{lam_min+kmax}(eta)
 */
func Coulomb_CL_array(lam_min float64, kmax int, eta float64, cl []float64) err.GSLError {
	var cl_0 Result
	Coulomb_CL_e(lam_min, eta, &cl_0)
	cl[0] = cl_0.val

	for k := 1; k <= kmax; k++ {
		L := lam_min + float64(k)
		cl[k] = cl[k-1] * math.Hypot(L, eta) / (L * (2.0*L + 1.0))
	}

	return nil
}

/* Determine the connection phase, phi_lambda.
 * See coulomb_FG_series() below. We have
 * to be careful about math.Sin(phi).0. Note that
 * there is an underflow condition for large
 * positive eta in any case.
 */
func coulomb_connection(lam, eta float64, cos_phi, sin_phi *float64) err.GSLError {
	if eta > -gsl.LnMinFloat64/2.0*gsl.Pi-1.0 {
		*cos_phi = 1.0
		*sin_phi = 0.0
		return err.ERROR("error", err.EUNDRFLW)
	} else if eta > -gsl.LnFloat64Eps/(4.0*gsl.Pi) {
		eps := 2.0 * math.Exp(-2.0*gsl.Pi*eta)
		tpl := math.Tan(gsl.Pi * lam)
		dth := eps * tpl / (tpl*tpl + 1.0)
		*cos_phi = -1.0 + 0.5*dth*dth
		*sin_phi = -dth
		return nil
	} else {
		X := math.Tanh(gsl.Pi*eta) / math.Tan(gsl.Pi*lam)
		phi := -math.Atan(X) - (lam+0.5)*gsl.Pi
		*cos_phi = math.Cos(phi)
		*sin_phi = math.Sin(phi)
		return nil
	}
}

/* Evaluate the Frobenius series for F_lam(eta,x) and G_lam(eta,x).
 * Homegrown algebra. Evaluates the series for F_{lam} and
 * F_{-lam-1}, then uses
 *    G_{lam} = (F_{lam} math.Cos(phi) - F_{-lam-1}) / math.Sin(phi)
 * where
 *    phi = Arg[Gamma[1+lam+I eta]] - Arg[Gamma[-lam + I eta]] - (lam+1/2)Pi
 *        = Arg[Sin[Pi(-lam+I eta)] - (lam+1/2)Pi
 *        = atan2(-math.Cos(lam Pi)sinh(eta Pi), -math.Sin(lam Pi)cosh(eta Pi)) - (lam+1/2)Pi
 *
 *        = -amath.Tan(X) - (lam+1/2) Pi,  X = tanh(eta Pi)/math.Tan(lam Pi)
 *
 * Not appropriate for lam <= -1/2, lam = 0, or lam >= 1/2.
 */
func coulomb_FG_series(lam, eta, x float64, F, G *Result) err.GSLError {

	var (
		max_iter                 = 800
		ClamA, ClamB             Result
		stat_A                   = cLeta(lam, eta, &ClamA)
		stat_B                   = cLeta(-lam-1.0, eta, &ClamB)
		tlp1                     = 2.0*lam + 1.0
		pow_x                    = math.Pow(x, lam)
		cos_phi_lam, sin_phi_lam float64
		uA_mm2                   = 1.0 /* uA sum is for F_{lam} */
		uA_mm1                   = x * eta / (lam + 1.0)
		uA_m                     float64
		uB_mm2                   = 1.0 /* uB sum is for F_{-lam-1} */
		uB_mm1                   = -x * eta / lam
		uB_m                     float64
		A_sum                    = uA_mm2 + uA_mm1
		B_sum                    = uB_mm2 + uB_mm1
		A_abs_del_prev           = math.Abs(A_sum)
		B_abs_del_prev           = math.Abs(B_sum)
		FA, FB                   Result
		m                        = 2
		stat_conn                = coulomb_connection(lam, eta, &cos_phi_lam, &sin_phi_lam)
	)

	if stat_conn.Status() == err.EUNDRFLW {
		F.val = 0.0 /* FIXME: should this be set to Inf too like G? */
		F.err = 0.0
		return OverflowError(G)
	}

	for m < max_iter {
		var abs_dA, abs_dB float64
		fm := float64(m)
		uA_m = x * (2.0*eta*uA_mm1 - x*uA_mm2) / (fm * (fm + tlp1))
		uB_m = x * (2.0*eta*uB_mm1 - x*uB_mm2) / (fm * (fm - tlp1))
		A_sum += uA_m
		B_sum += uB_m
		abs_dA = math.Abs(uA_m)
		abs_dB = math.Abs(uB_m)
		if m > 15 {
			/* Don't bother checking until we have gone out a little ways;
			 * a minor optimization. Also make sure to check both the
			 * current and the previous increment because the odd and even
			 * terms of the sum can have very different behaviour, depending
			 * on the value of eta.
			 */
			max_abs_dA := math.Max(abs_dA, A_abs_del_prev)
			max_abs_dB := math.Max(abs_dB, B_abs_del_prev)
			abs_A := math.Abs(A_sum)
			abs_B := math.Abs(B_sum)
			if max_abs_dA/(max_abs_dA+abs_A) < 4.0*gsl.Float64Eps && max_abs_dB/(max_abs_dB+abs_B) < 4.0*gsl.Float64Eps {
				break
			}
		}
		A_abs_del_prev = abs_dA
		B_abs_del_prev = abs_dB
		uA_mm2 = uA_mm1
		uA_mm1 = uA_m
		uB_mm2 = uB_mm1
		uB_mm1 = uB_m
		m++
	}

	FA.val = A_sum * ClamA.val * pow_x * x
	FA.err = math.Abs(A_sum)*ClamA.err*pow_x*x + 2.0*gsl.Float64Eps*math.Abs(FA.val)
	FB.val = B_sum * ClamB.val / pow_x
	FB.err = math.Abs(B_sum)*ClamB.err/pow_x + 2.0*gsl.Float64Eps*math.Abs(FB.val)

	F.val = FA.val
	F.err = FA.err

	G.val = (FA.val*cos_phi_lam - FB.val) / sin_phi_lam
	G.err = (FA.err*math.Abs(cos_phi_lam) + FB.err) / math.Abs(sin_phi_lam)

	if m >= max_iter {
		return err.ERROR("error", err.EMAXITER)
	}

	return err.ErrorSelect(stat_A, stat_B)
}

/* Evaluate the Frobenius series for F_0(eta,x) and G_0(eta,x).
 * See [Bardin et al., CPC 3, 73 (1972), (14)-(17)];
 * note the misprint in (17): nu_0=1 is correct, not nu_0=0.
 */
func coulomb_FG0_series(eta, x float64, F, G *Result) err.GSLError {
	var (
		max_iter       = 800
		x2             = x * x
		tex            = 2.0 * eta * x
		C0             Result
		stat_CL        = cLeta(0.0, eta, &C0)
		r1pie          Result
		psi_stat       = Psi_1piy_e(eta, &r1pie)
		u_mm2          = 0.0 /* u_0 */
		u_mm1          = x   /* u_1 */
		u_m            float64
		v_mm2          = 1.0                                     /* nu_0 */
		v_mm1          = tex * (2.0*gsl.Euler - 1.0 + r1pie.val) /* nu_1 */
		v_m            float64
		u_sum          = u_mm2 + u_mm1
		v_sum          = v_mm2 + v_mm1
		u_abs_del_prev = math.Abs(u_sum)
		v_abs_del_prev = math.Abs(v_sum)
		m              = 2
		u_sum_err      = 2.0 * gsl.Float64Eps * math.Abs(u_sum)
		v_sum_err      = 2.0 * gsl.Float64Eps * math.Abs(v_sum)
		ln2x           = math.Log(2.0 * x)
	)

	for m < max_iter {
		var abs_du, abs_dv float64
		fm := float64(m)
		m_mm1 := fm * (fm - 1.0)
		u_m = (tex*u_mm1 - x2*u_mm2) / m_mm1
		v_m = (tex*v_mm1 - x2*v_mm2 - 2.0*eta*(2*fm-1)*u_m) / m_mm1
		u_sum += u_m
		v_sum += v_m
		abs_du = math.Abs(u_m)
		abs_dv = math.Abs(v_m)
		u_sum_err += 2.0 * gsl.Float64Eps * abs_du
		v_sum_err += 2.0 * gsl.Float64Eps * abs_dv
		if m > 15 {
			/* Don't bother checking until we have gone out a little ways;
			 * a minor optimization. Also make sure to check both the
			 * current and the previous increment because the odd and even
			 * terms of the sum can have very different behaviour, depending
			 * on the value of eta.
			 */
			max_abs_du := math.Max(abs_du, u_abs_del_prev)
			max_abs_dv := math.Max(abs_dv, v_abs_del_prev)
			abs_u := math.Abs(u_sum)
			abs_v := math.Abs(v_sum)
			if max_abs_du/(max_abs_du+abs_u) < 40.0*gsl.Float64Eps && max_abs_dv/(max_abs_dv+abs_v) < 40.0*gsl.Float64Eps {
				break
			}
		}
		u_abs_del_prev = abs_du
		v_abs_del_prev = abs_dv
		u_mm2 = u_mm1
		u_mm1 = u_m
		v_mm2 = v_mm1
		v_mm1 = v_m
		m++
	}

	F.val = C0.val * u_sum
	F.err = C0.err * math.Abs(u_sum)
	F.err += math.Abs(C0.val) * u_sum_err
	F.err += 2.0 * gsl.Float64Eps * math.Abs(F.val)

	G.val = (v_sum + 2.0*eta*u_sum*ln2x) / C0.val
	G.err = (math.Abs(v_sum) + math.Abs(2.0*eta*u_sum*ln2x)) / math.Abs(C0.val) * math.Abs(C0.err/C0.val)
	G.err += (v_sum_err + math.Abs(2.0*eta*u_sum_err*ln2x)) / math.Abs(C0.val)
	G.err += 2.0 * gsl.Float64Eps * math.Abs(G.val)

	if m == max_iter {
		return err.ERROR("error", err.EMAXITER)
	}

	return err.ErrorSelect(psi_stat, stat_CL)
}

/* Evaluate the Frobenius series for F_{-1/2}(eta,x) and G_{-1/2}(eta,x).
 * Homegrown algebra.
 */
func coulomb_FGmhalf_series(eta, x float64, F, G *Result) err.GSLError {
	var (
		max_iter            = 800
		rx                  = math.Sqrt(x)
		x2                  = x * x
		tex                 = 2.0 * eta * x
		Cmhalf              Result
		stat_CL             = cLeta(-0.5, eta, &Cmhalf)
		u_mm2               = 1.0         /* u_0 */
		u_mm1               = tex * u_mm2 /* u_1 */
		u_m                 float64
		v_mm2, v_mm1, v_m   float64
		f_sum, g_sum        float64
		tmp1                float64
		rpsi_1pe, rpsi_1p2e Result
		m                   = 2
	)

	Psi_1piy_e(eta, &rpsi_1pe)
	Psi_1piy_e(2.0*eta, &rpsi_1p2e)

	v_mm2 = 2.0*gsl.Euler - gsl.Ln2 - rpsi_1pe.val + 2.0*rpsi_1p2e.val
	v_mm1 = tex * (v_mm2 - 2.0*u_mm2)

	f_sum = u_mm2 + u_mm1
	g_sum = v_mm2 + v_mm1

	for m < max_iter {
		fm := float64(m)
		m2 := fm * fm
		u_m = (tex*u_mm1 - x2*u_mm2) / m2
		v_m = (tex*v_mm1 - x2*v_mm2 - 2.0*fm*u_m) / m2
		f_sum += u_m
		g_sum += v_m
		if f_sum != 0.0 && g_sum != 0.0 && (math.Abs(u_m/f_sum)+math.Abs(v_m/g_sum) < 10.0*gsl.Float64Eps) {
			break
		}
		u_mm2 = u_mm1
		u_mm1 = u_m
		v_mm2 = v_mm1
		v_mm1 = v_m
		m++
	}

	F.val = Cmhalf.val * rx * f_sum
	F.err = Cmhalf.err*math.Abs(rx*f_sum) + 2.0*gsl.Float64Eps*math.Abs(F.val)

	tmp1 = f_sum * math.Log(x)
	G.val = -rx * (tmp1 + g_sum) / Cmhalf.val
	G.err = math.Abs(rx) * (math.Abs(tmp1) + math.Abs(g_sum)) / math.Abs(Cmhalf.val) * math.Abs(Cmhalf.err/Cmhalf.val)

	if m == max_iter {
		return err.ERROR("error", err.EMAXITER)
	}

	return stat_CL
}

/* Evolve the backwards recurrence for F,F'.
 *
 *    F_{lam-1}  = (S_lam F_lam + F_lam') / R_lam
 *    F_{lam-1}' = (S_lam F_{lam-1} - R_lam F_lam)
 * where
 *    R_lam = math.Sqrt(1 + (eta/lam)^2)
 *    S_lam = lam/x + eta/lam
 *
 */
func coulomb_F_recur(lam_min float64, kmax int, eta, x, F_lam_max, Fp_lam_max float64, F_lam_min, Fp_lam_min *float64) err.GSLError {
	var (
		x_inv   = 1.0 / x
		fcl     = F_lam_max
		fpl     = Fp_lam_max
		lam_max = lam_min + float64(kmax)
		lam     = lam_max
		k       int
	)

	for k = kmax - 1; k >= 0; k-- {
		el := eta / lam
		rl := math.Hypot(1.0, el)
		sl := el + lam*x_inv
		var fc_lm1 float64
		fc_lm1 = (fcl*sl + fpl) / rl
		fpl = fc_lm1*sl - fcl*rl
		fcl = fc_lm1
		lam -= 1.0
	}

	*F_lam_min = fcl
	*Fp_lam_min = fpl
	return nil
}

/* Evolve the forward recurrence for G,G'.
 *
 *   G_{lam+1}  = (S_lam G_lam - G_lam')/R_lam
 *   G_{lam+1}' = R_{lam+1} G_lam - S_lam G_{lam+1}
 *
 * where S_lam and R_lam are as above in the F recursion.
 */
func coulomb_G_recur(lam_min float64, kmax int, eta, x, G_lam_min, Gp_lam_min float64, G_lam_max, Gp_lam_max *float64) err.GSLError {
	var (
		x_inv = 1.0 / x
		gcl   = G_lam_min
		gpl   = Gp_lam_min
		lam   = lam_min + 1.0
		k     int
	)

	for k = 1; k <= kmax; k++ {
		el := eta / lam
		rl := math.Hypot(1.0, el)
		sl := el + lam*x_inv
		gcl1 := (sl*gcl - gpl) / rl
		gpl = rl*gcl - sl*gcl1
		gcl = gcl1
		lam += 1.0
	}

	*G_lam_max = gcl
	*Gp_lam_max = gpl
	return nil
}

/* Evaluate the first continued fraction, giving
 * the ratio F'/F at the upper lambda value.
 * We also determine the sign of F at that point,
 * since it is the sign of the last denominator
 * in the continued fraction.
 */
func coulomb_CF1(lambda, eta, x float64, fcl_sign, result *float64, count *int) err.GSLError {
	var (
		CF1_small = 1.e-30
		CF1_abort = 1.0e+05
		CF1_acc   = 2.0 * gsl.Float64Eps
		x_inv     = 1.0 / x
		px        = lambda + 1.0 + CF1_abort

		pk   = lambda + 1.0
		F    = eta/pk + pk*x_inv
		D, C float64
		df   float64
	)

	*fcl_sign = 1.0
	*count = 0

	if math.Abs(F) < CF1_small {
		F = CF1_small
	}
	D = 0.0
	C = F

	for math.Abs(df-1.0) > CF1_acc {
		pk1 := pk + 1.0
		ek := eta / pk
		rk2 := 1.0 + ek*ek
		tk := (pk + pk1) * (x_inv + ek/pk1)
		D = tk - rk2*D
		C = tk - rk2/C
		if math.Abs(C) < CF1_small {
			C = CF1_small
		}
		if math.Abs(D) < CF1_small {
			D = CF1_small
		}
		D = 1.0 / D
		df = D * C
		F = F * df
		if D < 0.0 {
			/* sign of result depends on sign of denominator */
			*fcl_sign = -*fcl_sign
		}
		pk = pk1
		if pk > px {
			*result = F
			return err.ERROR("error", err.ERUNAWAY)
		}
		*count++
	}

	*result = F
	return nil
}

/* Evaluate the second continued fraction to
 * obtain the ratio
 *    (G' + i F')/(G + i F) := P + i Q
 * at the specified lambda value.
 */
func coulomb_CF2(lambda, eta, x float64, result_P, result_Q *float64, count *int) err.GSLError {
	var (
		status     err.GSLError
		CF2_acc    = 4.0 * gsl.Float64Eps
		CF2_abort  = 2.0e+05
		wi         = 2.0 * eta
		x_inv      = 1.0 / x
		e2mm1      = eta*eta + lambda*(lambda+1.0)
		ar         = -e2mm1
		ai         = eta
		br         = 2.0 * (x - eta)
		bi         = 2.0
		dr         = br / (br*br + bi*bi)
		di         = -bi / (br*br + bi*bi)
		dp         = -x_inv * (ar*di + ai*dr)
		dq         = x_inv * (ar*dr - ai*di)
		A, B, C, D float64
		pk         = 0.0
		P          = 0.0
		Q          = 1.0 - eta*x_inv
	)

	*count = 0

	for math.Abs(dp)+math.Abs(dq) > (math.Abs(P)+math.Abs(Q))*CF2_acc {
		P += dp
		Q += dq
		pk += 2.0
		ar += pk
		ai += wi
		bi += 2.0
		D = ar*dr - ai*di + br
		di = ai*dr + ar*di + bi
		C = 1.0 / (D*D + di*di)
		dr = C * D
		di = -C * di
		A = br*dr - bi*di - 1.
		B = bi*dr + br*di
		C = dp*A - dq*B
		dq = dp*B + dq*A
		dp = C
		if pk > CF2_abort {
			status = err.RunAway()
			break
		}

		*count++
	}

	if Q < CF2_abort*gsl.Float64Eps*math.Abs(P) {
		status = err.Loss()
	}

	*result_P = P
	*result_Q = Q
	return status
}

/* WKB evaluation of F, G. Assumes  0 < x < turning point.
 * Overflows are trapped, GSL_EOVRFLW is signalled,
 * and an exponent is returned such that:
 *
 *   result_F = fjwkb * math.Exp(-exponent)
 *   result_G = gjwkb * math.Exp( exponent)
 *
 * See [Biedenharn et al. Phys. Rev. 97, 542-554 (1955), Section IV]
 *
 * Unfortunately, this is not very accurate in general. The
 * test cases typically have 3-4 digits of precision. One could
 * argue that this is ok for general use because, for instance,
 * F is exponentially small in this region and so the absolute
 * accuracy is still roughly acceptable. But it would be better
 * to have a systematic method for improving the precision. See
 * the Abad+Sesma method discussion below.
 */
func coulomb_jwkb(lam, eta, x float64, fjwkb, gjwkb *Result, exponent *float64) err.GSLError {
	var (
		llp1           = lam*(lam+1.0) + 6.0/35.0
		llp1_eff       = math.Max(llp1, 0.0)
		rho_ghalf      = math.Sqrt(x*(2.0*eta-x) + llp1_eff)
		sinh_arg       = math.Sqrt(llp1_eff/(eta*eta+llp1_eff)) * rho_ghalf / x
		sinh_inv       = math.Log(sinh_arg + math.Hypot(1.0, sinh_arg))
		phi            = math.Abs(rho_ghalf - eta*math.Atan2(rho_ghalf, x-eta) - math.Sqrt(llp1_eff)*sinh_inv)
		zeta_half      = math.Pow(3.0*phi/2.0, 1.0/3.0)
		prefactor      = math.Sqrt(gsl.Pi * phi * x / (6.0 * rho_ghalf))
		F              = prefactor * 3.0 / zeta_half
		G              = prefactor * 3.0 / zeta_half /* Note the math.Sqrt(3) from Bi normalization */
		F_exp          float64
		G_exp          float64
		airy_scale_exp = phi
		ai, bi         Result
	)

	Airy_Ai_scaled_e(zeta_half*zeta_half, gsl.MODE_DEFAULT, &ai)
	Airy_Bi_scaled_e(zeta_half*zeta_half, gsl.MODE_DEFAULT, &bi)
	F *= ai.val
	G *= bi.val
	F_exp = math.Log(F) - airy_scale_exp
	G_exp = math.Log(G) + airy_scale_exp

	if G_exp >= gsl.LnMaxFloat64 {
		fjwkb.val = F
		gjwkb.val = G
		fjwkb.err = 1.0e-3 * math.Abs(F) /* FIXME: real error here ... could be smaller */
		gjwkb.err = 1.0e-3 * math.Abs(G)
		*exponent = airy_scale_exp
		return err.ERROR("error", err.EOVRFLW)
	} else {
		fjwkb.val = math.Exp(F_exp)
		gjwkb.val = math.Exp(G_exp)
		fjwkb.err = 1.0e-3 * math.Abs(fjwkb.val)
		gjwkb.err = 1.0e-3 * math.Abs(gjwkb.val)
		*exponent = 0.0
		return nil
	}
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Coulomb_wave_FG_e(eta, x, lam_F float64, k_lam_G int, F, Fp, G, Gp *Result, exp_F, exp_G *float64) err.GSLError {
	lam_G := lam_F - float64(k_lam_G)

	if x < 0.0 || lam_F <= -0.5 || lam_G <= -0.5 {
		F.val = 0.0
		F.err = 0.0
		Fp.val = 0.0
		Fp.err = 0.0
		G.val = 0.0
		G.err = 0.0
		Gp.val = 0.0
		Gp.err = 0.0
		*exp_F = 0.0
		*exp_G = 0.0
		return err.ERROR("domain error", err.EDOM)
	} else if x == 0.0 {
		var C0 Result
		cLeta(0.0, eta, &C0)
		F.val = 0.0
		F.err = 0.0
		Fp.val = 0.0
		Fp.err = 0.0
		G.val = 0.0 /* FIXME: should be Inf */
		G.err = 0.0
		Gp.val = 0.0 /* FIXME: should be Inf */
		Gp.err = 0.0
		*exp_F = 0.0
		*exp_G = 0.0
		if lam_F == 0.0 {
			Fp.val = C0.val
			Fp.err = C0.err
		}
		if lam_G == 0.0 {
			Gp.val = 1.0 / C0.val
			Gp.err = math.Abs(C0.err/C0.val) / math.Abs(C0.val)
		}

		return err.ERROR("domain error", err.EDOM)
		/* After all, since we are asking for G, this is a domain error... */
	} else if x < 1.2 && 2.0*gsl.Pi*eta < 0.9*(-gsl.LnMinFloat64) && math.Abs(eta*x) < 10.0 {
		/* Reduce to a small lambda value and use the series
		 * representations for F and G. We cannot allow eta to
		 * be large and positive because the connection formula
		 * for G_lam is badly behaved due to an underflow in math.Sin(phi_lam)
		 * [see coulomb_FG_series() and coulomb_connection() above].
		 * Note that large negative eta is ok however.
		 */
		var (
			SMALL                               = gsl.SqrtFloat64Eps
			N                                   = int(lam_F + 0.5)
			span                                = int(gsl.Max(k_lam_G, N))
			lam_min                             = lam_F - float64(N) /* -1/2 <= lam_min < 1/2 */
			F_lam_F, Fp_lam_F                   float64
			G_lam_G                             = 0.0
			Gp_lam_G                            = 0.0
			F_lam_F_err, Fp_lam_F_err           float64
			Fp_over_F_lam_F                     float64
			F_sign_lam_F                        float64
			F_lam_min_unnorm, Fp_lam_min_unnorm float64
			Fp_over_F_lam_min                   float64
			F_lam_min                           Result
			G_lam_min, Gp_lam_min               Result
			F_scale                             float64
			Gerr_frac                           float64
			F_scale_frac_err                    float64
			F_unnorm_frac_err                   float64

			/* Determine F'/F at lam_F. */
			CF1_count int
			stat_CF1  = coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F, &CF1_count)

			stat_ser err.GSLError
			stat_Fr  err.GSLError
			stat_Gr  err.GSLError
		)

		/* Recurse down with unnormalized F,F' values. */
		F_lam_F = SMALL
		Fp_lam_F = Fp_over_F_lam_F * F_lam_F
		if span != 0 {
			stat_Fr = coulomb_F_recur(lam_min, span, eta, x, F_lam_F, Fp_lam_F, &F_lam_min_unnorm, &Fp_lam_min_unnorm)
		} else {
			F_lam_min_unnorm = F_lam_F
			Fp_lam_min_unnorm = Fp_lam_F
			stat_Fr = nil
		}

		/* Determine F and G at lam_min. */
		if lam_min == -0.5 {
			stat_ser = coulomb_FGmhalf_series(eta, x, &F_lam_min, &G_lam_min)
		} else if lam_min == 0.0 {
			stat_ser = coulomb_FG0_series(eta, x, &F_lam_min, &G_lam_min)
		} else if lam_min == 0.5 {
			/* This cannot happen. */
			F.val = F_lam_F
			F.err = 2.0 * gsl.Float64Eps * math.Abs(F.val)
			Fp.val = Fp_lam_F
			Fp.err = 2.0 * gsl.Float64Eps * math.Abs(Fp.val)
			G.val = G_lam_G
			G.err = 2.0 * gsl.Float64Eps * math.Abs(G.val)
			Gp.val = Gp_lam_G
			Gp.err = 2.0 * gsl.Float64Eps * math.Abs(Gp.val)
			*exp_F = 0.0
			*exp_G = 0.0
			return err.ERROR("error", err.ESANITY)
		} else {
			stat_ser = coulomb_FG_series(lam_min, eta, x, &F_lam_min, &G_lam_min)
		}

		/* Determine remaining quantities. */
		Fp_over_F_lam_min = Fp_lam_min_unnorm / F_lam_min_unnorm
		Gp_lam_min.val = Fp_over_F_lam_min*G_lam_min.val - 1.0/F_lam_min.val
		Gp_lam_min.err = math.Abs(Fp_over_F_lam_min) * G_lam_min.err
		Gp_lam_min.err += math.Abs(1.0/F_lam_min.val) * math.Abs(F_lam_min.err/F_lam_min.val)
		F_scale = F_lam_min.val / F_lam_min_unnorm

		/* Apply scale to the original F,F' values. */
		F_scale_frac_err = math.Abs(F_lam_min.err / F_lam_min.val)
		F_unnorm_frac_err = 2.0 * gsl.Float64Eps * float64(CF1_count+span+1)
		F_lam_F *= F_scale
		F_lam_F_err = math.Abs(F_lam_F) * (F_unnorm_frac_err + F_scale_frac_err)
		Fp_lam_F *= F_scale
		Fp_lam_F_err = math.Abs(Fp_lam_F) * (F_unnorm_frac_err + F_scale_frac_err)

		/* Recurse up to get the required G,G' values. */
		stat_Gr = coulomb_G_recur(lam_min, int(gsl.Max(N-k_lam_G, 0)), eta, x, G_lam_min.val, Gp_lam_min.val, &G_lam_G, &Gp_lam_G)

		F.val = F_lam_F
		F.err = F_lam_F_err
		F.err += 2.0 * gsl.Float64Eps * math.Abs(F_lam_F)

		Fp.val = Fp_lam_F
		Fp.err = Fp_lam_F_err
		Fp.err += 2.0 * gsl.Float64Eps * math.Abs(Fp_lam_F)

		Gerr_frac = math.Abs(G_lam_min.err/G_lam_min.val) + math.Abs(Gp_lam_min.err/Gp_lam_min.val)

		G.val = G_lam_G
		G.err = Gerr_frac * math.Abs(G_lam_G)
		G.err += 2.0 * float64(CF1_count+1) * gsl.Float64Eps * math.Abs(G.val)

		Gp.val = Gp_lam_G
		Gp.err = Gerr_frac * math.Abs(Gp.val)
		Gp.err += 2.0 * float64(CF1_count+1) * gsl.Float64Eps * math.Abs(Gp.val)

		*exp_F = 0.0
		*exp_G = 0.0

		return err.ErrorSelect(stat_ser, stat_CF1, stat_Fr, stat_Gr)
	} else if x < 2.0*eta {
		/* Use WKB approximation to obtain F and G at the two
		 * lambda values, and use the Wronskian and the
		 * continued fractions for F'/F to obtain F' and G'.
		 */
		var (
			F_lam_F, G_lam_F     Result
			F_lam_G, G_lam_G     Result
			exp_lam_F, exp_lam_G float64
			stat_lam_F           err.GSLError
			stat_lam_G           err.GSLError
			// stat_CF1_lam_F       error
			// stat_CF1_lam_G       error
			CF1_count       int
			Fp_over_F_lam_F float64
			Fp_over_F_lam_G float64
			F_sign_lam_F    float64
			F_sign_lam_G    float64
		)

		stat_lam_F = coulomb_jwkb(lam_F, eta, x, &F_lam_F, &G_lam_F, &exp_lam_F)
		if k_lam_G == 0 {
			stat_lam_G = stat_lam_F
			F_lam_G = F_lam_F
			G_lam_G = G_lam_F
			exp_lam_G = exp_lam_F
		} else {
			stat_lam_G = coulomb_jwkb(lam_G, eta, x, &F_lam_G, &G_lam_G, &exp_lam_G)
		}

		coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F, &CF1_count)
		if k_lam_G == 0 {
			// stat_CF1_lam_G = stat_CF1_lam_F
			F_sign_lam_G = F_sign_lam_F
			Fp_over_F_lam_G = Fp_over_F_lam_F
		} else {
			coulomb_CF1(lam_G, eta, x, &F_sign_lam_G, &Fp_over_F_lam_G, &CF1_count)
		}

		F.val = F_lam_F.val
		F.err = F_lam_F.err

		G.val = G_lam_G.val
		G.err = G_lam_G.err

		Fp.val = Fp_over_F_lam_F * F_lam_F.val
		Fp.err = math.Abs(Fp_over_F_lam_F) * F_lam_F.err
		Fp.err += 2.0 * gsl.Float64Eps * math.Abs(Fp.val)

		Gp.val = Fp_over_F_lam_G*G_lam_G.val - 1.0/F_lam_G.val
		Gp.err = math.Abs(Fp_over_F_lam_G) * G_lam_G.err
		Gp.err += math.Abs(1.0/F_lam_G.val) * math.Abs(F_lam_G.err/F_lam_G.val)

		*exp_F = exp_lam_F
		*exp_G = exp_lam_G

		if stat_lam_F.Status() == err.EOVRFLW {
			return err.ERROR("overflow", err.EOVRFLW)
		}

		if stat_lam_G.Status() == err.EOVRFLW {
			return err.ERROR("overflow", err.EOVRFLW)
		}

		return err.ErrorSelect(stat_lam_F, stat_lam_G)
	} else {
		/* x > 2 eta, so we know that we can find a lambda value such
		 * that x is above the turning point. We do this, evaluate
		 * using Steed's method at that oscillatory point, then
		 * use recursion on F and G to obtain the required values.
		 *
		 * lam_0   = a value of lambda such that x is below the turning point
		 * lam_min = minimum of lam_0 and the requested lam_G, since
		 *           we must go at least as low as lam_G
		 */
		var (
			SMALL                               = gsl.SqrtFloat64Eps
			C                                   = math.Sqrt(1.0 + 4.0*x*(x-2.0*eta))
			N                                   = int(math.Ceil(lam_F - C + 0.5))
			lam_0                               = lam_F - gsl.Max(N, 0)
			lam_min                             = math.Min(lam_0, lam_G)
			F_lam_F, Fp_lam_F                   float64
			G_lam_G, Gp_lam_G                   float64
			F_lam_min_unnorm, Fp_lam_min_unnorm float64
			F_lam_min                           float64 // , Fp_lam_min
			G_lam_min, Gp_lam_min               float64
			Fp_over_F_lam_F                     float64
			Fp_over_F_lam_min                   float64
			F_sign_lam_F, F_sign_lam_min        float64
			P_lam_min, Q_lam_min                float64
			alpha                               float64
			gamma                               float64
			F_scale                             float64
			CF1_count                           int
			CF2_count                           int
			stat_CF1                            = coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F, &CF1_count)
			stat_CF2                            err.GSLError
			stat_Fr                             err.GSLError
			stat_Gr                             err.GSLError
			F_recur_count                       int
			G_recur_count                       int
			err_amplify                         float64
		)

		F_lam_F = F_sign_lam_F * SMALL /* unnormalized */
		Fp_lam_F = Fp_over_F_lam_F * F_lam_F

		/* Backward recurrence to get F,Fp at lam_min */
		F_recur_count = int(gsl.Max(k_lam_G, N))
		stat_Fr = coulomb_F_recur(lam_min, F_recur_count, eta, x, F_lam_F, Fp_lam_F, &F_lam_min_unnorm, &Fp_lam_min_unnorm)
		Fp_over_F_lam_min = Fp_lam_min_unnorm / F_lam_min_unnorm

		/* Steed evaluation to complete evaluation of F,Fp,G,Gp at lam_min */
		stat_CF2 = coulomb_CF2(lam_min, eta, x, &P_lam_min, &Q_lam_min, &CF2_count)
		alpha = Fp_over_F_lam_min - P_lam_min
		gamma = alpha / Q_lam_min

		F_sign_lam_min = gsl.Sign(F_lam_min_unnorm)

		F_lam_min = F_sign_lam_min / math.Sqrt(alpha*alpha/Q_lam_min+Q_lam_min)
		// Fp_lam_min = Fp_over_F_lam_min * F_lam_min
		G_lam_min = gamma * F_lam_min
		Gp_lam_min = (P_lam_min*gamma - Q_lam_min) * F_lam_min

		/* Apply scale to values of F,Fp at lam_F (the top). */
		F_scale = F_lam_min / F_lam_min_unnorm
		F_lam_F *= F_scale
		Fp_lam_F *= F_scale

		/* Forward recurrence to get G,Gp at lam_G (the top). */
		G_recur_count = int(gsl.Max(N-k_lam_G, 0))
		stat_Gr = coulomb_G_recur(lam_min, G_recur_count, eta, x, G_lam_min, Gp_lam_min, &G_lam_G, &Gp_lam_G)

		err_amplify = float64(CF1_count + CF2_count + F_recur_count + G_recur_count + 1)

		F.val = F_lam_F
		F.err = 8.0 * err_amplify * gsl.Float64Eps * math.Abs(F.val)

		Fp.val = Fp_lam_F
		Fp.err = 8.0 * err_amplify * gsl.Float64Eps * math.Abs(Fp.val)

		G.val = G_lam_G
		G.err = 8.0 * err_amplify * gsl.Float64Eps * math.Abs(G.val)

		Gp.val = Gp_lam_G
		Gp.err = 8.0 * err_amplify * gsl.Float64Eps * math.Abs(Gp.val)

		*exp_F = 0.0
		*exp_G = 0.0

		return err.ErrorSelect(stat_CF1, stat_CF2, stat_Fr, stat_Gr)
	}
}

func Coulomb_wave_F_array(lam_min float64, kmax int, eta, x float64, fc_array []float64, F_exp *float64) err.GSLError {
	if x == 0.0 {
		var k int
		*F_exp = 0.0
		for k = 0; k <= kmax; k++ {
			fc_array[k] = 0.0
		}
		if lam_min == 0.0 {
			var f_0 Result
			cLeta(0.0, eta, &f_0)
			fc_array[0] = f_0.val
		}
		return nil
	} else {
		x_inv := 1.0 / x
		lam_max := lam_min + float64(kmax)
		var F, Fp, G, Gp Result
		var G_exp float64
		stat_FG := Coulomb_wave_FG_e(eta, x, lam_max, 0, &F, &Fp, &G, &Gp, F_exp, &G_exp)
		fcl := F.val
		fpl := Fp.val
		lam := lam_max
		var k int

		fc_array[kmax] = F.val

		for k = kmax - 1; k >= 0; k-- {
			el := eta / lam
			rl := math.Hypot(1.0, el)
			sl := el + lam*x_inv
			fc_lm1 := (fcl*sl + fpl) / rl
			fc_array[k] = fc_lm1
			fpl = fc_lm1*sl - fcl*rl
			fcl = fc_lm1
			lam -= 1.0
		}

		return stat_FG
	}
}

func Coulomb_wave_FG_array(lam_min float64, kmax int, eta, x float64, fc_array, gc_array []float64, F_exp, G_exp *float64) err.GSLError {
	x_inv := 1.0 / x
	lam_max := lam_min + float64(kmax)
	var F, Fp, G, Gp Result
	stat_FG := Coulomb_wave_FG_e(eta, x, lam_max, kmax, &F, &Fp, &G, &Gp, F_exp, G_exp)
	fcl := F.val
	fpl := Fp.val
	lam := lam_max
	var k int
	var gcl, gpl float64

	fc_array[kmax] = F.val

	for k = kmax - 1; k >= 0; k-- {
		el := eta / lam
		rl := math.Hypot(1.0, el)
		sl := el + lam*x_inv
		var fc_lm1 float64
		fc_lm1 = (fcl*sl + fpl) / rl
		fc_array[k] = fc_lm1
		fpl = fc_lm1*sl - fcl*rl
		fcl = fc_lm1
		lam -= 1.0
	}

	gcl = G.val
	gpl = Gp.val
	lam = lam_min + 1.0

	gc_array[0] = G.val

	for k = 1; k <= kmax; k++ {
		el := eta / lam
		rl := math.Hypot(1.0, el)
		sl := el + lam*x_inv
		gcl1 := (sl*gcl - gpl) / rl
		gc_array[k] = gcl1
		gpl = rl*gcl - sl*gcl1
		gcl = gcl1
		lam += 1.0
	}

	return stat_FG
}

func Coulomb_wave_FGp_array(lam_min float64, kmax int, eta, x float64, fc_array, fcp_array, gc_array, gcp_array []float64, F_exp, G_exp *float64) err.GSLError {
	var (
		x_inv        = 1.0 / x
		lam_max      = lam_min + float64(kmax)
		F, Fp, G, Gp Result
		stat_FG      = Coulomb_wave_FG_e(eta, x, lam_max, kmax, &F, &Fp, &G, &Gp, F_exp, G_exp)
		fcl          = F.val
		fpl          = Fp.val
		lam          = lam_max
		k            int
		gcl, gpl     float64
	)

	fc_array[kmax] = F.val
	fcp_array[kmax] = Fp.val

	for k = kmax - 1; k >= 0; k-- {
		el := eta / lam
		rl := math.Hypot(1.0, el)
		sl := el + lam*x_inv
		var fc_lm1 float64
		fc_lm1 = (fcl*sl + fpl) / rl
		fc_array[k] = fc_lm1
		fpl = fc_lm1*sl - fcl*rl
		fcp_array[k] = fpl
		fcl = fc_lm1
		lam -= 1.0
	}

	gcl = G.val
	gpl = Gp.val
	lam = lam_min + 1.0

	gc_array[0] = G.val
	gcp_array[0] = Gp.val

	for k = 1; k <= kmax; k++ {
		el := eta / lam
		rl := math.Hypot(1.0, el)
		sl := el + lam*x_inv
		gcl1 := (sl*gcl - gpl) / rl
		gc_array[k] = gcl1
		gpl = rl*gcl - sl*gcl1
		gcp_array[k] = gpl
		gcl = gcl1
		lam += 1.0
	}

	return stat_FG
}

func Coulomb_wave_sphF_array(lam_min float64, kmax int, eta, x float64, fc_array []float64, F_exp *float64) err.GSLError {
	if x < 0.0 || lam_min < -0.5 {
		return err.ERROR("domain error", err.EDOM)
	} else if x < 10.0/gsl.MaxFloat64 {
		var k int
		for k = 0; k <= kmax; k++ {
			fc_array[k] = 0.0
		}
		if lam_min == 0.0 {
			fc_array[0] = math.Sqrt(c0sq(eta))
		}
		*F_exp = 0.0
		if x == 0.0 {
			return nil
		} else {
			return err.ERROR("underflow", err.EUNDRFLW)
		}
	} else {
		var k int
		stat_F := Coulomb_wave_F_array(lam_min, kmax, eta, x, fc_array, F_exp)

		for k = 0; k <= kmax; k++ {
			fc_array[k] = fc_array[k] / x
		}

		return stat_F
	}
}
