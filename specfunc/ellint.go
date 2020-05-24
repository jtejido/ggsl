/* specfunc/ellint.c
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
	"github.com/jtejido/ggsl/sys"
	"math"
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/* based on Carlson's algorithms:
   [B. C. Carlson Numer. Math. 33, 1 (1979)]

   see also:
   [B.C. Carlson, Special Functions of Applied Mathematics (1977)]
*/

/* According to Carlson's algorithm, the errtol parameter
   typically effects the relative error in the following way:

   relative error < 16 errtol^6 / (1 - 2 errtol)

     errtol     precision
     ------     ----------
     0.001       1.0e-17
     0.003       2.0e-14
     0.01        2.0e-11
     0.03        2.0e-8
     0.1         2.0e-5
*/
func Ellint_RC_e(x, y float64, mode gsl.MODE_T, result *Result) err.GSLError {

	lolim := 5.0 * gsl.MinFloat64
	uplim := 0.2 * gsl.MaxFloat64
	errtol := 0.03
	if mode == gsl.PREC_DOUBLE {
		errtol = 0.001
	}

	prec := sys.PrecEps[mode]
	nmax := 10000

	if x < 0.0 || y < 0.0 || x+y < lolim {
		return DomainError(result)
	} else if gsl.Max(x, y) < uplim {
		c1 := 1.0 / 7.0
		c2 := 9.0 / 22.0
		xn := x
		yn := y
		var mu, sn, lamda, s float64
		n := 0
		for {
			mu = (xn + yn + yn) / 3.0
			sn = (yn+mu)/mu - 2.0
			if math.Abs(sn) < errtol {
				break
			}
			lamda = 2.0*math.Sqrt(xn)*math.Sqrt(yn) + yn
			xn = (xn + lamda) * 0.25
			yn = (yn + lamda) * 0.25
			n++
			if n == nmax {
				return MaxIterError(result)
			}
		}
		s = sn * sn * (0.3 + sn*(c1+sn*(0.375+sn*c2)))
		result.val = (1.0 + s) / math.Sqrt(mu)
		result.err = prec * math.Abs(result.val)
		return nil
	} else {
		return DomainError(result)
	}
}

func Ellint_RD_e(x, y, z float64, mode gsl.MODE_T, result *Result) err.GSLError {
	errtol := 0.03
	if mode == gsl.PREC_DOUBLE {
		errtol = 0.001
	}
	prec := sys.PrecEps[mode]
	lolim := 2.0 / math.Pow(gsl.MaxFloat64, 2.0/3.0)
	uplim := math.Pow(0.1*errtol/gsl.MinFloat64, 2.0/3.0)
	nmax := 10000

	if gsl.Min(x, y) < 0.0 || gsl.Min(x+y, z) < lolim {
		return DomainError(result)
	} else if gsl.Max(x, y, z) < uplim {
		c1 := 3.0 / 14.0
		c2 := 1.0 / 6.0
		c3 := 9.0 / 22.0
		c4 := 3.0 / 26.0
		xn := x
		yn := y
		zn := z
		sigma := 0.0
		power4 := 1.0
		var ea, eb, ec, ed, ef, s1, s2 float64
		var mu, xndev, yndev, zndev float64
		n := 0
		for {
			var xnroot, ynroot, znroot, lamda float64
			var epslon float64
			mu = (xn + yn + 3.0*zn) * 0.2
			xndev = (mu - xn) / mu
			yndev = (mu - yn) / mu
			zndev = (mu - zn) / mu
			epslon = gsl.Max(math.Abs(xndev), math.Abs(yndev), math.Abs(zndev))
			if epslon < errtol {
				break
			}
			xnroot = math.Sqrt(xn)
			ynroot = math.Sqrt(yn)
			znroot = math.Sqrt(zn)
			lamda = xnroot*(ynroot+znroot) + ynroot*znroot
			sigma += power4 / (znroot * (zn + lamda))
			power4 *= 0.25
			xn = (xn + lamda) * 0.25
			yn = (yn + lamda) * 0.25
			zn = (zn + lamda) * 0.25
			n++
			if n == nmax {
				return MaxIterError(result)
			}
		}
		ea = xndev * yndev
		eb = zndev * zndev
		ec = ea - eb
		ed = ea - 6.0*eb
		ef = ed + ec + ec
		s1 = ed * (-c1 + 0.25*c3*ed - 1.5*c4*zndev*ef)
		s2 = zndev * (c2*ef + zndev*(-c3*ec+zndev*c4*ea))
		result.val = 3.0*sigma + power4*(1.0+s1+s2)/(mu*math.Sqrt(mu))
		result.err = prec * math.Abs(result.val)
		return nil
	} else {
		return DomainError(result)
	}
}

func Ellint_RF_e(x, y, z float64, mode gsl.MODE_T, result *Result) err.GSLError {
	lolim := 5.0 * gsl.MinFloat64
	uplim := 0.2 * gsl.MaxFloat64
	errtol := 0.03
	if mode == gsl.PREC_DOUBLE {
		errtol = 0.001
	}
	prec := sys.PrecEps[mode]
	nmax := 10000

	if x < 0.0 || y < 0.0 || z < 0.0 {
		return DomainError(result)
	} else if x+y < lolim || x+z < lolim || y+z < lolim {
		return DomainError(result)
	} else if gsl.Max(x, y, z) < uplim {
		c1 := 1.0 / 24.0
		c2 := 3.0 / 44.0
		c3 := 1.0 / 14.0
		xn := x
		yn := y
		zn := z
		var mu, xndev, yndev, zndev, e2, e3, s float64
		n := 0
		for {
			var epslon, lamda float64
			var xnroot, ynroot, znroot float64
			mu = (xn + yn + zn) / 3.0
			xndev = 2.0 - (mu+xn)/mu
			yndev = 2.0 - (mu+yn)/mu
			zndev = 2.0 - (mu+zn)/mu
			epslon = gsl.Max(math.Abs(xndev), math.Abs(yndev), math.Abs(zndev))
			if epslon < errtol {
				break
			}
			xnroot = math.Sqrt(xn)
			ynroot = math.Sqrt(yn)
			znroot = math.Sqrt(zn)
			lamda = xnroot*(ynroot+znroot) + ynroot*znroot
			xn = (xn + lamda) * 0.25
			yn = (yn + lamda) * 0.25
			zn = (zn + lamda) * 0.25
			n++
			if n == nmax {
				return MaxIterError(result)
			}
		}
		e2 = xndev*yndev - zndev*zndev
		e3 = xndev * yndev * zndev
		s = 1.0 + (c1*e2-0.1-c2*e3)*e2 + c3*e3
		result.val = s / math.Sqrt(mu)
		result.err = prec * math.Abs(result.val)
		return nil
	} else {
		return DomainError(result)
	}
}

func Ellint_RJ_e(x, y, z, p float64, mode gsl.MODE_T, result *Result) err.GSLError {
	errtol := 0.03
	if mode == gsl.PREC_DOUBLE {
		errtol = 0.001
	}
	prec := sys.PrecEps[mode]
	lolim := math.Pow(5.0*gsl.MinFloat64, 1.0/3.0)
	uplim := 0.3 * math.Pow(0.2*gsl.MaxFloat64, 1.0/3.0)
	nmax := 10000

	if x < 0.0 || y < 0.0 || z < 0.0 {
		return DomainError(result)
	} else if x+y < lolim || x+z < lolim || y+z < lolim || p < lolim {
		return DomainError(result)
	} else if gsl.Max(x, y, z, p) < uplim {
		c1 := 3.0 / 14.0
		c2 := 1.0 / 3.0
		c3 := 3.0 / 22.0
		c4 := 3.0 / 26.0
		xn := x
		yn := y
		zn := z
		pn := p
		sigma := 0.0
		power4 := 1.0
		var mu, xndev, yndev, zndev, pndev float64
		var ea, eb, ec, e2, e3, s1, s2, s3 float64
		n := 0
		for {
			var xnroot, ynroot, znroot float64
			var lamda, alfa, beta float64
			var epslon float64
			var rcresult Result
			var rcstatus err.GSLError
			mu = (xn + yn + zn + pn + pn) * 0.2
			xndev = (mu - xn) / mu
			yndev = (mu - yn) / mu
			zndev = (mu - zn) / mu
			pndev = (mu - pn) / mu
			epslon = gsl.Max(math.Abs(xndev), math.Abs(yndev), math.Abs(zndev), math.Abs(pndev))
			if epslon < errtol {
				break
			}
			xnroot = math.Sqrt(xn)
			ynroot = math.Sqrt(yn)
			znroot = math.Sqrt(zn)
			lamda = xnroot*(ynroot+znroot) + ynroot*znroot
			alfa = pn*(xnroot+ynroot+znroot) + xnroot*ynroot*znroot
			alfa = alfa * alfa
			beta = pn * (pn + lamda) * (pn + lamda)
			rcstatus = Ellint_RC_e(alfa, beta, mode, &rcresult)
			if rcstatus != nil {
				result.val = 0.0
				result.err = 0.0
				return rcstatus
			}
			sigma += power4 * rcresult.val
			power4 *= 0.25
			xn = (xn + lamda) * 0.25
			yn = (yn + lamda) * 0.25
			zn = (zn + lamda) * 0.25
			pn = (pn + lamda) * 0.25
			n++
			if n == nmax {
				return MaxIterError(result)
			}
		}

		ea = xndev*(yndev+zndev) + yndev*zndev
		eb = xndev * yndev * zndev
		ec = pndev * pndev
		e2 = ea - 3.0*ec
		e3 = eb + 2.0*pndev*(ea-ec)
		s1 = 1.0 + e2*(-c1+0.75*c3*e2-1.5*c4*e3)
		s2 = eb * (0.5*c2 + pndev*(-c3-c3+pndev*c4))
		s3 = pndev*ea*(c2-pndev*c3) - c2*pndev*ec
		result.val = 3.0*sigma + power4*(s1+s2+s3)/(mu*math.Sqrt(mu))
		result.err = prec * math.Abs(result.val)
		return nil
	} else {
		return DomainError(result)
	}
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.1)] */
func Ellint_F_e(phi, k float64, mode gsl.MODE_T, result *Result) err.GSLError {
	/* Angular reduction to -pi/2 < phi < pi/2 (we should really use an
	   exact reduction but this will have to do for now) BJG */
	nc := math.Floor(phi/gsl.Pi + 0.5)
	phi_red := phi - nc*gsl.Pi
	phi = phi_red

	{
		sin_phi := math.Sin(phi)
		sin2_phi := sin_phi * sin_phi
		x := 1.0 - sin2_phi
		y := 1.0 - k*k*sin2_phi
		var rf Result
		status := Ellint_RF_e(x, y, 1.0, mode, &rf)
		result.val = sin_phi * rf.val
		result.err = gsl.Float64Eps*math.Abs(result.val) + math.Abs(sin_phi*rf.err)
		if nc == 0 {
			return status
		} else {
			var rk Result /* add extra terms from periodicity */
			rkstatus := Ellint_Kcomp_e(k, mode, &rk)
			result.val += 2 * nc * rk.val
			result.err += 2 * math.Abs(nc) * rk.err
			return err.ErrorSelect(status, rkstatus)
		}
	}
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.2)] */
func Ellint_E_e(phi, k float64, mode gsl.MODE_T, result *Result) err.GSLError {
	/* Angular reduction to -pi/2 < phi < pi/2 (we should really use an
	   exact reduction but this will have to do for now) BJG */

	nc := math.Floor(phi/gsl.Pi + 0.5)
	phi_red := phi - nc*gsl.Pi
	phi = phi_red

	{
		sin_phi := math.Sin(phi)
		sin2_phi := sin_phi * sin_phi
		x := 1.0 - sin2_phi
		y := 1.0 - k*k*sin2_phi

		if x < gsl.Float64Eps {
			var re Result
			status := Ellint_Ecomp_e(k, mode, &re)
			/* could use A&S 17.4.14 to improve the value below */
			result.val = 2*nc*re.val + gsl.Sign(sin_phi)*re.val
			result.err = 2*math.Abs(nc)*re.err + re.err
			return status
		} else {
			var rf, rd Result
			sin3_phi := sin2_phi * sin_phi
			rfstatus := Ellint_RF_e(x, y, 1.0, mode, &rf)
			rdstatus := Ellint_RD_e(x, y, 1.0, mode, &rd)
			result.val = sin_phi*rf.val - k*k/3.0*sin3_phi*rd.val
			result.err = gsl.Float64Eps * math.Abs(sin_phi*rf.val)
			result.err += math.Abs(sin_phi * rf.err)
			result.err += k * k / 3.0 * gsl.Float64Eps * math.Abs(sin3_phi*rd.val)
			result.err += k * k / 3.0 * math.Abs(sin3_phi*rd.err)
			if nc == 0 {
				return err.ErrorSelect(rfstatus, rdstatus)
			} else {
				var re Result /* add extra terms from periodicity */
				restatus := Ellint_Ecomp_e(k, mode, &re)
				result.val += 2 * nc * re.val
				result.err += 2 * math.Abs(nc) * re.err
				return err.ErrorSelect(rfstatus, rdstatus, restatus)
			}
		}
	}
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.3)] */
func Ellint_P_e(phi, k, n float64, mode gsl.MODE_T, result *Result) err.GSLError {
	/* Angular reduction to -pi/2 < phi < pi/2 (we should really use an
	   exact reduction but this will have to do for now) BJG */

	nc := math.Floor(phi/gsl.Pi + 0.5)
	phi_red := phi - nc*gsl.Pi
	phi = phi_red

	/* FIXME: need to handle the case of small x, as for E,F */

	{
		sin_phi := math.Sin(phi)
		sin2_phi := sin_phi * sin_phi
		sin3_phi := sin2_phi * sin_phi
		x := 1.0 - sin2_phi
		y := 1.0 - k*k*sin2_phi
		var rf, rj Result
		rfstatus := Ellint_RF_e(x, y, 1.0, mode, &rf)
		rjstatus := Ellint_RJ_e(x, y, 1.0, 1.0+n*sin2_phi, mode, &rj)
		result.val = sin_phi*rf.val - n/3.0*sin3_phi*rj.val
		result.err = gsl.Float64Eps * math.Abs(sin_phi*rf.val)
		result.err += math.Abs(sin_phi * rf.err)
		result.err += n / 3.0 * gsl.Float64Eps * math.Abs(sin3_phi*rj.val)
		result.err += n / 3.0 * math.Abs(sin3_phi*rj.err)
		if nc == 0 {
			return err.ErrorSelect(rfstatus, rjstatus)
		} else {
			var rp Result /* add extra terms from periodicity */
			rpstatus := Ellint_Pcomp_e(k, n, mode, &rp)
			result.val += 2 * nc * rp.val
			result.err += 2 * math.Abs(nc) * rp.err
			return err.ErrorSelect(rfstatus, rjstatus, rpstatus)
		}
	}
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.4)] */
func Ellint_D_e(phi, k float64, mode gsl.MODE_T, result *Result) err.GSLError {
	/* Angular reduction to -pi/2 < phi < pi/2 (we should really use an
	   exact reduction but this will have to do for now) BJG */

	nc := math.Floor(phi/gsl.Pi + 0.5)
	phi_red := phi - nc*gsl.Pi
	phi = phi_red

	/* FIXME: need to handle the case of small x, as for E,F */
	{
		sin_phi := math.Sin(phi)
		sin2_phi := sin_phi * sin_phi
		sin3_phi := sin2_phi * sin_phi
		x := 1.0 - sin2_phi
		y := 1.0 - k*k*sin2_phi
		var rd Result
		status := Ellint_RD_e(x, y, 1.0, mode, &rd)
		result.val = sin3_phi / 3.0 * rd.val
		result.err = gsl.Float64Eps*math.Abs(result.val) + math.Abs(sin3_phi/3.0*rd.err)
		if nc == 0 {
			return status
		} else {
			var rd Result /* add extra terms from periodicity */
			rdstatus := Ellint_Dcomp_e(k, mode, &rd)
			result.val += 2 * nc * rd.val
			result.err += 2 * math.Abs(nc) * rd.err
			return err.ErrorSelect(status, rdstatus)
		}
	}
}

func Ellint_Dcomp_e(k float64, mode gsl.MODE_T, result *Result) err.GSLError {
	if k*k >= 1.0 {
		return DomainError(result)
	} else {
		y := 1.0 - k*k /* FIXME: still need to handle k~=~1 */
		var rd Result
		status := Ellint_RD_e(0.0, y, 1.0, mode, &rd)
		result.val = (1.0 / 3.0) * rd.val
		result.err = gsl.Float64Eps*math.Abs(result.val) + math.Abs(1.0/3.0*rd.err)
		return status
	}
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.5)] */
func Ellint_Kcomp_e(k float64, mode gsl.MODE_T, result *Result) err.GSLError {
	if k*k >= 1.0 {
		return DomainError(result)
	} else if k*k >= 1.0-gsl.SqrtFloat64Eps {
		/* [Abramowitz+Stegun, 17.3.34] */
		y := 1.0 - k*k
		a := []float64{1.38629436112, 0.09666344259, 0.03590092383}
		b := []float64{0.5, 0.12498593597, 0.06880248576}
		ta := a[0] + y*(a[1]+y*a[2])
		tb := -math.Log(y) * (b[0] + y*(b[1]+y*b[2]))
		result.val = ta + tb
		result.err = 2.0 * gsl.Float64Eps * (math.Abs(result.val) + math.Abs(k/y))
		return nil
	} else {
		/* This was previously computed as,

		     return gsl_sf_ellint_RF_e(0.0, 1.0 - k*k, 1.0, mode, result);

		   but this underestimated the total error for small k, since the
		   argument y=1-k^2 is not exact (there is an absolute error of
		   gsl.Float64Eps near y=0 due to cancellation in the subtraction).
		   Taking the singular behavior of -log(y) above gives an error
		   of 0.5*epsilon/y near y=0. (BJG) */

		y := 1.0 - k*k
		status := Ellint_RF_e(0.0, y, 1.0, mode, result)
		result.err += 0.5 * gsl.Float64Eps / y
		return status
	}
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.6)] */
func Ellint_Ecomp_e(k float64, mode gsl.MODE_T, result *Result) err.GSLError {
	if k*k >= 1.0 {
		return DomainError(result)
	} else if k*k >= 1.0-gsl.SqrtFloat64Eps {
		/* [Abramowitz+Stegun, 17.3.36] */
		y := 1.0 - k*k
		a := []float64{0.44325141463, 0.06260601220, 0.04757383546}
		b := []float64{0.24998368310, 0.09200180037, 0.04069697526}
		ta := 1.0 + y*(a[0]+y*(a[1]+a[2]*y))
		tb := -y * math.Log(y) * (b[0] + y*(b[1]+b[2]*y))
		result.val = ta + tb
		result.err = 2.0 * gsl.Float64Eps * result.val
		return nil
	} else {
		var rf, rd Result
		y := 1.0 - k*k
		rfstatus := Ellint_RF_e(0.0, y, 1.0, mode, &rf)
		rdstatus := Ellint_RD_e(0.0, y, 1.0, mode, &rd)
		result.val = rf.val - k*k/3.0*rd.val
		result.err = rf.err + k*k/3.0*rd.err
		return err.ErrorSelect(rfstatus, rdstatus)
	}
}

/* [Carlson, Numer. Math. 33 (1979) 1, (4.3) phi=pi/2] */
func Ellint_Pcomp_e(k, n float64, mode gsl.MODE_T, result *Result) err.GSLError {
	if k*k >= 1.0 {
		return DomainError(result)
	} else { /* FIXME: need to handle k ~=~ 1  cancellations */
		var rf, rj Result
		y := 1.0 - k*k
		rfstatus := Ellint_RF_e(0.0, y, 1.0, mode, &rf)
		rjstatus := Ellint_RJ_e(0.0, y, 1.0, 1.0+n, mode, &rj)
		result.val = rf.val - (n/3.0)*rj.val
		result.err = rf.err + math.Abs(n/3.0)*rj.err
		return err.ErrorSelect(rfstatus, rjstatus)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Ellint_Kcomp(k float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_Kcomp_e(k, mode, result)
	return EvalResult(result, status)
}

func Ellint_Ecomp(k float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_Ecomp_e(k, mode, result)
	return EvalResult(result, status)
}

func Ellint_Pcomp(k, n float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_Pcomp_e(k, n, mode, result)
	return EvalResult(result, status)
}

func Ellint_Dcomp(k float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_Dcomp_e(k, mode, result)
	return EvalResult(result, status)
}

func Ellint_F(phi, k float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_F_e(phi, k, mode, result)
	return EvalResult(result, status)
}

func Ellint_E(phi, k float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_E_e(phi, k, mode, result)
	return EvalResult(result, status)
}

func Ellint_P(phi, k, n float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_P_e(phi, k, n, mode, result)
	return EvalResult(result, status)
}

func Ellint_D(phi, k float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_D_e(phi, k, mode, result)
	return EvalResult(result, status)
}

func Ellint_RC(x, y float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_RC_e(x, y, mode, result)
	return EvalResult(result, status)
}

func Ellint_RD(x, y, z float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_RD_e(x, y, z, mode, result)
	return EvalResult(result, status)
}

func Ellint_RF(x, y, z float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_RF_e(x, y, z, mode, result)
	return EvalResult(result, status)
}

func Ellint_RJ(x, y, z, p float64, mode gsl.MODE_T) float64 {
	result := new(Result)
	status := Ellint_RJ_e(x, y, z, p, mode, result)
	return EvalResult(result, status)
}
