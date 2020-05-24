/* specfunc/elljac.c
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

/* GJ: See [Thompson, Atlas for Computing Mathematical Functions] */

/* BJG 2005-07: New algorithm based on Algorithm 5 from Numerische
   Mathematik 7, 78-90 (1965) "Numerical Calculation of Elliptic
   Integrals and Elliptic Functions" R. Bulirsch.

   Minor tweak is to avoid division by zero when sin(x u_l) = 0 by
   computing reflected values sn(K-u) cn(K-u) dn(K-u) and using
   transformation from Abramowitz & Stegun table 16.8 column "K-u"*/
func Elljac_e(u, m float64, sn, cn, dn *float64) err.GSLError {
	if math.Abs(m) > 1.0 {
		*sn = 0.0
		*cn = 0.0
		*dn = 0.0
		return err.ERROR("|m| > 1.0", err.EDOM)
	} else if math.Abs(m) < 2.0*gsl.Float64Eps {
		*sn = math.Sin(u)
		*cn = math.Cos(u)
		*dn = 1.0
		return nil
	} else if math.Abs(m-1.0) < 2.0*gsl.Float64Eps {
		*sn = math.Tanh(u)
		*cn = 1.0 / math.Cosh(u)
		*dn = *cn
		return nil
	} else {
		var status err.GSLError
		N := 16
		var mu, nu, c, d [16]float64
		var sin_umu, cos_umu, t, r float64
		n := 0

		mu[0] = 1.0
		nu[0] = math.Sqrt(1.0 - m)

		for math.Abs(mu[n]-nu[n]) > 4.0*gsl.Float64Eps*math.Abs(mu[n]+nu[n]) {
			mu[n+1] = 0.5 * (mu[n] + nu[n])
			nu[n+1] = math.Sqrt(mu[n] * nu[n])
			n++
			if n >= N-1 {
				status = err.MaxIteration()
				break
			}
		}

		sin_umu = math.Sin(u * mu[n])
		cos_umu = math.Cos(u * mu[n])

		/* Since sin(u*mu(n)) can be zero we switch to computing sn(K-u),
		   cn(K-u), dn(K-u) when |sin| < |cos| */

		if math.Abs(sin_umu) < math.Abs(cos_umu) {
			t = sin_umu / cos_umu

			c[n] = mu[n] * t
			d[n] = 1.0

			for n > 0 {
				n--
				c[n] = d[n+1] * c[n+1]
				r = (c[n+1] * c[n+1]) / mu[n+1]
				d[n] = (r + nu[n]) / (r + mu[n])
			}

			*dn = math.Sqrt(1.0-m) / d[n]
			*cn = (*dn) * gsl.Sign(cos_umu) / Hypot(1.0, c[n])
			*sn = (*cn) * c[n] / math.Sqrt(1.0-m)
		} else {
			t = cos_umu / sin_umu

			c[n] = mu[n] * t
			d[n] = 1.0

			for n > 0 {
				n--
				c[n] = d[n+1] * c[n+1]
				r = (c[n+1] * c[n+1]) / mu[n+1]
				d[n] = (r + nu[n]) / (r + mu[n])
			}

			*dn = d[n]
			*sn = gsl.Sign(sin_umu) / Hypot(1.0, c[n])
			*cn = c[n] * (*sn)
		}

		return status
	}
}
