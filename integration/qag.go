/* integration/qag.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
package integration

import (
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"math"
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Qag(f gsl.Function, a, b, epsabs, epsrel float64, limit int, workspace *Workspace, result, abserr *float64, q QKFunction) err.GSLError {
	if q == nil {
		return err.ERROR("q does not specify a known integration rule", err.EINVAL)
	}

	var (
		area, errsum                       float64
		result0, abserr0, resabs0, resasc0 float64
		tolerance                          float64
		iteration                          int
		roundoff_type1                     int
		roundoff_type2                     int
		error_type                         int
		round_off                          float64
	)

	// Initialize results
	workspace.Initialize(a, b)

	*result = 0
	*abserr = 0

	if limit > workspace.limit {
		return err.ERROR("iteration limit exceeds available workspace", err.EINVAL)
	}

	if epsabs <= 0 && (epsrel < 50*gsl.Float64Eps || epsrel < 0.5e-28) {
		return err.ERROR("tolerance cannot be acheived with given epsabs and epsrel", err.EBADTOL)
	}

	// perform the first integration
	q(f, a, b, &result0, &abserr0, &resabs0, &resasc0)
	workspace.SetInitialResult(result0, abserr0)

	// Test on accuracy
	tolerance = gsl.Max(epsabs, epsrel*math.Abs(result0))
	round_off = 50 * gsl.Float64Eps * resabs0
	if abserr0 <= round_off && abserr0 > tolerance {
		*result = result0
		*abserr = abserr0
		return err.ERROR("cannot reach tolerance because of roundoff error on first attempt", err.EROUND)
	} else if (abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0.0 {
		*result = result0
		*abserr = abserr0
		return nil
	} else if limit == 1 {
		*result = result0
		*abserr = abserr0
		return err.ERROR("a maximum of one iteration was insufficient", err.EMAXITER)
	}

	area = result0
	errsum = abserr0

	iteration = 1
	for iteration < limit && errsum > tolerance && error_type == 0 {
		var (
			a1, b1, a2, b2          float64
			a_i, b_i, r_i, e_i      float64
			area1, area2, area12    float64
			error1, error2, error12 float64
			resasc1, resasc2        float64
			resabs1, resabs2        float64
		)

		workspace.Retrieve(&a_i, &b_i, &r_i, &e_i)
		a1 = a_i
		b1 = 0.5 * (a_i + b_i)
		a2 = b1
		b2 = b_i

		q(f, a1, b1, &area1, &error1, &resabs1, &resasc1)
		q(f, a2, b2, &area2, &error2, &resabs2, &resasc2)

		area12 = area1 + area2
		error12 = error1 + error2

		errsum += (error12 - e_i)
		area += area12 - r_i

		if resasc1 != error1 && resasc2 != error2 {
			delta := r_i - area12

			if math.Abs(delta) <= 1.0e-5*math.Abs(area12) && error12 >= 0.99*e_i {
				roundoff_type1++
			}
			if iteration >= 10 && error12 > e_i {
				roundoff_type2++
			}
		}

		tolerance = math.Max(epsabs, epsrel*math.Abs(area))

		if errsum > tolerance {
			if roundoff_type1 >= 6 || roundoff_type2 >= 20 {
				error_type = 2 /* round off error */
			}

			/* set error flag in the case of bad integrand behaviour at
			   a point of the integration range */

			if subinterval_too_small(a1, a2, b2) {
				error_type = 3
			}
		}

		workspace.Update(a1, b1, area1, error1, a2, b2, area2, error2)
		workspace.Retrieve(&a_i, &b_i, &r_i, &e_i)

		iteration++
	}

	*result = workspace.SumResults()
	*abserr = errsum

	if errsum <= tolerance {
		return nil
	} else if error_type == 2 {
		return err.ERROR("roundoff error prevents tolerance from being achieved", err.EROUND)
	} else if error_type == 3 {
		return err.ERROR("bad integrand behavior found in the integration interval", err.ESING)
	} else if iteration == limit {
		return err.ERROR("maximum number of subdivisions reached", err.EMAXITER)
	}

	return err.ERROR("could not integrate function", err.EFAILED)
}
