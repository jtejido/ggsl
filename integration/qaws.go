/* integration/qaws.c
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

func Qaws(f gsl.Function, a, b float64, t *QAWSTable, epsabs, epsrel float64, limit int, workspace *Workspace, result, abserr *float64) err.GSLError {
	var (
		area, errsum                               float64
		result0, abserr0                           float64
		tolerance                                  float64
		iteration                                  int
		roundoff_type1, roundoff_type2, error_type int
	)

	/* Initialize results */

	workspace.Initialize(a, b)

	*result = 0
	*abserr = 0

	if limit > workspace.limit {
		err.ERROR("iteration limit exceeds available workspace", err.EINVAL)
	}

	if b <= a {
		err.ERROR("limits must form an ascending sequence, a < b", err.EINVAL)
	}

	if epsabs <= 0 && (epsrel < 50*gsl.Float64Eps || epsrel < 0.5e-28) {
		err.ERROR("tolerance cannot be achieved with given epsabs and epsrel", err.EBADTOL)
	}

	/* perform the first integration */

	{
		var (
			area1, area2                 float64
			error1, error2               float64
			err_reliable1, err_reliable2 int
			a1                           float64 = a
			b1                           float64 = 0.5 * (a + b)
			a2                           float64 = b1
			b2                           float64 = b
		)

		qc25s(f, a, b, a1, b1, t, &area1, &error1, &err_reliable1)
		qc25s(f, a, b, a2, b2, t, &area2, &error2, &err_reliable2)

		if error1 > error2 {
			workspace.AppendInterval(a1, b1, area1, error1)
			workspace.AppendInterval(a2, b2, area2, error2)
		} else {
			workspace.AppendInterval(a2, b2, area2, error2)
			workspace.AppendInterval(a1, b1, area1, error1)
		}

		result0 = area1 + area2
		abserr0 = error1 + error2
	}

	/* Test on accuracy */

	tolerance = math.Max(epsabs, epsrel*math.Abs(result0))

	/* Test on accuracy, use 0.01 relative error as an extra safety
	   margin on the first iteration (ignored for subsequent iterations) */

	if abserr0 < tolerance && abserr0 < 0.01*math.Abs(result0) {
		*result = result0
		*abserr = abserr0

		return nil
	} else if limit == 1 {
		*result = result0
		*abserr = abserr0

		err.ERROR("a maximum of one iteration was insufficient", err.EMAXITER)
	}

	area = result0
	errsum = abserr0

	iteration = 2

	for iteration < limit && error_type == 0 && errsum > tolerance {
		var (
			a1, b1, a2, b2               float64
			a_i, b_i, r_i, e_i           float64
			area1, area2, area12         float64
			error1, error2, error12      float64
			err_reliable1, err_reliable2 int
		)
		/* Bisect the subinterval with the largest error estimate */

		workspace.Retrieve(&a_i, &b_i, &r_i, &e_i)

		a1 = a_i
		b1 = 0.5 * (a_i + b_i)
		a2 = b1
		b2 = b_i

		qc25s(f, a, b, a1, b1, t, &area1, &error1, &err_reliable1)
		qc25s(f, a, b, a2, b2, t, &area2, &error2, &err_reliable2)

		area12 = area1 + area2
		error12 = error1 + error2

		errsum += (error12 - e_i)
		area += area12 - r_i

		if err_reliable1 != 0 && err_reliable2 != 0 {
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
		err.ERROR("roundoff error prevents tolerance from being achieved", err.EROUND)
	} else if error_type == 3 {
		err.ERROR("bad integrand behavior found in the integration interval", err.ESING)
	} else if iteration == limit {
		err.ERROR("maximum number of subdivisions reached", err.EMAXITER)
	}

	return err.ERROR("could not integrate function", err.EFAILED)

}
