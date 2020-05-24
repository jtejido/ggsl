/* integration/qawo.c
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
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"math"
)

func Qawo(f gsl.Function, a, epsabs, epsrel float64, limit int, workspace *Workspace, wf *QAWOTable, result, abserr *float64) err.GSLError {
	var (
		area, errsum                                   float64
		res_ext, err_ext                               float64
		result0, abserr0, resabs0, resasc0             float64
		tolerance                                      float64
		ertest                                         float64
		error_over_large_intervals                     float64
		reseps, abseps, correc                         float64
		ktmin                                          int
		roundoff_type1, roundoff_type2, roundoff_type3 int
		error_type                                     int
		error_type2                                    bool
		iteration                                      int
		positive_integrand                             bool
		extrapolate                                    bool
		extall                                         bool
		disallow_extrapolation                         bool
		table                                          = NewExtrapolationTable()
		b                                              = a + wf.L
		abs_omega                                      = math.Abs(wf.omega)
	)

	/* Initialize results */
	workspace.Initialize(a, b)

	*result = 0
	*abserr = 0

	if limit > workspace.limit {
		return err.ERROR("iteration limit exceeds available workspace", err.EINVAL)
	}

	/* Test on accuracy */

	if epsabs <= 0 && (epsrel < 50*gsl.Float64Eps || epsrel < 0.5e-28) {
		return err.ERROR("tolerance cannot be achieved with given epsabs and epsrel", err.EBADTOL)
	}

	/* Perform the first integration */

	qc25f(f, a, b, wf, 0, &result0, &abserr0, &resabs0, &resasc0)

	workspace.SetInitialResult(result0, abserr0)

	tolerance = math.Max(epsabs, epsrel*math.Abs(result0))

	if abserr0 <= 100*gsl.Float64Eps*resabs0 && abserr0 > tolerance {
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

	/* Initialization */

	if 0.5*abs_omega*math.Abs(b-a) <= 2 {
		table.AppendTable(result0)
		extall = true
	}

	area = result0
	errsum = abserr0

	res_ext = result0
	err_ext = gsl.MaxFloat64

	positive_integrand = test_positivity(result0, resabs0)

	iteration = 1

	for iteration < limit {
		var (
			current_level           int
			a1, b1, a2, b2          float64
			a_i, b_i, r_i, e_i      float64
			area1, area2, area12    float64
			error1, error2, error12 float64
			resasc1, resasc2        float64
			resabs1, resabs2        float64
			last_e_i                float64
		)

		/* Bisect the subinterval with the largest error estimate */

		workspace.Retrieve(&a_i, &b_i, &r_i, &e_i)

		current_level = workspace.level[workspace.i] + 1

		if current_level >= wf.n {
			error_type = -1 /* exceeded limit of table */
			break
		}

		a1 = a_i
		b1 = 0.5 * (a_i + b_i)
		a2 = b1
		b2 = b_i

		iteration++

		qc25f(f, a1, b1, wf, current_level, &area1, &error1, &resabs1, &resasc1)
		qc25f(f, a2, b2, wf, current_level, &area2, &error2, &resabs2, &resasc2)

		area12 = area1 + area2
		error12 = error1 + error2
		last_e_i = e_i

		/* Improve previous approximations to the integral and test for
		   accuracy.

		   We write these expressions in the same way as the original
		   QUADPACK code so that the rounding errors are the same, which
		   makes testing easier. */

		errsum = errsum + error12 - e_i
		area = area + area12 - r_i

		tolerance = math.Max(epsabs, epsrel*math.Abs(area))

		if resasc1 != error1 && resasc2 != error2 {
			delta := r_i - area12

			if math.Abs(delta) <= 1.0e-5*math.Abs(area12) && error12 >= 0.99*e_i {
				if !extrapolate {
					roundoff_type1++
				} else {
					roundoff_type2++
				}
			}
			if iteration > 10 && error12 > e_i {
				roundoff_type3++
			}
		}

		/* Test for roundoff and eventually set error flag */

		if roundoff_type1+roundoff_type2 >= 10 || roundoff_type3 >= 20 {
			error_type = 2 /* round off error */
		}

		if roundoff_type2 >= 5 {
			error_type2 = true
		}

		/* set error flag in the case of bad integrand behaviour at
		   a point of the integration range */

		if subinterval_too_small(a1, a2, b2) {
			error_type = 4
		}

		/* append the newly-created intervals to the list */

		workspace.Update(a1, b1, area1, error1, a2, b2, area2, error2)

		if errsum <= tolerance {
			goto compute_result
		}

		if error_type != 0 {
			break
		}

		if iteration >= limit-1 {
			error_type = 1
			break
		}

		/* set up variables on first iteration */

		if iteration == 2 && extall {
			error_over_large_intervals = errsum
			ertest = tolerance
			table.AppendTable(area)
			continue
		}

		if disallow_extrapolation {
			continue
		}

		if extall {
			error_over_large_intervals += -last_e_i

			if current_level < workspace.maximum_level {
				error_over_large_intervals += error12
			}

			if extrapolate {
				goto label70
			}
		}

		if workspace.IsLargeInterval() {
			continue
		}

		if extall {
			extrapolate = true
			workspace.nrmax = 1
		} else {
			/* test whether the interval to be bisected next is the
			   smallest interval. */
			i := workspace.i
			width := workspace.blist[i] - workspace.alist[i]

			if 0.25*math.Abs(width)*abs_omega > 2 {
				continue
			}

			extall = true
			error_over_large_intervals = errsum
			ertest = tolerance
			continue
		}

	label70:
		if !error_type2 && error_over_large_intervals > ertest {
			if workspace.IsIncreaseNrmax() {
				continue
			}
		}

		/* Perform extrapolation */

		table.AppendTable(area)

		if table.n < 3 {
			workspace.ResetNrmax()
			extrapolate = false
			error_over_large_intervals = errsum
			continue
		}

		table.Qelg(&reseps, &abseps)

		ktmin++

		if ktmin > 5 && err_ext < 0.001*errsum {
			error_type = 5
		}

		if abseps < err_ext {
			ktmin = 0
			err_ext = abseps
			res_ext = reseps
			correc = error_over_large_intervals
			ertest = math.Max(epsabs, epsrel*math.Abs(reseps))
			if err_ext <= ertest {
				break
			}
		}

		/* Prepare bisection of the smallest interval. */

		if table.n == 1 {
			disallow_extrapolation = true
		}

		if error_type == 5 {
			break
		}

		/* work on interval with largest error */

		workspace.ResetNrmax()
		extrapolate = false
		error_over_large_intervals = errsum

	}

	*result = res_ext
	*abserr = err_ext

	if err_ext == gsl.MaxFloat64 {
		goto compute_result
	}

	if error_type != 0 || error_type2 {
		if error_type2 {
			err_ext += correc
		}

		if error_type == 0 {
			error_type = 3
		}

		if *result != 0 && area != 0 {
			if err_ext/math.Abs(res_ext) > errsum/math.Abs(area) {
				goto compute_result
			}
		} else if err_ext > errsum {
			goto compute_result
		} else if area == 0.0 {
			goto return_error
		}
	}

	/*  Test on divergence. */

	{
		max_area := math.Max(math.Abs(res_ext), math.Abs(area))

		if !positive_integrand && max_area < 0.01*resabs0 {
			goto return_error
		}
	}

	{
		ratio := res_ext / area

		if ratio < 0.01 || ratio > 100 || errsum > math.Abs(area) {
			error_type = 6
		}
	}

	goto return_error

compute_result:

	*result = workspace.SumResults()
	*abserr = errsum

return_error:

	if error_type > 2 {
		error_type--
	}

	if error_type == 0 {
		return nil
	} else if error_type == 1 {
		return err.ERROR("number of iterations was insufficient", err.EMAXITER)
	} else if error_type == 2 {
		return err.ERROR("cannot reach tolerance because of roundoff error", err.EROUND)
	} else if error_type == 3 {
		return err.ERROR("bad integrand behavior found in the integration interval", err.ESING)
	} else if error_type == 4 {
		return err.ERROR("roundoff error detected in extrapolation table", err.EROUND)
	} else if error_type == 5 {
		return err.ERROR("integral is divergent, or slowly convergent", err.EDIVERGE)
	} else if error_type == -1 {
		return err.ERROR("exceeded limit of trigonometric table", err.ETABLE)
	}

	return err.ERROR("could not integrate function", err.EFAILED)

}
