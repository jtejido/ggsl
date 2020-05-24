/* integration/qags.c
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

/*-*-*-*-*-*-*-*-*-* Private Structs *-*-*-*-*-*-*-*-*-*-*/
type i_transform struct {
	f gsl.Function
}

func (it *i_transform) Evaluate(t float64) float64 {
	x := (1 - t) / t
	y := it.f.Evaluate(x) + it.f.Evaluate(-x)
	return (y / t) / t
}

type il_transform struct {
	f gsl.Function
	b float64
}

func (ilt *il_transform) Evaluate(t float64) float64 {
	x := ilt.b - (1-t)/t
	y := ilt.f.Evaluate(x)
	return (y / t) / t
}

type iu_transform struct {
	f gsl.Function
	a float64
}

func (iut *iu_transform) Evaluate(t float64) float64 {
	x := iut.a + (1-t)/t
	y := iut.f.Evaluate(x)
	return (y / t) / t
}

func qags(f gsl.Function, a, b, epsabs, epsrel float64, limit int, workspace *Workspace, result, abserr *float64, q QKFunction) err.GSLError {
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
		disallow_extrapolation                         bool
		table                                          = NewExtrapolationTable()
	)

	/* Initialize results */
	workspace.Initialize(a, b)

	*result = 0
	*abserr = 0

	if limit > workspace.limit {
		err.ERROR("iteration limit exceeds available workspace", err.EINVAL)
	}

	/* Test on accuracy */
	if epsabs <= 0 && (epsrel < 50*gsl.Float64Eps || epsrel < 0.5e-28) {
		return err.ERROR("tolerance cannot be achieved with given epsabs and epsrel", err.EBADTOL)
	}

	/* Perform the first integration */
	q(f, a, b, &result0, &abserr0, &resabs0, &resasc0)
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
	table.AppendTable(result0)
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

		a1 = a_i
		b1 = 0.5 * (a_i + b_i)
		a2 = b1
		b2 = b_i

		iteration++

		q(f, a1, b1, &area1, &error1, &resabs1, &resasc1)
		q(f, a2, b2, &area2, &error2, &resabs2, &resasc2)

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

		if iteration == 2 {
			error_over_large_intervals = errsum /* set up variables on first iteration */
			ertest = tolerance
			table.AppendTable(area)
			continue
		}

		if disallow_extrapolation {
			continue
		}

		error_over_large_intervals += -last_e_i

		if current_level < workspace.maximum_level {
			error_over_large_intervals += error12
		}

		if !extrapolate {
			/* test whether the interval to be bisected next is the
			   smallest interval. */

			if workspace.IsLargeInterval() {
				continue
			}

			extrapolate = true
			workspace.nrmax = 1
		}

		if !error_type2 && error_over_large_intervals > ertest {
			if workspace.IsIncreaseNrmax() {
				continue
			}
		}

		/* Perform extrapolation */

		table.AppendTable(area)

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

		if res_ext != 0.0 && area != 0.0 {
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

		if ratio < 0.01 || ratio > 100.0 || errsum > math.Abs(area) {
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
		return err.ERROR("roundoff error detected in the extrapolation table", err.EROUND)
	} else if error_type == 5 {
		return err.ERROR("integral is divergent, or slowly convergent", err.EDIVERGE)
	}

	return err.ERROR("could not integrate function", err.EFAILED)

}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Qags(f gsl.Function, a, b, epsabs, epsrel float64, limit int, workspace *Workspace, result, abserr *float64) err.GSLError {
	return qags(f, a, b, epsabs, epsrel, limit, workspace, result, abserr, Qk21)
}

func Qagi(f gsl.Function, epsabs, epsrel float64, limit int, workspace *Workspace, result, abserr *float64) err.GSLError {
	f_transform := &i_transform{f}
	return qags(f_transform, 0.0, 1.0, epsabs, epsrel, limit, workspace, result, abserr, Qk15)
}

func Qagil(f gsl.Function, b, epsabs, epsrel float64, limit int, workspace *Workspace, result, abserr *float64) err.GSLError {
	f_transform := &il_transform{f, b}
	return qags(f_transform, 0.0, 1.0, epsabs, epsrel, limit, workspace, result, abserr, Qk15)
}

func Qagiu(f gsl.Function, a, epsabs, epsrel float64, limit int, workspace *Workspace, result, abserr *float64) err.GSLError {
	f_transform := &iu_transform{f, a}
	return qags(f_transform, 0.0, 1.0, epsabs, epsrel, limit, workspace, result, abserr, Qk15)
}
