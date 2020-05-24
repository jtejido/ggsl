/* integration/qawf.c
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

func Qawf(f gsl.Function, a, epsabs float64, limit int, workspace, cycle_workspace *Workspace, wf *QAWOTable, result, abserr *float64) err.GSLError {
	var (
		area, errsum                          float64
		res_ext, err_ext                      float64
		correc, total_error, truncation_error float64
		ktmin                                 int
		iteration                             int
		table                                 = NewExtrapolationTable()
		cycle                                 float64
		omega                                 float64 = wf.omega
		p                                     float64 = 0.9
		factor                                float64 = 1
		initial_eps, eps                      float64
		error_type                            int
	)
	/* Initialize results */

	workspace.Initialize(a, a)

	*result = 0
	*abserr = 0

	if limit > workspace.limit {
		return err.ERROR("iteration limit exceeds available workspace", err.EINVAL)
	}

	/* Test on accuracy */

	if epsabs <= 0 {
		return err.ERROR("absolute tolerance epsabs must be positive", err.EBADTOL)
	}

	if omega == 0.0 {
		if wf.sine == GSL_INTEG_SINE {
			/* The function sin(w x) f(x) is always zero for w = 0 */

			*result = 0
			*abserr = 0

			return nil
		} else {
			/* The function cos(w x) f(x) is always f(x) for w = 0 */

			return Qagiu(f, a, epsabs, 0.0, cycle_workspace.limit, cycle_workspace, result, abserr)
		}
	}

	if epsabs > gsl.MinFloat64/(1-p) {
		eps = epsabs * (1 - p)
	} else {
		eps = epsabs
	}

	initial_eps = eps

	area = 0
	errsum = 0

	res_ext = 0
	err_ext = gsl.MaxFloat64
	correc = 0

	cycle = (2*math.Floor(math.Abs(omega)) + 1) * gsl.Pi / math.Abs(omega)

	wf.SetLength(cycle)

	for iteration = 0; iteration < limit; iteration++ {
		var (
			area1, error1, reseps, erreps float64
			a1                            float64 = a + float64(iteration)*cycle
			b1                            float64 = a1 + cycle
			epsabs1                       float64 = eps * factor
		)

		status := Qawo(f, a1, epsabs1, 0.0, limit, cycle_workspace, wf, &area1, &error1)
		workspace.AppendInterval(a1, b1, area1, error1)
		factor *= p
		area = area + area1
		errsum = errsum + error1

		/* estimate the truncation error as 50 times the final term */

		truncation_error = 50 * math.Abs(area1)

		total_error = errsum + truncation_error

		if total_error < epsabs && iteration > 4 {
			goto compute_result
		}

		if error1 > correc {
			correc = error1
		}

		if status != nil {
			eps = math.Max(initial_eps, correc*(1.0-p))
		}

		if status != nil && total_error < 10*correc && iteration > 3 {
			goto compute_result
		}

		table.AppendTable(area)

		if table.n < 2 {
			continue
		}

		table.Qelg(&reseps, &erreps)

		ktmin++

		if ktmin >= 15 && err_ext < 0.001*total_error {
			error_type = 4
		}

		if erreps < err_ext {
			ktmin = 0
			err_ext = erreps
			res_ext = reseps

			if err_ext+10*correc <= epsabs {
				break
			}
			if err_ext <= epsabs && 10*correc >= epsabs {
				break
			}
		}

	}

	if iteration == limit {
		error_type = 1
	}

	if err_ext == math.MaxFloat64 {
		goto compute_result
	}

	err_ext = err_ext + 10*correc

	*result = res_ext
	*abserr = err_ext

	if error_type == 0 {
		return nil
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

	if error_type == 4 {
		err_ext = err_ext + truncation_error
	}

	goto return_error

compute_result:

	*result = area
	*abserr = total_error

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
