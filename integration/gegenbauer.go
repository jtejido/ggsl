/* integration/gegenbauer.c
 *
 * Copyright (C) 2017 Konrad Griessinger, Patrick Alken
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

/*
 * The code in this module is based on IQPACK, specifically the LGPL
 * implementation found in HERMITE_RULE:
 * https://people.sc.fsu.edu/~jburkardt/c_src/hermite_rule/hermite_rule.html
 */
package integration

import (
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"github.com/jtejido/ggsl/specfunc"
	"math"
)

type Gegenbauer struct {
}

func (g Gegenbauer) check(n int, params *fixedParams) err.GSLError {
	if math.Abs(params.b-params.a) <= gsl.Float64Eps {
		return err.ERROR("|b - a| too small", err.EDOM)
	} else if params.a >= params.b {
		return err.ERROR("lower integration limit must be smaller than upper limit", err.EDOM)
	} else if params.alpha <= -1.0 {
		return err.ERROR("alpha must be > -1", err.EDOM)
	}

	return nil

}

func (g Gegenbauer) init(n int, diag, subdiag []float64, params *fixedParams) err.GSLError {
	/* construct the diagonal and subdiagonal elements of Jacobi matrix */
	diag[0] = 0.0
	subdiag[0] = math.Sqrt(1.0 / (2.0*params.alpha + 3.0))
	for i := 1; i < n; i++ {
		diag[i] = 0.0
		subdiag[i] = math.Sqrt((float64(i) + 1.0) * (2.0*params.alpha + float64(i) + 1.0) / (4.0*(params.alpha+float64(i)+1.0)*(params.alpha+float64(i)+1.0) - 1.0))
	}

	params.zemu = math.SqrtPi * specfunc.Gamma(params.alpha+1.0) / specfunc.Gamma(params.alpha+1.5)
	params.shft = 0.5 * (params.b + params.a)
	params.slp = 0.5 * (params.b - params.a)
	params.al = params.alpha
	params.be = params.alpha

	return nil
}
