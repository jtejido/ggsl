/* integration/chebyshev.c
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
	"math"
)

type ChebyshevType1 struct {
}

func (ct1 ChebyshevType1) check(n int, params *fixedParams) err.GSLError {
	if math.Abs(params.b-params.a) <= gsl.Float64Eps {
		err.ERROR("|b - a| too small", err.EDOM)
	} else if params.a >= params.b {
		err.ERROR("lower integration limit must be smaller than upper limit", err.EDOM)
	}

	return nil

}

func (ct1 ChebyshevType1) init(n int, diag, subdiag []float64, params *fixedParams) err.GSLError {
	/* construct the diagonal and subdiagonal elements of Jacobi matrix */
	diag[0] = 0.0
	subdiag[0] = gsl.SqrtHalf
	for i := 1; i < n; i++ {
		diag[i] = 0.0
		subdiag[i] = 0.5
	}

	params.zemu = gsl.Pi
	params.shft = 0.5 * (params.b + params.a)
	params.slp = 0.5 * (params.b - params.a)
	params.al = -0.5
	params.be = -0.5

	return nil
}
