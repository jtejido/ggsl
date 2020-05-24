/* integration/hermite.c
 *
 * Copyright (C) 2017 Patrick Alken
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
	"github.com/jtejido/ggsl/err"
	"github.com/jtejido/ggsl/specfunc"
	"math"
)

type Hermite struct {
}

func (l Hermite) check(n int, params *fixedParams) err.GSLError {
	if params.b <= 0.0 {
		return err.ERROR("b must be positive", err.EDOM)
	} else if params.alpha <= -1.0 {
		return err.ERROR("alpha must be > -1", err.EDOM)
	}

	return nil
}

func (l Hermite) init(n int, diag, subdiag []float64, params *fixedParams) err.GSLError {
	/* construct the diagonal and subdiagonal elements of Jacobi matrix */
	for i := 1; i <= n; i++ {
		diag[i-1] = 0.0
		subdiag[i-1] = math.Sqrt(0.5 * (float64(i) + params.alpha*float64(i%2)))
	}

	params.zemu = specfunc.Gamma(0.5 * (params.alpha + 1.0))
	params.shft = params.a
	params.slp = 1.0 / math.Sqrt(params.b)
	params.al = params.alpha
	params.be = 0.0

	return nil
}
