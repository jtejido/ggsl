/* integration/rational.c
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
	"github.com/lucky-se7en/ggsl/specfunc"
	"math"
)

type Rational struct {
}

func (l Rational) check(n int, params *fixedParams) err.GSLError {
	if math.Abs(params.b-params.a) <= gsl.Float64Eps {
		return err.ERROR("|b - a| too small", err.EDOM)
	} else if params.alpha <= -1.0 {
		return err.ERROR("alpha must be > -1", err.EDOM)
	} else if params.beta >= 0.0 || params.alpha+params.beta+2*float64(n) >= 0.0 || 0.0 >= params.alpha+2*float64(n) {
		return err.ERROR("beta < alpha + beta + 2n < 0 is required", err.EDOM)
	} else if params.a+params.b <= 0.0 {
		return err.ERROR("a + b <= 0 is not allowed", err.EDOM)
	}

	return nil
}

func (l Rational) init(n int, diag, subdiag []float64, params *fixedParams) err.GSLError {
	absum := params.beta + params.alpha
	a1 := params.alpha + 1.0
	aba1 := absum * a1
	ab2i := absum + 2.0

	/* construct the diagonal and subdiagonal elements of Jacobi matrix */

	diag[0] = -a1 / (absum + 2.0)
	subdiag[0] = math.Sqrt(-diag[0] * (params.beta + 1.0) / ((absum + 2.0) * (absum + 3.0)))

	for i := 1; i < n-1; i++ {
		ab2i += 2.0
		diag[i] = (-aba1 - 2.0*float64(i)*(absum+float64(i)+1.0)) / (ab2i * (ab2i - 2.0))
		subdiag[i] = math.Sqrt((float64(i) + 1.0) * (params.alpha + float64(i) + 1.0) / (ab2i - 1.0) * (params.beta + float64(i) + 1.0) / (ab2i * ab2i) * (absum + float64(i) + 1.0) / (ab2i + 1.0))
	}

	diag[n-1] = (-aba1 - 2.0*(float64(n)-1.0)*(absum+float64(n))) / ((absum + 2.0*float64(n)) * (absum + 2.0*float64(n) - 2.0))
	subdiag[n-1] = 0.0

	params.zemu = specfunc.Gamma(params.alpha+1.0) * specfunc.Gamma(-absum-1.0) / specfunc.Gamma(-params.beta)
	params.shft = params.a
	params.slp = params.b + params.a
	params.al = params.alpha
	params.be = params.beta

	return nil
}
