/* integration/jacobi.c
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
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"github.com/lucky-se7en/ggsl/specfunc"
	"math"
)

type Jacobi struct {
}

func (j Jacobi) check(n int, params *fixedParams) err.GSLError {
	if math.Abs(params.b-params.a) <= gsl.Float64Eps {
		return err.ERROR("|b - a| too small", err.EDOM)
	} else if params.a >= params.b {
		return err.ERROR("lower integration limit must be smaller than upper limit", err.EDOM)
	} else if params.alpha <= -1.0 || params.beta <= -1.0 {
		return err.ERROR("alpha and beta must be > -1", err.EDOM)
	}

	return nil
}

func (j Jacobi) init(n int, diag, subdiag []float64, params *fixedParams) err.GSLError {
	absum := params.beta + params.alpha
	abdiff := params.beta - params.alpha
	a2b2 := absum * abdiff /* beta^2 - alpha^2 */

	/* construct the diagonal and subdiagonal elements of Jacobi matrix */
	diag[0] = abdiff / (absum + 2.0)
	subdiag[0] = 2.0 * math.Sqrt((params.alpha+1.0)*(params.beta+1.0)/(absum+3.0)) / (absum + 2.0)
	for i := 1; i < n; i++ {
		diag[i] = a2b2 / ((absum + 2.0*float64(i)) * (absum + 2.0*float64(i) + 2.0))
		subdiag[i] = math.Sqrt(4.0*(float64(i)+1.0)*(params.alpha+float64(i)+1.0)*(params.beta+float64(i)+1.0)*(absum+float64(i)+1.0)/(math.Pow((absum+2.0*float64(i)+2.0), 2.0)-1.0)) / (absum + 2.0*float64(i) + 2.0)
	}

	params.zemu = math.Pow(2.0, absum+1.0) * specfunc.Gamma(params.alpha+1.0) * specfunc.Gamma(params.beta+1.0) / specfunc.Gamma(absum+2.0)
	params.shft = 0.5 * (params.b + params.a)
	params.slp = 0.5 * (params.b - params.a)
	params.al = params.alpha
	params.be = params.beta

	return nil
}
