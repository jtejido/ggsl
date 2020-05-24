/* integration/qc25c.c
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
	"math"
)

type fnCauchy struct {
	f           gsl.Function
	singularity float64
}

func (fnc *fnCauchy) Evaluate(x float64) float64 {
	return fnc.f.Evaluate(x) / (x - fnc.singularity)
}

func qc25c(f gsl.Function, a, b, c float64, result, abserr *float64, err_reliable *int) {
	cc := (2*c - b - a) / (b - a)

	if math.Abs(cc) > 1.1 {
		var resabs, resasc float64
		weighted_function := &fnCauchy{f, c}

		Qk15(weighted_function, a, b, result, abserr, &resabs, &resasc)

		if *abserr == resasc {
			*err_reliable = 0
		} else {
			*err_reliable = 1
		}

		return
	} else {
		cheb12 := make([]float64, 13)
		cheb24 := make([]float64, 25)
		moment := make([]float64, 25)
		var res12, res24 float64

		Qcheb(f, a, b, cheb12, cheb24)
		computeMoments(cc, moment)

		for i := 0; i < 13; i++ {
			res12 += cheb12[i] * moment[i]
		}

		for i := 0; i < 25; i++ {
			res24 += cheb24[i] * moment[i]
		}

		*result = res24
		*abserr = math.Abs(res24 - res12)
		*err_reliable = 0

	}
}

func computeMoments(cc float64, moment []float64) {
	a0 := math.Log(math.Abs((1.0 - cc) / (1.0 + cc)))
	a1 := 2 + a0*cc

	moment[0] = a0
	moment[1] = a1

	for k := 2; k < 25; k++ {
		var a2 float64

		if (k % 2) == 0 {
			a2 = 2.0*cc*a1 - a0
		} else {
			km1 := float64(k) - 1.0
			a2 = 2.0*cc*a1 - a0 - 4.0/(km1*km1-1.0)
		}

		moment[k] = a2

		a0 = a1
		a1 = a2
	}
}
