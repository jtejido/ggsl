/* integration/qelg.c
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
	"math"
)

type extrapolationTable struct {
	n      int
	rlist2 []float64
	nres   int
	res3la []float64
}

func NewExtrapolationTable() *extrapolationTable {
	return &extrapolationTable{n: 0, rlist2: make([]float64, 52), nres: 0, res3la: make([]float64, 3)}
}

func (table *extrapolationTable) AppendTable(y float64) {
	n := table.n
	table.rlist2[n] = y
	table.n++
}

func (table *extrapolationTable) Qelg(result, abserr *float64) {
	epstab := table.rlist2
	res3la := table.res3la
	n := table.n - 1

	current := epstab[n]

	absolute := gsl.MaxFloat64
	relative := 5 * gsl.Float64Eps * math.Abs(current)

	newelm := n / 2
	n_orig := n
	n_final := n

	nres_orig := table.nres

	*result = current
	*abserr = gsl.MaxFloat64

	if n < 2 {
		*result = current
		*abserr = math.Max(absolute, relative)
		return
	}

	epstab[n+2] = epstab[n]
	epstab[n] = gsl.MaxFloat64

	for i := 0; i < newelm; i++ {
		res := epstab[n-2*i+2]
		e0 := epstab[n-2*i-2]
		e1 := epstab[n-2*i-1]
		e2 := res

		e1abs := math.Abs(e1)
		delta2 := e2 - e1
		err2 := math.Abs(delta2)
		tol2 := math.Max(math.Abs(e2), e1abs) * gsl.Float64Eps
		delta3 := e1 - e0
		err3 := math.Abs(delta3)
		tol3 := math.Max(e1abs, math.Abs(e0)) * gsl.Float64Eps

		var e3, delta1, err1, tol1, ss float64

		if err2 <= tol2 && err3 <= tol3 {
			/* If e0, e1 and e2 are equal to within machine accuracy,
			   convergence is assumed.  */

			*result = res
			absolute = err2 + err3
			relative = 5 * gsl.Float64Eps * math.Abs(res)
			*abserr = math.Max(absolute, relative)
			return
		}

		e3 = epstab[n-2*i]
		epstab[n-2*i] = e1
		delta1 = e1 - e3
		err1 = math.Abs(delta1)
		tol1 = math.Max(e1abs, math.Abs(e3)) * gsl.Float64Eps

		/* If two elements are very close to each other, omit a part of
		   the table by adjusting the value of n */

		if err1 <= tol1 || err2 <= tol2 || err3 <= tol3 {
			n_final = 2 * i
			break
		}

		ss = (1/delta1 + 1/delta2) - 1/delta3

		/* Test to detect irregular behaviour in the table, and
		   eventually omit a part of the table by adjusting the value of
		   n. */

		if math.Abs(ss*e1) <= 0.0001 {
			n_final = 2 * i
			break
		}

		/* Compute a new element and eventually adjust the value of
		   result. */

		res = e1 + 1/ss
		epstab[n-2*i] = res

		error := err2 + math.Abs(res-e2) + err3

		if error <= *abserr {
			*abserr = error
			*result = res
		}

	}

	/* Shift the table */
	limexp := 50 - 1

	if n_final == limexp {
		n_final = 2 * (limexp / 2)
	}

	if n_orig%2 == 1 {
		for i := 0; i <= newelm; i++ {
			epstab[1+i*2] = epstab[i*2+3]
		}
	} else {
		for i := 0; i <= newelm; i++ {
			epstab[i*2] = epstab[i*2+2]
		}
	}

	if n_orig != n_final {
		for i := 0; i <= n_final; i++ {
			epstab[i] = epstab[n_orig-n_final+i]
		}
	}

	table.n = n_final + 1

	if nres_orig < 3 {
		res3la[nres_orig] = *result
		*abserr = gsl.MaxFloat64
	} else { /* Compute error estimate */
		*abserr = (math.Abs(*result-res3la[2]) + math.Abs(*result-res3la[1]) + math.Abs(*result-res3la[0]))
		res3la[0] = res3la[1]
		res3la[1] = res3la[2]
		res3la[2] = *result
	}

	/* In QUADPACK the variable table->nres is incremented at the top of
	   qelg, so it increases on every call. This leads to the array
	   res3la being accessed when its elements are still undefined, so I
	   have moved the update to this point so that its value more
	   useful. */

	table.nres = nres_orig + 1

	*abserr = math.Max(*abserr, 5*gsl.Float64Eps*math.Abs(*result))

	return
}
