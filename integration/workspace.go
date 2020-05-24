/* integration/workspace.c
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
	"github.com/lucky-se7en/ggsl/err"
)

type Workspace struct {
	limit         int
	size          int
	nrmax         int
	i             int
	maximum_level int
	alist         []float64
	blist         []float64
	rlist         []float64
	elist         []float64
	order         []int
	level         []int
}

func NewWorkspace(n int) (*Workspace, err.GSLError) {
	if n == 0 {
		return nil, err.ERROR("workspace length n must be positive integer", err.EDOM)
	}

	w := &Workspace{
		size:          0,
		limit:         n,
		maximum_level: 0,
		alist:         make([]float64, n),
		blist:         make([]float64, n),
		elist:         make([]float64, n),
		rlist:         make([]float64, n),
		order:         make([]int, n),
		level:         make([]int, n),
	}

	return w, nil
}

func (workspace *Workspace) Initialize(a, b float64) {
	workspace.size = 0
	workspace.nrmax = 0
	workspace.i = 0
	workspace.alist[0] = a
	workspace.blist[0] = b
	workspace.rlist[0] = 0.0
	workspace.elist[0] = 0.0
	workspace.order[0] = 0
	workspace.level[0] = 0

	workspace.maximum_level = 0
}

func (workspace *Workspace) SetInitialResult(result, reserr float64) {
	workspace.size = 1
	workspace.rlist[0] = result
	workspace.elist[0] = reserr
}

func (workspace *Workspace) Retrieve(a, b, r, e *float64) {
	i := workspace.i
	*a = workspace.alist[i]
	*b = workspace.blist[i]
	*r = workspace.rlist[i]
	*e = workspace.elist[i]
}

func (workspace *Workspace) Update(a1, b1, area1, error1, a2, b2, area2, error2 float64) {
	i_max := workspace.i
	i_new := workspace.size

	new_level := workspace.level[i_max] + 1

	/* append the newly-created intervals to the list */

	if error2 > error1 {
		workspace.alist[i_max] = a2
		workspace.rlist[i_max] = area2
		workspace.elist[i_max] = error2
		workspace.level[i_max] = new_level

		workspace.alist[i_new] = a1
		workspace.blist[i_new] = b1
		workspace.rlist[i_new] = area1
		workspace.elist[i_new] = error1
		workspace.level[i_new] = new_level
	} else {
		workspace.blist[i_max] = b1 /* alist[maxerr] is already == a1 */
		workspace.rlist[i_max] = area1
		workspace.elist[i_max] = error1
		workspace.level[i_max] = new_level

		workspace.alist[i_new] = a2
		workspace.blist[i_new] = b2
		workspace.rlist[i_new] = area2
		workspace.elist[i_new] = error2
		workspace.level[i_new] = new_level
	}

	workspace.size++
	if new_level > workspace.maximum_level {
		workspace.maximum_level = new_level
	}
	workspace.Qpsrt()
}

func (workspace *Workspace) AppendInterval(a1, b1, area1, error1 float64) {
	i_new := workspace.size
	workspace.alist[i_new] = a1
	workspace.blist[i_new] = b1
	workspace.rlist[i_new] = area1
	workspace.elist[i_new] = error1
	workspace.order[i_new] = i_new
	workspace.level[i_new] = 0

	workspace.size++
}

func (workspace *Workspace) Sort() {
	elist := workspace.elist
	order := workspace.order

	nint := workspace.size

	for i := 0; i < nint; i++ {
		i1 := order[i]
		e1 := elist[i1]
		i_max := i1
		for j := i + 1; j < nint; j++ {
			i2 := order[j]
			e2 := elist[i2]

			if e2 >= e1 {
				i_max = i2
				e1 = e2
			}
		}

		if i_max != i1 {
			order[i] = order[i_max]
			order[i_max] = i1
		}
	}

	workspace.i = order[0]
}

func (workspace *Workspace) IsLargeInterval() bool {
	i := workspace.i
	level := workspace.level

	if level[i] < workspace.maximum_level {
		return true
	} else {
		return false
	}
}

func (workspace *Workspace) Qpsrt() {
	last := workspace.size - 1
	limit := workspace.limit

	var errmax float64
	var errmin float64
	var i, k, top int

	i_nrmax := workspace.nrmax
	i_maxerr := workspace.order[i_nrmax]

	/* Check whether the list contains more than two error estimates */

	if last < 2 {
		workspace.order[0] = 0
		workspace.order[1] = 1
		workspace.i = i_maxerr
		return
	}

	errmax = workspace.elist[i_maxerr]

	/* This part of the routine is only executed if, due to a difficult
	   integrand, subdivision increased the error estimate. In the normal
	   case the insert procedure should start after the nrmax-th largest
	   error estimate. */

	for i_nrmax > 0 && errmax > workspace.elist[workspace.order[i_nrmax-1]] {
		workspace.order[i_nrmax] = workspace.order[i_nrmax-1]
		i_nrmax--
	}

	/* Compute the number of elements in the list to be maintained in
	   descending order. This number depends on the number of
	   subdivisions still allowed. */

	if last < (limit/2 + 2) {
		top = last
	} else {
		top = limit - last + 1
	}

	/* Insert errmax by traversing the list top-down, starting
	   comparison from the element elist(order(i_nrmax+1)). */

	i = i_nrmax + 1

	/* The order of the tests in the following line is important to
	   prevent a segmentation fault */

	for i < top && errmax < workspace.elist[workspace.order[i]] {
		workspace.order[i-1] = workspace.order[i]
		i++
	}

	workspace.order[i-1] = i_maxerr

	/* Insert errmin by traversing the list bottom-up */

	errmin = workspace.elist[last]

	k = top - 1

	for k > i-2 && errmin >= workspace.elist[workspace.order[k]] {
		workspace.order[k+1] = workspace.order[k]
		k--
	}

	workspace.order[k+1] = last

	/* Set i_max and e_max */

	i_maxerr = workspace.order[i_nrmax]

	workspace.i = i_maxerr
	workspace.nrmax = i_nrmax
}

/* The smallest interval has the largest error.  Before bisecting
   decrease the sum of the errors over the larger intervals
   (error_over_large_intervals) and perform extrapolation. */
func (workspace *Workspace) IsIncreaseNrmax() bool {
	id := workspace.nrmax
	var jupbnd int

	level := workspace.level
	order := workspace.order

	limit := workspace.limit
	last := workspace.size - 1

	if last > (1 + limit/2) {
		jupbnd = limit + 1 - last
	} else {
		jupbnd = last
	}

	for k := id; k <= jupbnd; k++ {
		i_max := order[workspace.nrmax]

		workspace.i = i_max

		if level[i_max] < workspace.maximum_level {
			return true
		}

		workspace.nrmax++

	}
	return false
}

func (workspace *Workspace) ResetNrmax() {
	workspace.nrmax = 0
	workspace.i = workspace.order[0]
}

func (workspace *Workspace) SumResults() (sum float64) {
	n := workspace.size

	for k := 0; k < n; k++ {
		sum += workspace.rlist[k]
	}

	return
}
