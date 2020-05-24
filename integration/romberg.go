/* integration/romberg.c
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

// This function allocates a workspace for Romberg integration, specifying a maximum of n iterations, or divisions of the interval.
// Since the number of divisions is 2^n + 1, n can be kept relatively small (i.e. 10 or 20). It is capped at a maximum value of 30 to prevent overflow.
// The size of the workspace is O(2n).
type RombergWorkspace struct {
	n            int
	work1, work2 []float64
}

func NewRombergWorkspace(n int) (*RombergWorkspace, err.GSLError) {

	/* check inputs */
	if n < 1 {
		return nil, err.ERROR("workspace size n must be at least 1", err.EDOM)
	}
	w := new(RombergWorkspace)
	/* ceiling on n, since the number of points is 2^n + 1 */
	w.n = int(gsl.Min(n, 30))

	w.work1 = make([]float64, n)
	w.work2 = make([]float64, n)
	return w, nil
}

// This function integrates f(x), specified by f, from a to b, storing the answer in result. At each step in the iteration, convergence is tested by checking:
// | I_k - I_{k-1} | \le \textrm{max} \left( epsabs, epsrel \times |I_k| \right)
// where I_k is the current approximation and I_{k-1} is the approximation of the previous iteration.
// If the method does not converge within the previously specified n iterations,
// the function stores the best current estimate in result and returns GSL_EMAXITER.
// If the method converges, the function returns GSL_SUCCESS. The total number of function evaluations is returned in neval.
func Romberg(f gsl.Function, a, b, epsabs, epsrel float64, result *float64, neval *int, w *RombergWorkspace) err.GSLError {
	if epsabs < 0.0 {
		return err.ERROR("epsabs must be non-negative", err.EDOM)
	} else if epsrel < 0.0 {
		return err.ERROR("epsrel must be non-negative", err.EDOM)
	} else {
		n := w.n
		Rp := w.work1 /* previous row */
		Rc := w.work2 /* current row */
		var Rtmp []float64
		h := 0.5 * (b - a) /* step size */

		/* R(0,0) */
		Rp[0] = h * (f.Evaluate(a) + f.Evaluate(b))
		*neval = 2

		/*ROMBERG_PRINT_ROW((size_t) 0, Rp);*/
		for i := 1; i < n; i++ {
			sum := 0.0
			two_i := 1 << i /* 2^i */

			for j := 1; j < two_i; j += 2 {
				sum += f.Evaluate(a + float64(j)*h)
				*neval++
			}

			/* R(i,0) */
			Rc[0] = sum*h + 0.5*Rp[0]

			four_j := 4.0
			for j := 1; j <= i; j++ {
				Rc[j] = (four_j*Rc[j-1] - Rp[j-1]) / (four_j - 1.0)
				four_j *= 4.0
			}

			/*ROMBERG_PRINT_ROW(i, Rc);*/

			/*
			 * compute difference between current and previous result and
			 * check for convergence
			 */
			err := math.Abs(Rc[i] - Rp[i-1])
			if (err < epsabs) || (err < epsrel*math.Abs(Rc[i])) {
				*result = Rc[i]
				return nil
			}

			/* swap Rp and Rc */
			Rtmp = Rp
			Rp = Rc
			Rc = Rtmp

			h *= 0.5
		}

		/* did not converge - return best guess */
		*result = Rp[n-1]

		return err.MaxIteration()
	}
}
