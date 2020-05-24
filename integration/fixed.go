/* integration/fixed.c
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

/* the code in this module performs fixed-point quadrature calculations for
 * integrands and is based on IQPACK */
package integration

import (
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"math"
)

type fixedParams struct {
	alpha, beta, a, b, zemu, shft, slp, al, be float64
}

type FixedType interface {
	check(n int, params *fixedParams) err.GSLError
	init(n int, diag, subdiag []float64, params *fixedParams) err.GSLError
}

// This workspace is used for fixed point quadrature rules.
type FixedWorkspace struct {
	n       int       /* number of nodes/weights */
	weights []float64 /* quadrature weights */
	x       []float64 /* quadrature nodes */
	diag    []float64 /* diagonal of Jacobi matrix */
	subdiag []float64 /* subdiagonal of Jacobi matrix */
	t       FixedType
}

func NewFixedWorkspace(t FixedType, n int, a, b, alpha, beta float64) (*FixedWorkspace, err.GSLError) {
	/* check inputs */
	if n < 1 {
		return nil, err.ERROR("workspace size n must be at least 1", err.EDOM)
	}

	w := new(FixedWorkspace)
	w.weights = make([]float64, n)
	w.x = make([]float64, n)
	w.diag = make([]float64, n)
	w.subdiag = make([]float64, n)
	w.n = n
	w.t = t

	/* compute quadrature weights and nodes */
	if w.compute(a, b, alpha, beta) != nil {
		return nil, err.ERROR("error in integration parameters", err.EDOM)
	}

	return w, nil
}

func (w *FixedWorkspace) N() int {
	return w.n
}

func (w *FixedWorkspace) Nodes() []float64 {
	return w.x
}

func (w *FixedWorkspace) Weights() []float64 {
	return w.weights
}

/*
  fixed_compute()
  Compute quadrature weights and nodes
*/
func (w *FixedWorkspace) compute(a, b, alpha, beta float64) err.GSLError {
	var er err.GSLError
	n := w.n
	params := new(fixedParams)
	params.a = a
	params.b = b
	params.alpha = alpha
	params.beta = beta

	/* check input parameters */
	er = w.t.check(n, params)
	if er != nil {
		return er
	}

	/* initialize Jacobi matrix */
	er = w.t.init(n, w.diag, w.subdiag, params)
	if er != nil {
		return er
	}

	if params.zemu <= 0.0 {
		return err.ERROR("zeroth moment must be positive", err.EINVAL)
	}

	for i := 0; i < n; i++ {
		w.x[i] = w.diag[i]
	}

	w.weights[0] = math.Sqrt(params.zemu)

	for i := 1; i < n; i++ {
		w.weights[i] = 0.0
	}

	/* diagonalize the Jacobi matrix */
	er = imtqlx(n, w.x, w.subdiag, w.weights)
	if er != nil {
		return er
	}

	for i := 0; i < n; i++ {
		w.weights[i] = w.weights[i] * w.weights[i]
	}

	/*
	 * The current weights and nodes are valid for a = 0, b = 1.
	 * Now scale them for arbitrary a,b
	 */
	{
		p := math.Pow(params.slp, params.al+params.be+1.0)
		for k := 0; k < n; k++ {
			w.x[k] = params.shft + params.slp*w.x[k]
			w.weights[k] = w.weights[k] * p
		}
	}

	return nil
}

/******************************************************************************/
/*
  Purpose:

    IMTQLX diagonalizes a symmetric tridiagonal matrix.

  Discussion:

    This routine is a slightly modified version of the EISPACK routine to
    perform the implicit QL algorithm on a symmetric tridiagonal matrix.

    The authors thank the authors of EISPACK for permission to use this
    routine.

    It has been modified to produce the product Q' * Z, where Z is an input
    vector and Q is the orthogonal matrix diagonalizing the input matrix.
    The changes consist (essentially) of applying the orthogonal transformations
    directly to Z as they are generated.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

    Roger Martin, James Wilkinson,
    The Implicit QL Algorithm,
    Numerische Mathematik,
    Volume 12, Number 5, December 1968, pages 377-383.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D(N), the diagonal entries of the matrix.
    On output, the information in D has been overwritten.

    Input/output, double E(N), the subdiagonal entries of the
    matrix, in entries E(1) through E(N-1).  On output, the information in
    E has been overwritten.

    Input/output, double Z(N).  On input, a vector.  On output,
    the value of Q' * Z, where Q is the matrix that diagonalizes the
    input symmetric tridiagonal matrix.
*/

func imtqlx(n int, d, e, z []float64) err.GSLError {
	var b, c, f, g float64
	var i, ii int
	itn := 30
	var j, k, l, m, mml int
	var p, r, s float64

	if n == 1 {
		return nil
	}

	e[n-1] = 0.0

	for l = 1; l <= n; l++ {
		j = 0
		for {
			for m = l; m <= n; m++ {
				if m == n {
					break
				}

				if math.Abs(e[m-1]) <= gsl.Float64Eps*(math.Abs(d[m-1])+math.Abs(d[m])) {
					break
				}
			}
			p = d[l-1]
			if m == l {
				break
			}
			if itn <= j {
				return err.MaxIteration()
			}
			j = j + 1
			g = (d[l] - p) / (2.0 * e[l-1])
			r = math.Sqrt(g*g + 1.0)
			g = d[m-1] - p + e[l-1]/(g+math.Abs(r)*gsl.Sign(g))
			s = 1.0
			c = 1.0
			p = 0.0
			mml = m - l

			for ii = 1; ii <= mml; ii++ {
				i = m - ii
				f = s * e[i-1]
				b = c * e[i-1]

				if math.Abs(g) <= math.Abs(f) {
					c = g / f
					r = math.Sqrt(c*c + 1.0)
					e[i] = f * r
					s = 1.0 / r
					c = c * s
				} else {
					s = f / g
					r = math.Sqrt(s*s + 1.0)
					e[i] = g * r
					c = 1.0 / r
					s = s * c
				}
				g = d[i] - p
				r = (d[i-1]-g)*s + 2.0*c*b
				p = s * r
				d[i] = g + p
				g = c*r - b
				f = z[i]
				z[i] = s*z[i-1] + c*f
				z[i-1] = c*z[i-1] - s*f
			}
			d[l-1] = d[l-1] - p
			e[l-1] = g
			e[m-1] = 0.0
		}
	}
	/*
	   Sorting.
	*/
	for ii = 2; ii <= m; ii++ {
		i = ii - 1
		k = i
		p = d[i-1]

		for j = ii; j <= n; j++ {
			if d[j-1] < p {
				k = j
				p = d[j-1]
			}
		}

		if k != i {
			d[k-1] = d[i-1]
			d[i-1] = p
			p = z[i-1]
			z[i-1] = z[k-1]
			z[k-1] = p
		}
	}

	return nil
}

// This function integrates the function f(x) provided in func using previously computed fixed quadrature rules.
// The integral is approximated as
//  \sum_{i=1}^n w_i f(x_i)
// where w_i are the quadrature weights and x_i are the quadrature nodes computed previously by gsl_integration_fixed_alloc().
// The sum is stored in result on output.
func Fixed(f gsl.Function, result *float64, w *FixedWorkspace) err.GSLError {
	n := w.n
	sum := 0.0

	for i := 0; i < n; i++ {
		fi := f.Evaluate(w.x[i])
		sum += w.weights[i] * fi
	}

	*result = sum

	return nil
}
