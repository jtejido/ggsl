/* integration/qmomof.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
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
	"math"
)

type QAWO int

const (
	GSL_INTEG_COSINE QAWO = iota
	GSL_INTEG_SINE
)

type slice struct {
	data   []float64
	offset int
}

type QAWOTable struct {
	n             int
	omega, L, par float64
	sine          QAWO
	chebmo        *slice
}

func NewQAWOTable(omega, L float64, sine QAWO, n int) (*QAWOTable, err.GSLError) {
	if n == 0 {
		return nil, err.ERROR("table length n must be positive integer", err.EDOM)
	}

	t := new(QAWOTable)
	t.n = n
	t.sine = sine
	t.omega = omega
	t.L = L
	t.par = 0.5 * omega * L
	t.chebmo = &slice{make([]float64, 25*n), 0}

	/* precompute the moments */

	{
		scale := 1.0

		for i := 0; i < t.n; i++ {
			t.chebmo.offset = 25 * i
			computeMomentsQAWO(t.par*scale, t.chebmo)
			scale *= 0.5
		}
	}

	return t, nil
}

func (t *QAWOTable) Set(omega, L float64, sine QAWO) err.GSLError {
	t.omega = omega
	t.sine = sine
	t.L = L
	t.par = 0.5 * omega * L

	/* recompute the moments */

	{
		scale := 1.0

		for i := 0; i < t.n; i++ {
			t.chebmo.offset = 25 * i
			computeMomentsQAWO(t.par*scale, t.chebmo)
			scale *= 0.5
		}
	}

	return nil
}

func (t *QAWOTable) SetLength(L float64) err.GSLError {
	/* return immediately if the length is the same as the old length */

	if L == t.L {
		return nil
	}

	/* otherwise reset the table and compute the new parameters */

	t.L = L
	t.par = 0.5 * t.omega * L

	/* recompute the moments */

	{
		scale := 1.0

		for i := 0; i < t.n; i++ {
			t.chebmo.offset = 25 * i
			computeMomentsQAWO(t.par*scale, t.chebmo)
			scale *= 0.5
		}
	}

	return nil
}

func computeMomentsQAWO(par float64, chebmo *slice) {
	var (
		v      = &slice{make([]float64, 28), 0}
		d      = make([]float64, 25)
		d1     = make([]float64, 25)
		d2     = make([]float64, 25)
		noeq   = 25
		par2   = par * par
		par4   = par2 * par2
		par22  = par2 + 2.0
		sinpar = math.Sin(par)
		cospar = math.Cos(par)
	)

	/* compute the chebyschev moments with respect to cosine */
	ac := 8 * cospar
	as := 24 * par * sinpar

	v.data[v.offset+0] = 2 * sinpar / par
	v.data[v.offset+1] = (8*cospar + (2*par2-8)*sinpar/par) / par2
	v.data[v.offset+2] = (32*(par2-12)*cospar + (2*((par2-80)*par2+192)*sinpar)/par) / par4

	if math.Abs(par) <= 24 {
		/* compute the moments as the solution of a boundary value
		   problem using the asyptotic expansion as an endpoint */
		var an2, ass, asap float64
		an := 6.

		for k := 0; k < noeq-1; k++ {
			an2 = an * an
			d[k] = -2 * (an2 - 4) * (par22 - 2*an2)
			d2[k] = (an - 1) * (an - 2) * par2
			d1[k+1] = (an + 3) * (an + 4) * par2
			v.data[v.offset+k+3] = as - (an2-4)*ac
			an = an + 2.0
		}

		an2 = an * an

		d[noeq-1] = -2 * (an2 - 4) * (par22 - 2*an2)
		v.data[v.offset+noeq+2] = as - (an2-4)*ac
		v.data[v.offset+3] = v.data[v.offset+3] - 56*par2*v.data[v.offset+2]

		ass = par * sinpar
		asap = (((((210*par2-1)*cospar-(105*par2-63)*ass)/an2-(1-15*par2)*cospar+15*ass)/an2-cospar+3*ass)/an2 - cospar) / an2
		v.data[v.offset+noeq+2] = v.data[v.offset+noeq+2] - 2*asap*par2*(an-1)*(an-2)
		v.offset = 3
		dgtsl(noeq, d1, d, d2, v)

	} else {
		/* compute the moments by forward recursion */
		an := 4.

		for k := 3; k < 13; k++ {
			an2 := an * an
			v.data[v.offset+k] = ((an2-4)*(2*(par22-2*an2)*v.data[v.offset+k-1]-ac) + as - par2*(an+1)*(an+2)*v.data[v.offset+k-2]) / (par2 * (an - 1) * (an - 2))
			an = an + 2.0
		}
	}

	for i := 0; i < 13; i++ {
		chebmo.data[chebmo.offset+2*i] = v.data[v.offset+i]
	}

	/* compute the chebyschev moments with respect to sine */

	v.data[v.offset+0] = 2 * (sinpar - par*cospar) / par2
	v.data[v.offset+1] = (18-48/par2)*sinpar/par2 + (-2+48/par2)*cospar/par

	ac = -24 * par * cospar
	as = -8 * sinpar

	if math.Abs(par) <= 24 {
		/* compute the moments as the solution of a boundary value
		   problem using the asyptotic expansion as an endpoint */
		var an2, ass, asap float64
		an := 5.

		for k := 0; k < noeq-1; k++ {
			an2 = an * an
			d[k] = -2 * (an2 - 4) * (par22 - 2*an2)
			d2[k] = (an - 1) * (an - 2) * par2
			d1[k+1] = (an + 3) * (an + 4) * par2
			v.data[v.offset+k+2] = ac + (an2-4)*as
			an = an + 2.0
		}

		an2 = an * an

		d[noeq-1] = -2 * (an2 - 4) * (par22 - 2*an2)
		v.data[v.offset+noeq+1] = ac + (an2-4)*as
		v.data[v.offset+2] = v.data[v.offset+2] - 42*par2*v.data[v.offset+1]

		ass = par * cospar
		asap = (((((105*par2-63)*ass-(210*par2-1)*sinpar)/an2+(15*par2-1)*sinpar-15*ass)/an2-sinpar-3*ass)/an2 - sinpar) / an2
		v.data[v.offset+noeq+1] = v.data[v.offset+noeq+1] - 2*asap*par2*(an-1)*(an-2)
		v.offset = 2
		dgtsl(noeq, d1, d, d2, v)

	} else {
		/* compute the moments by forward recursion */
		an := 3.
		for k := 2; k < 12; k++ {
			an2 := an * an
			v.data[v.offset+k] = ((an2-4)*(2*(par22-2*an2)*v.data[v.offset+k-1]+as) + ac - par2*(an+1)*(an+2)*v.data[v.offset+k-2]) / (par2 * (an - 1) * (an - 2))
			an = an + 2.0
		}
	}

	for i := 0; i < 12; i++ {
		chebmo.data[chebmo.offset+2*i+1] = v.data[v.offset+i]
	}

	chebmo.offset = 0
}

func dgtsl(n int, c, d, e []float64, b *slice) err.GSLError {
	/* solves a tridiagonal matrix A x = b

	   c[1 .. n - 1]   subdiagonal of the matrix A
	   d[0 .. n - 1]   diagonal of the matrix A
	   e[0 .. n - 2]   superdiagonal of the matrix A

	   b[0 .. n - 1]   right hand side, replaced by the solution vector x */

	c[0] = d[0]

	if n == 0 {
		return nil
	}

	if n == 1 {
		b.data[b.offset+0] = b.data[b.offset+0] / d[0]
		return nil
	}

	d[0] = e[0]
	e[0] = 0
	e[n-1] = 0

	for k := 0; k < n-1; k++ {
		k1 := k + 1

		if math.Abs(c[k1]) >= math.Abs(c[k]) {
			{
				t := c[k1]
				c[k1] = c[k]
				c[k] = t
			}
			{
				t := d[k1]
				d[k1] = d[k]
				d[k] = t
			}
			{
				t := e[k1]
				e[k1] = e[k]
				e[k] = t
			}
			{
				t := b.data[b.offset+k1]
				b.data[b.offset+k1] = b.data[b.offset+k]
				b.data[b.offset+k] = t
			}
		}

		if c[k] == 0 {
			return err.Failure()
		}

		{
			t := -c[k1] / c[k]

			c[k1] = d[k1] + t*d[k]
			d[k1] = e[k1] + t*e[k]
			e[k1] = 0
			b.data[b.offset+k1] = b.data[b.offset+k1] + t*b.data[b.offset+k]
		}

	}

	if c[n-1] == 0 {
		return err.Failure()
	}

	b.data[b.offset+n-1] = b.data[b.offset+n-1] / c[n-1]

	b.data[b.offset+n-2] = (b.data[b.offset+n-2] - d[n-2]*b.data[b.offset+n-1]) / c[n-2]

	for k := n; k > 2; k-- {
		kb := k - 3
		b.data[b.offset+kb] = (b.data[b.offset+kb] - d[kb]*b.data[b.offset+kb+1] - e[kb]*b.data[b.offset+kb+2]) / c[kb]
	}

	b.offset = 0
	return nil
}
