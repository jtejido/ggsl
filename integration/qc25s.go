/* integration/qc25s.c
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

type fnQAWSParams struct {
	f     gsl.Function
	a, b  float64
	table *QAWSTable
}

type fnQAWSR struct {
	p *fnQAWSParams
}

func (fnqr *fnQAWSR) Evaluate(x float64) float64 {
	f := fnqr.p.f
	t := fnqr.p.table

	factor := 1.0

	if t.beta != 0.0 {
		factor *= math.Pow(fnqr.p.b-x, t.beta)
	}

	if t.nu == 1 {
		factor *= math.Log(fnqr.p.b - x)
	}

	return factor * f.Evaluate(x)
}

type fnQAWSL struct {
	p *fnQAWSParams
}

func (fnql *fnQAWSL) Evaluate(x float64) float64 {
	f := fnql.p.f
	t := fnql.p.table

	factor := 1.0

	if t.alpha != 0.0 {
		factor *= math.Pow(x-fnql.p.a, t.alpha)
	}

	if t.mu == 1 {
		factor *= math.Log(x - fnql.p.a)
	}

	return factor * f.Evaluate(x)
}

type fnQAWS struct {
	p *fnQAWSParams
}

func (fnq *fnQAWS) Evaluate(x float64) float64 {
	f := fnq.p.f
	t := fnq.p.table

	factor := 1.0

	if t.alpha != 0.0 {
		factor *= math.Pow(x-fnq.p.a, t.alpha)
	}

	if t.beta != 0.0 {
		factor *= math.Pow(fnq.p.b-x, t.beta)
	}
	if t.mu == 1 {
		factor *= math.Log(x - fnq.p.a)
	}
	if t.nu == 1 {
		factor *= math.Log(fnq.p.b - x)
	}
	return factor * f.Evaluate(x)
}

func qc25s(f gsl.Function, a, b, a1, b1 float64, t *QAWSTable, result, abserr *float64, err_reliable *int) {
	params := &fnQAWSParams{f, a, b, t}

	if a1 == a && (t.alpha != 0.0 || t.mu != 0) {
		cheb12 := make([]float64, 13)
		cheb24 := make([]float64, 25)

		factor := math.Pow(0.5*(b1-a1), t.alpha+1.0)

		weighted_function := &fnQAWSR{params}
		Qcheb(weighted_function, a1, b1, cheb12, cheb24)

		if t.mu == 0 {
			res12 := 0.
			res24 := 0.
			u := factor

			computeResult(t.ri, cheb12, cheb24, &res12, &res24)

			*result = u * res24
			*abserr = math.Abs(u * (res24 - res12))
		} else {
			var res12a, res24a float64
			var res12b, res24b float64

			u := factor * math.Log(b1-a1)
			v := factor

			computeResult(t.ri, cheb12, cheb24, &res12a, &res24a)
			computeResult(t.rg, cheb12, cheb24, &res12b, &res24b)

			*result = u*res24a + v*res24b
			*abserr = math.Abs(u*(res24a-res12a)) + math.Abs(v*(res24b-res12b))
		}

		*err_reliable = 0

		return
	} else if b1 == b && (t.beta != 0.0 || t.nu != 0) {
		cheb12 := make([]float64, 13)
		cheb24 := make([]float64, 25)
		factor := math.Pow(0.5*(b1-a1), t.beta+1.0)

		weighted_function := &fnQAWSL{params}

		Qcheb(weighted_function, a1, b1, cheb12, cheb24)

		if t.nu == 0 {
			var res12, res24 float64
			u := factor

			computeResult(t.rj, cheb12, cheb24, &res12, &res24)

			*result = u * res24
			*abserr = math.Abs(u * (res24 - res12))
		} else {
			var res12a, res24a float64
			var res12b, res24b float64

			u := factor * math.Log(b1-a1)
			v := factor

			computeResult(t.rj, cheb12, cheb24, &res12a, &res24a)
			computeResult(t.rh, cheb12, cheb24, &res12b, &res24b)

			*result = u*res24a + v*res24b
			*abserr = math.Abs(u*(res24a-res12a)) + math.Abs(v*(res24b-res12b))
		}

		*err_reliable = 0

		return
	} else {
		var resabs, resasc float64

		weighted_function := &fnQAWS{params}

		Qk15(weighted_function, a1, b1, result, abserr, &resabs, &resasc)

		if *abserr == resasc {
			*err_reliable = 0
		} else {
			*err_reliable = 1
		}

		return
	}

}

func computeResult(r, cheb12, cheb24 []float64, result12, result24 *float64) {

	res12 := 0.
	res24 := 0.

	for i := 0; i < 13; i++ {
		res12 += r[i] * cheb12[i]
	}

	for i := 0; i < 25; i++ {
		res24 += r[i] * cheb24[i]
	}

	*result12 = res12
	*result24 = res24
}
