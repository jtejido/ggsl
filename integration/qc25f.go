/* integration/qc25f.c
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

type fnSin struct {
	f     gsl.Function
	omega float64
}

func (fs *fnSin) Evaluate(x float64) float64 {
	w := fs.omega
	wx := w * x
	sinwx := math.Sin(wx)
	return fs.f.Evaluate(x) * sinwx
}

type fnCos struct {
	f     gsl.Function
	omega float64
}

func (fc *fnCos) Evaluate(x float64) float64 {
	w := fc.omega
	wx := w * x
	coswx := math.Cos(wx)
	return fc.f.Evaluate(x) * coswx
}

func qc25f(f gsl.Function, a, b float64, wf *QAWOTable, level int, result, abserr, resabs, resasc *float64) {
	center := 0.5 * (a + b)
	half_length := 0.5 * (b - a)
	omega := wf.omega
	par := omega * half_length

	if math.Abs(par) < 2 {
		var weighted_function gsl.Function

		if wf.sine == GSL_INTEG_SINE {
			weighted_function = &fnSin{f, omega}
		} else {
			weighted_function = &fnCos{f, omega}
		}

		Qk15(weighted_function, a, b, result, abserr, resabs, resasc)

		return
	} else {
		cheb12 := make([]float64, 13)
		cheb24 := make([]float64, 25)
		var (
			result_abs, res12_cos, res12_sin, res24_cos, res24_sin float64
			est_cos, est_sin                                       float64
			c, s                                                   float64
		)

		Qcheb(f, a, b, cheb12, cheb24)

		if level >= wf.n {
			/* table overflow should not happen, check before calling */
			err.ERROR_VOID("table overflow in internal function", err.ESANITY)
		}

		/* obtain moments from the table */

		moment := &slice{wf.chebmo.data, 25 * level}

		res12_cos = cheb12[12] * moment.data[moment.offset+12]
		res12_sin = 0

		for i := 0; i < 6; i++ {
			k := 10 - 2*i
			res12_cos += cheb12[k] * moment.data[moment.offset+k]
			res12_sin += cheb12[k+1] * moment.data[moment.offset+k+1]
		}

		res24_cos = cheb24[24] * moment.data[moment.offset+24]
		res24_sin = 0

		result_abs = math.Abs(cheb24[24])

		for i := 0; i < 12; i++ {
			k := 22 - 2*i
			res24_cos += cheb24[k] * moment.data[moment.offset+k]
			res24_sin += cheb24[k+1] * moment.data[moment.offset+k+1]
			result_abs += math.Abs(cheb24[k]) + math.Abs(cheb24[k+1])
		}

		est_cos = math.Abs(res24_cos - res12_cos)
		est_sin = math.Abs(res24_sin - res12_sin)

		c = half_length * math.Cos(center*omega)
		s = half_length * math.Sin(center*omega)

		if wf.sine == GSL_INTEG_SINE {
			*result = c*res24_sin + s*res24_cos
			*abserr = math.Abs(c*est_sin) + math.Abs(s*est_cos)
		} else {
			*result = c*res24_cos - s*res24_sin
			*abserr = math.Abs(c*est_cos) + math.Abs(s*est_sin)
		}

		*resabs = result_abs * half_length
		*resasc = gsl.MaxFloat64

		return
	}
}
