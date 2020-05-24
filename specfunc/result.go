/* specfunc/result.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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
package specfunc

import (
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"math"
)

type Result struct {
	val float64
	err float64
}

func (r *Result) Value() float64 {
	return r.val
}

func (r *Result) Err() float64 {
	return r.err
}

type Result_e10 struct {
	val float64
	err float64
	e10 int
}

func (r *Result_e10) Value() float64 {
	return r.val
}

func (r *Result_e10) Err() float64 {
	return r.err
}

func (r *Result_e10) E10() int {
	return r.e10
}

func Result_smash_e(re *Result_e10, r *Result) err.GSLError {
	if re.e10 == 0 {
		/* nothing to smash */
		r.val = re.val
		r.err = re.err
		return nil
	} else {
		av := math.Abs(re.val)
		ae := math.Abs(re.err)

		if gsl.SqrtMinFloat64 < av && av < gsl.SqrtMaxFloat64 && gsl.SqrtMinFloat64 < ae && ae < gsl.SqrtMaxFloat64 && 0.49*gsl.LnMinFloat64 < float64(re.e10) && float64(re.e10) < 0.49*gsl.LnMaxFloat64 {
			scale := math.Exp(float64(re.e10) * math.Ln10)
			r.val = re.val * scale
			r.err = re.err * scale
			return nil
		} else {
			return Exp_mult_err_e(float64(re.e10)*math.Ln10, 0.0, re.val, re.err, r)
		}
	}
}
