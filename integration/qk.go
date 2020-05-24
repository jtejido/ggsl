/* integration/qk.c
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

func Qk(n int, xgk, wg, wgk, fv1, fv2 []float64, f gsl.Function, a, b float64, result, abserr, resabs, resasc *float64) {
	center := 0.5 * (a + b)
	half_length := 0.5 * (b - a)
	abs_half_length := math.Abs(half_length)
	f_center := f.Evaluate(center)

	result_gauss := 0.0
	result_kronrod := f_center * wgk[n-1]

	result_abs := math.Abs(result_kronrod)
	result_asc := 0.0
	mean := 0.0
	err := 0.0

	if n%2 == 0 {
		result_gauss = f_center * wg[n/2-1]
	}

	for j := 0; j < (n-1)/2; j++ {
		jtw := j*2 + 1 /* in original fortran j=1,2,3 jtw=2,4,6 */
		abscissa := half_length * xgk[jtw]
		fval1 := f.Evaluate(center - abscissa)
		fval2 := f.Evaluate(center + abscissa)
		fsum := fval1 + fval2
		fv1[jtw] = fval1
		fv2[jtw] = fval2
		result_gauss += wg[j] * fsum
		result_kronrod += wgk[jtw] * fsum
		result_abs += wgk[jtw] * (math.Abs(fval1) + math.Abs(fval2))
	}

	for j := 0; j < n/2; j++ {
		jtwm1 := j * 2
		abscissa := half_length * xgk[jtwm1]
		fval1 := f.Evaluate(center - abscissa)
		fval2 := f.Evaluate(center + abscissa)
		fv1[jtwm1] = fval1
		fv2[jtwm1] = fval2
		result_kronrod += wgk[jtwm1] * (fval1 + fval2)
		result_abs += wgk[jtwm1] * (math.Abs(fval1) + math.Abs(fval2))
	}

	mean = result_kronrod * 0.5

	result_asc = wgk[n-1] * math.Abs(f_center-mean)

	for j := 0; j < n-1; j++ {
		result_asc += wgk[j] * (math.Abs(fv1[j]-mean) + math.Abs(fv2[j]-mean))
	}

	/* scale by the width of the integration region */

	err = (result_kronrod - result_gauss) * half_length

	result_kronrod *= half_length
	result_abs *= abs_half_length
	result_asc *= abs_half_length

	*result = result_kronrod
	*resabs = result_abs
	*resasc = result_asc
	*abserr = RescaleError(err, result_abs, result_asc)
}
