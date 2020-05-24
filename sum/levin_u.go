/* sum/levin_u.c
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
package sum

import (
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"math"
)

func SumLevinUAccel(array []float64, array_size int, w *SumLevinUWorkspace, sum_accel, abserr *float64) err.GSLError {
	return SumLevinUMinMax(array, array_size, 0, array_size-1, w, sum_accel, abserr)
}

func SumLevinUMinMax(array []float64, array_size, min_terms, max_terms int, w *SumLevinUWorkspace, sum_accel, abserr *float64) err.GSLError {
	/* Ignore any trailing zeros in the array */
	size := array_size

	for size > 0 && array[size-1] == 0 {
		size--
	}

	if size == 0 {
		*sum_accel = 0.0
		*abserr = 0.0
		w.sum_plain = 0.0
		w.terms_used = 0
		return nil
	} else if size == 1 {
		*sum_accel = array[0]
		*abserr = 0.0
		w.sum_plain = array[0]
		w.terms_used = 1
		return nil
	} else {
		var (
			SMALL = 0.01
			nmax  = int(math.Max(float64(max_terms), float64(array_size))) - 1
		)
		noise_n := 0.0
		// noise_nm1 := 0.0
		trunc_n := 0.0
		trunc_nm1 := 0.0
		actual_trunc_n := 0.0
		actual_trunc_nm1 := 0.0
		result_n := 0.0
		result_nm1 := 0.0
		variance := 0.0
		n := 0
		i := 0
		better := false
		before := false
		converging := false
		least_trunc := gsl.MaxFloat64
		least_trunc_noise := gsl.MaxFloat64
		least_trunc_result := 0.0

		/* Calculate specified minimum number of terms.  No convergence
		   tests are made, and no truncation information is stored.  */

		for n = 0; n < min_terms; n++ {
			t := array[n]
			result_nm1 = result_n
			SumLevinUStep(t, n, nmax, w, &result_n)
		}

		least_trunc_result = result_n

		for i = 0; i < n; i++ {
			dn := w.dsum[i] * gsl.MachEps * array[i]
			variance += dn * dn
		}

		noise_n = math.Sqrt(variance)

		/* Calculate up to maximum number of terms.  Check truncation
		   condition.  */

		for ; n <= nmax; n++ {
			t := array[n]

			result_nm1 = result_n
			SumLevinUStep(t, n, nmax, w, &result_n)

			/* Compute the truncation error directly */

			actual_trunc_nm1 = actual_trunc_n
			actual_trunc_n = math.Abs(result_n - result_nm1)

			/* Average results to make a more reliable estimate of the
			   real truncation error */

			trunc_nm1 = trunc_n
			trunc_n = 0.5 * (actual_trunc_n + actual_trunc_nm1)

			// noise_nm1 = noise_n
			variance = 0.0

			for i = 0; i <= n; i++ {
				dn := w.dsum[i] * gsl.MachEps * array[i]
				variance += dn * dn
			}

			noise_n = math.Sqrt(variance)

			/* Determine if we are in the convergence region.  */

			better = (trunc_n < trunc_nm1 || trunc_n < SMALL*math.Abs(result_n))
			converging = converging || (better && before)
			before = better

			if converging {
				if trunc_n < least_trunc {
					/* Found a low truncation point in the convergence
					   region. Save it. */

					least_trunc_result = result_n
					least_trunc = trunc_n
					least_trunc_noise = noise_n
				}

				if noise_n > trunc_n/3.0 {
					break
				}

				if trunc_n < 10.0*gsl.MachEps*math.Abs(result_n) {
					break
				}

			}

		}

		if converging {
			/* Stopped in the convergence region.  Return result and
			   error estimate.  */

			*sum_accel = least_trunc_result
			*abserr = math.Max(least_trunc, least_trunc_noise)
			w.terms_used = n
			return nil
		} else {
			/* Never reached the convergence region.  Use the last
			   calculated values.  */

			*sum_accel = result_n
			*abserr = math.Max(trunc_n, noise_n)
			w.terms_used = n
			return nil
		}
	}
}

func SumLevinUStep(term float64, n, nmax int, w *SumLevinUWorkspace, sum_accel *float64) err.GSLError {

	I := func(ii, jj int) int {
		return (ii)*(nmax+1) + (jj)
	}

	if n == 0 {
		*sum_accel = term
		w.sum_plain = term
		w.q_den[0] = 1.0 / term
		w.q_num[0] = 1.0
		w.dq_den[I(0, 0)] = -1.0 / (term * term)
		w.dq_num[I(0, 0)] = 0.0
		w.dsum[0] = 1.0
		return nil
	} else {
		result := 0.0
		factor := 1.0
		ratio := float64(n / (n + 1))
		var i, j int

		w.sum_plain += term
		w.q_den[n] = 1.0 / (term * float64(n+1) * float64(n+1))
		w.q_num[n] = w.sum_plain * w.q_den[n]

		for i = 0; i < n; i++ {
			w.dq_den[I(i, n)] = 0
			w.dq_num[I(i, n)] = w.q_den[n]
		}

		w.dq_den[I(n, n)] = -w.q_den[n] / term
		w.dq_num[I(n, n)] = (w.q_den[n] + w.sum_plain) * w.dq_den[I(n, n)]

		for j = n - 1; j >= 0; j-- {
			c := factor * float64(j+1) / float64(n+1)
			factor *= ratio
			w.q_den[j] = w.q_den[j+1] - c*w.q_den[j]
			w.q_num[j] = w.q_num[j+1] - c*w.q_num[j]

			for i = 0; i < n; i++ {
				w.dq_den[I(i, j)] = w.dq_den[I(i, j+1)] - c*w.dq_den[I(i, j)]
				w.dq_num[I(i, j)] = w.dq_num[I(i, j+1)] - c*w.dq_num[I(i, j)]
			}

			w.dq_den[I(n, j)] = w.dq_den[I(n, j+1)]
			w.dq_num[I(n, j)] = w.dq_num[I(n, j+1)]
		}

		result = w.q_num[0] / w.q_den[0]

		*sum_accel = result

		for i = 0; i <= n; i++ {
			w.dsum[i] = (w.dq_num[I(i, 0)] - result*w.dq_den[I(i, 0)]) / w.q_den[0]
		}

		return nil
	}
}
