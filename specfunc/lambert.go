/* specfunc/lambert.c
 *
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Gerard Jungman
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

/* Started with code donated by K. Briggs; added
 * error estimates, GSL foo, and minor tweaks.
 * Some Lambert-ology from
 *  [Corless, Gonnet, Hare, and Jeffrey, "On Lambert's W Function".]
 */

/* Halley iteration (eqn. 5.12, Corless et al) */
func halley_iteration(x, w_initial float64, max_iters uint, result *Result) err.GSLError {
	w := w_initial
	var i uint

	for i = 0; i < max_iters; i++ {
		var tol float64
		e := math.Exp(w)
		p := w + 1.0
		t := w*e - x
		/* printf("FOO: %20.16g  %20.16g\n", w, t); */

		if w > 0 {
			t = (t / p) / e /* Newton iteration */
		} else {
			t /= e*p - 0.5*(p+1.0)*t/p /* Halley iteration */
		}

		w -= t

		tol = 10 * gsl.Float64Eps * math.Max(math.Abs(w), 1.0/(math.Abs(p)*e))

		if math.Abs(t) < tol {
			result.val = w
			result.err = 2.0 * tol
			return nil
		}
	}

	/* should never get here */
	result.val = w
	result.err = math.Abs(w)
	return err.MaxIteration()
}

/* series which appears for q near zero;
 * only the argument is different for the different branches
 */
func series_eval(r float64) float64 {
	c := []float64{
		-1.0,
		2.331643981597124203363536062168,
		-1.812187885639363490240191647568,
		1.936631114492359755363277457668,
		-2.353551201881614516821543561516,
		3.066858901050631912893148922704,
		-4.175335600258177138854984177460,
		5.858023729874774148815053846119,
		-8.401032217523977370984161688514,
		12.250753501314460424,
		-18.100697012472442755,
		27.029044799010561650,
	}
	t_8 := c[8] + r*(c[9]+r*(c[10]+r*c[11]))
	t_5 := c[5] + r*(c[6]+r*(c[7]+r*t_8))
	t_1 := c[1] + r*(c[2]+r*(c[3]+r*(c[4]+r*t_5)))
	return c[0] + r*t_1
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Lambert_W0_e(x float64, result *Result) err.GSLError {
	one_over_E := 1.0 / math.E
	q := x + one_over_E

	if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if q < 0.0 {
		/* Strictly speaking this is an error. But because of the
		 * arithmetic operation connecting x and q, I am a little
		 * lenient in case of some epsilon overshoot. The following
		 * answer is quite accurate in that case. Anyway, we have
		 * to return GSL_EDOM.
		 */
		result.val = -1.0
		result.err = math.Sqrt(-q)
		return err.Domain()
	} else if q == 0.0 {
		result.val = -1.0
		result.err = gsl.Float64Eps /* cannot error is zero, maybe q == 0 by "accident" */
		return nil
	} else if q < 1.0e-03 {
		/* series near -1/E in sqrt(q) */
		r := math.Sqrt(q)
		result.val = series_eval(r)
		result.err = 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		var MAX_ITERS uint = 10
		var w float64

		if x < 1.0 {
			/* obtain initial approximation from series near x=0;
			 * no need for extra care, since the Halley iteration
			 * converges nicely on this branch
			 */
			p := math.Sqrt(2.0 * math.E * q)
			w = -1.0 + p*(1.0+p*(-1.0/3.0+p*11.0/72.0))
		} else {
			/* obtain initial approximation from rough asymptotic */
			w = math.Log(x)
			if x > 3.0 {
				w -= math.Log(w)
			}
		}

		return halley_iteration(x, w, MAX_ITERS, result)
	}
}

func Lambert_Wm1_e(x float64, result *Result) err.GSLError {
	if x > 0.0 {
		return Lambert_W0_e(x, result)
	} else if x == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else {
		var MAX_ITERS uint = 32
		one_over_E := 1.0 / math.E
		q := x + one_over_E
		var w float64

		if q < 0.0 {
			/* As in the W0 branch above, return some reasonable answer anyway. */
			result.val = -1.0
			result.err = math.Sqrt(-q)
			return err.Domain()
		}

		if x < -1.0e-6 {
			/* Obtain initial approximation from series about q = 0,
			 * as long as we're not very close to x = 0.
			 * Use full series and try to bail out if q is too small,
			 * since the Halley iteration has bad convergence properties
			 * in finite arithmetic for q very small, because the
			 * increment alternates and p is near zero.
			 */
			r := -math.Sqrt(q)
			w = series_eval(r)
			if q < 3.0e-3 {
				/* this approximation is good enough */
				result.val = w
				result.err = 5.0 * gsl.Float64Eps * math.Abs(w)
				return nil
			}
		} else {
			/* Obtain initial approximation from asymptotic near zero. */
			L_1 := math.Log(-x)
			L_2 := math.Log(-L_1)
			w = L_1 - L_2 + L_2/L_1
		}

		return halley_iteration(x, w, MAX_ITERS, result)
	}
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Lambert_W0(x float64) float64 {
	result := new(Result)
	status := Lambert_W0_e(x, result)
	return EvalResult(result, status)
}

func Lambert_Wm1(x float64) float64 {
	result := new(Result)
	status := Lambert_Wm1_e(x, result)
	return EvalResult(result, status)
}
