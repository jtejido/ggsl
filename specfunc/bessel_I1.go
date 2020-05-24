/* specfunc/bessel_I1.c
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

const (
	ROOT_EIGHT = 2.0 * math.Sqrt2
)

var (
	bi1d = []float64{
		-0.001971713261099859,
		0.407348876675464810,
		0.034838994299959456,
		0.001545394556300123,
		0.000041888521098377,
		0.000000764902676483,
		0.000000010042493924,
		0.000000000099322077,
		0.000000000000766380,
		0.000000000000004741,
		0.000000000000000024,
	}

	bi1 = &chebyshevSeries{
		bi1d,
		10,
		-1, 1,
		10,
	}

	ai1d = []float64{
		-0.02846744181881479,
		-0.01922953231443221,
		-0.00061151858579437,
		-0.00002069971253350,
		0.00000858561914581,
		0.00000104949824671,
		-0.00000029183389184,
		-0.00000001559378146,
		0.00000001318012367,
		-0.00000000144842341,
		-0.00000000029085122,
		0.00000000012663889,
		-0.00000000001664947,
		-0.00000000000166665,
		0.00000000000124260,
		-0.00000000000027315,
		0.00000000000002023,
		0.00000000000000730,
		-0.00000000000000333,
		0.00000000000000071,
		-0.00000000000000006,
	}

	ai1 = &chebyshevSeries{
		ai1d,
		20,
		-1, 1,
		11,
	}

	ai12d = []float64{
		0.02857623501828014,
		-0.00976109749136147,
		-0.00011058893876263,
		-0.00000388256480887,
		-0.00000025122362377,
		-0.00000002631468847,
		-0.00000000383538039,
		-0.00000000055897433,
		-0.00000000001897495,
		0.00000000003252602,
		0.00000000001412580,
		0.00000000000203564,
		-0.00000000000071985,
		-0.00000000000040836,
		-0.00000000000002101,
		0.00000000000004273,
		0.00000000000001041,
		-0.00000000000000382,
		-0.00000000000000186,
		0.00000000000000033,
		0.00000000000000028,
		-0.00000000000000003,
	}

	ai12 = &chebyshevSeries{
		ai12d,
		21,
		-1, 1,
		9,
	}
)

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/
func Bessel_I1_scaled_e(x float64, result *Result) err.GSLError {
	xmin := 2.0 * gsl.MinFloat64
	x_small := ROOT_EIGHT * gsl.SqrtFloat64Eps
	y := math.Abs(x)

	if y == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if y < xmin {
		return UnderflowError(result)
	} else if y < x_small {
		result.val = 0.5 * x
		result.err = 0.0
		return nil
	} else if y <= 3.0 {
		ey := math.Exp(-y)
		c := new(Result)
		bi1.Evaluate(y*y/4.5-1.0, c)
		result.val = x * ey * (0.875 + c.val)
		result.err = ey*c.err + y*gsl.Float64Eps*math.Abs(result.val)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if y <= 8.0 {
		sy := math.Sqrt(y)
		c := new(Result)
		ai1.Evaluate((48.0/y-11.0)/5.0, c)
		b := (0.375 + c.val) / sy
		s := -1.
		if x > 0.0 {
			s = 1.
		}

		result.val = s * b
		result.err = c.err / sy
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}
	sy := math.Sqrt(y)
	c := new(Result)
	ai12.Evaluate(16.0/y-1.0, c)
	b := (0.375 + c.val) / sy
	s := -1.
	if x > 0.0 {
		s = 1.
	}

	result.val = s * b
	result.err = c.err / sy
	result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
	return nil

}

func Bessel_I1_e(x float64, result *Result) err.GSLError {
	xmin := 2.0 * gsl.MinFloat64
	x_small := ROOT_EIGHT * gsl.SqrtFloat64Eps
	y := math.Abs(x)

	if y == 0.0 {
		result.val = 0.0
		result.err = 0.0
		return nil
	} else if y < xmin {
		return UnderflowError(result)
	} else if y < x_small {
		result.val = 0.5 * x
		result.err = 0.0
		return nil
	} else if y <= 3.0 {
		c := new(Result)
		bi1.Evaluate(y*y/4.5-1.0, c)
		result.val = x * (0.875 + c.val)
		result.err = y * c.err
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if y < gsl.LnMaxFloat64 {
		ey := math.Exp(y)
		I1_scaled := new(Result)
		Bessel_I1_scaled_e(x, I1_scaled)
		result.val = ey * I1_scaled.val
		result.err = ey*I1_scaled.err + y*gsl.Float64Eps*math.Abs(result.val)
		result.err += 2.0 * gsl.Float64Eps * math.Abs(result.val)
		return nil
	}

	return OverflowError(result)

}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/
func Bessel_I1_scaled(x float64) float64 {
	result := new(Result)
	status := Bessel_I1_scaled_e(x, result)
	return EvalResult(result, status)
}

func Bessel_I1(x float64) float64 {
	result := new(Result)
	status := Bessel_I1_e(x, result)
	return EvalResult(result, status)
}
