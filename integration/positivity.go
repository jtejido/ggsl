/* integration/positivity.c
 *
 * Copyright (C) 2017 Konrad Griessinger, Patrick Alken
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

/* Compare the integral of f(x) with the integral of |f(x)|
   to determine if f(x) covers both positive and negative values */
func test_positivity(result, resabs float64) bool {
	return (math.Abs(result) >= (1-50*gsl.Float64Eps)*resabs)
}
