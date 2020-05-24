/* integration/err.c
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

func RescaleError(err, result_abs, result_asc float64) float64 {
	err = math.Abs(err)

	if result_asc != 0 && err != 0 {
		scale := math.Pow((200 * err / result_asc), 1.5)

		if scale < 1 {
			err = result_asc * scale
		} else {
			err = result_asc
		}
	}

	if result_abs > gsl.MinFloat64/(50*gsl.Float64Eps) {
		min_err := 50 * gsl.Float64Eps * result_abs

		if min_err > err {
			err = min_err
		}
	}

	return err
}
