/* sys/prec.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
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
package sys

import (
	gsl "github.com/jtejido/ggsl"
)

var (
	PrecEps = [gsl.PREC_T_NUM]float64{
		gsl.Float64Eps,
		gsl.Float32Eps,
		gsl.SfltEps,
	}

	PrecSqrtEps = [gsl.PREC_T_NUM]float64{
		gsl.SqrtFloat64Eps,
		gsl.SqrtFloat32Eps,
		gsl.SqrtSfltEps,
	}

	PrecRoot3Eps = [gsl.PREC_T_NUM]float64{
		gsl.Root3Float64Eps,
		gsl.Root3Float32Eps,
		gsl.Root3SfltEps,
	}

	PrecRoot4Eps = [gsl.PREC_T_NUM]float64{
		gsl.Root4Float64Eps,
		gsl.Root4Float32Eps,
		gsl.Root4SfltEps,
	}

	PrecRoot5Eps = [gsl.PREC_T_NUM]float64{
		gsl.Root5Float64Eps,
		gsl.Root5Float32Eps,
		gsl.Root5SfltEps,
	}

	PrecRoot6Eps = [gsl.PREC_T_NUM]float64{
		gsl.Root6Float64Eps,
		gsl.Root6Float32Eps,
		gsl.Root6SfltEps,
	}
)
