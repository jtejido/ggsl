/* integration/qmomo.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
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
	"github.com/jtejido/ggsl/err"
	"math"
)

type QAWSTable struct {
	alpha, beta    float64
	mu, nu         int
	ri, rj, rg, rh []float64
}

func NewQAWSTable(alpha, beta float64, mu, nu int) (*QAWSTable, err.GSLError) {
	t := new(QAWSTable)

	if alpha < -1.0 {
		return nil, err.ERROR("alpha must be greater than -1.0", err.EINVAL)
	}

	if beta < -1.0 {
		return nil, err.ERROR("beta must be greater than -1.0", err.EINVAL)
	}

	if mu != 0 && mu != 1 {
		return nil, err.ERROR("mu must be 0 or 1", err.EINVAL)
	}

	if nu != 0 && nu != 1 {
		return nil, err.ERROR("nu must be 0 or 1", err.EINVAL)
	}

	t.alpha = alpha
	t.beta = beta
	t.mu = mu
	t.nu = nu
	t.ri = make([]float64, 25)
	t.rj = make([]float64, 25)
	t.rg = make([]float64, 25)
	t.rh = make([]float64, 25)

	t.Initialize()

	return t, nil
}

func (t *QAWSTable) Initialize() {
	alpha_p1 := t.alpha + 1.0
	beta_p1 := t.beta + 1.0

	alpha_p2 := t.alpha + 2.0
	beta_p2 := t.beta + 2.0

	r_alpha := math.Pow(2.0, alpha_p1)
	r_beta := math.Pow(2.0, beta_p1)

	var an, anm1 float64

	t.ri[0] = r_alpha / alpha_p1
	t.ri[1] = t.ri[0] * t.alpha / alpha_p2

	an = 2.0
	anm1 = 1.0

	for i := 2; i < 25; i++ {
		t.ri[i] = -(r_alpha + an*(an-alpha_p2)*t.ri[i-1]) / (anm1 * (an + alpha_p1))
		anm1 = an
		an = an + 1.0
	}

	t.rj[0] = r_beta / beta_p1
	t.rj[1] = t.rj[0] * t.beta / beta_p2

	an = 2.0
	anm1 = 1.0

	for i := 2; i < 25; i++ {
		t.rj[i] = -(r_beta + an*(an-beta_p2)*t.rj[i-1]) / (anm1 * (an + beta_p1))
		anm1 = an
		an = an + 1.0
	}

	t.rg[0] = -t.ri[0] / alpha_p1
	t.rg[1] = -t.rg[0] - 2.0*r_alpha/(alpha_p2*alpha_p2)

	an = 2.0
	anm1 = 1.0

	for i := 2; i < 25; i++ {
		t.rg[i] = -(an*(an-alpha_p2)*t.rg[i-1] - an*t.ri[i-1] + anm1*t.ri[i]) / (anm1 * (an + alpha_p1))
		anm1 = an
		an = an + 1.0
	}

	t.rh[0] = -t.rj[0] / beta_p1
	t.rh[1] = -t.rh[0] - 2.0*r_beta/(beta_p2*beta_p2)

	an = 2.0
	anm1 = 1.0

	for i := 2; i < 25; i++ {
		t.rh[i] = -(an*(an-beta_p2)*t.rh[i-1] - an*t.rj[i-1] + anm1*t.rj[i]) / (anm1 * (an + beta_p1))
		anm1 = an
		an = an + 1.0
	}

	for i := 1; i < 25; i += 2 {
		t.rj[i] *= -1
		t.rh[i] *= -1
	}

}

func (t *QAWSTable) Set(alpha, beta float64, mu, nu int) err.GSLError {
	if alpha < -1.0 {
		return err.ERROR("alpha must be greater than -1.0", err.EINVAL)
	}

	if beta < -1.0 {
		return err.ERROR("beta must be greater than -1.0", err.EINVAL)
	}

	if mu != 0 && mu != 1 {
		return err.ERROR("mu must be 0 or 1", err.EINVAL)
	}

	if nu != 0 && nu != 1 {
		return err.ERROR("nu must be 0 or 1", err.EINVAL)
	}

	t.alpha = alpha
	t.beta = beta
	t.mu = mu
	t.nu = nu

	t.Initialize()
	return nil
}
