/* sum/work_u.c
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
	"github.com/lucky-se7en/ggsl/err"
)

type SumLevinUWorkspace struct {
	size                               int
	i                                  int /* position in array */
	terms_used                         int /* number of calls */
	sum_plain                          float64
	q_num, q_den, dq_num, dq_den, dsum []float64
}

func NewSumLevinUWorkspace(n int) (*SumLevinUWorkspace, err.GSLError) {
	if n == 0 {
		return nil, err.ERROR("length n must be positive integer", err.EDOM)
	}

	w := &SumLevinUWorkspace{
		q_num:  make([]float64, n),
		q_den:  make([]float64, n),
		dq_num: make([]float64, n*n),
		dq_den: make([]float64, n*n),
		dsum:   make([]float64, n),
		size:   n,
	}

	return w, nil
}

func (s *SumLevinUWorkspace) SumPlain() float64 {
	return s.sum_plain
}

func (s *SumLevinUWorkspace) TermsUsed() int {
	return s.terms_used
}
