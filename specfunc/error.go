/* specfunc/error.c
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
	"runtime"
)

const (
	domain_err    = "domain error"
	underflow_err = "underflow"
	overflow_err  = "overflow"
)

func OverflowError(result ...*Result) err.GSLError {
	for _, r := range result {
		r.val = math.Inf(1)
		r.err = math.Inf(1)
	}
	errno := err.EOVRFLW
	reason := overflow_err
	_, file, line, _ := runtime.Caller(1)
	err.Error(reason, file, line, errno)
	return err.New(errno, reason)
}

func UnderflowError(result ...*Result) err.GSLError {
	for _, r := range result {
		r.val = 0.0
		r.err = gsl.MinFloat64
	}
	errno := err.EUNDRFLW
	reason := underflow_err
	_, file, line, _ := runtime.Caller(1)
	err.Error(reason, file, line, errno)
	return err.New(errno, reason)
}

func InternalOverflowError(result *Result) err.GSLError {
	result.val = math.Inf(1)
	result.err = math.Inf(1)
	return err.Overflow()
}

func InternalUnderflowError(result *Result) err.GSLError {
	result.val = 0.0
	result.err = gsl.MinFloat64
	return err.Underflow()
}

func DomainError(result ...*Result) err.GSLError {
	for _, r := range result {
		r.val = math.NaN()
		r.err = math.NaN()
	}
	errno := err.EDOM
	reason := domain_err
	_, file, line, _ := runtime.Caller(1)
	err.Error(reason, file, line, errno)
	return err.New(errno, reason)
}

func DomainErrorMsg(msg string, result *Result) err.GSLError {
	result.val = math.NaN()
	result.err = math.NaN()
	errno := err.EDOM
	_, file, line, _ := runtime.Caller(1)
	err.Error(msg, file, line, errno)
	return err.New(errno, msg)
}

func DomainError_e10(result *Result_e10) err.GSLError {
	result.val = math.NaN()
	result.err = math.NaN()
	result.e10 = 0
	errno := err.EDOM
	reason := domain_err
	_, file, line, _ := runtime.Caller(1)
	err.Error(reason, file, line, errno)
	return err.New(errno, reason)
}

func OverflowError_e10(result *Result_e10) err.GSLError {
	result.val = math.Inf(1)
	result.err = math.Inf(1)
	result.e10 = 0
	errno := err.EOVRFLW
	reason := overflow_err
	_, file, line, _ := runtime.Caller(1)
	err.Error(reason, file, line, errno)
	return err.New(errno, reason)
}

func UnderflowError_e10(result *Result_e10) err.GSLError {
	result.val = 0.0
	result.err = gsl.MinFloat64
	result.e10 = 0
	errno := err.EUNDRFLW
	reason := underflow_err
	_, file, line, _ := runtime.Caller(1)
	err.Error(reason, file, line, errno)
	return err.New(errno, reason)
}

func MaxIterError(result *Result) err.GSLError {
	result.val = math.NaN()
	result.err = math.NaN()
	errno := err.EMAXITER
	reason := "too many iterations error"
	_, file, line, _ := runtime.Caller(1)
	err.Error(reason, file, line, errno)
	return err.New(errno, reason)
}
