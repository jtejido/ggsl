/* integration/integration.c
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
	"fmt"
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"github.com/lucky-se7en/ggsl/sys"
	"math"
	"os"
)

type QKFunction func(f gsl.Function, a, b float64, result, abserr, defabs, resabs *float64)

/****************** TESTING **************************/
const (
	sqrt15 = 3.8729833462074168852
	sqrt30 = 5.4772255750516611346
	sqrt70 = 8.3666002653407554798
	const1 = 0.86113631159405257522 /* sqrt((3+2*sqrt(6./5))/7) */
	const2 = 0.33998104358485626480 /* sqrt((3-2*sqrt(6./5))/7) */
	const3 = 0.90617984593866399280 /* sqrt((5+2*sqrt(10./7)))/3 */
	const4 = 0.53846931010568309104 /* sqrt((5-2*sqrt(10./7)))/3 */
)

func initTesting() {
	err.SetErrorHandler(myErrorHandler) // to avoid hitting the default handler
	os.Setenv("GSL_TEST_VERBOSE", "1")  // to print results
}

func myErrorHandler(reason, file string, line, gsl_errno int) {
	fmt.Printf("(caught [%s:%d: %s (%d)])\n", file, line, reason, gsl_errno)
}

type counter struct {
	f     gsl.Function
	neval int
}

func (c *counter) Evaluate(x float64) float64 {
	c.neval++
	return c.f.Evaluate(x)
}

type float struct {
	value     float64
	precision float64
}

type f1 struct {
	alpha float64
}

func (f f1) Evaluate(x float64) float64 {
	return math.Pow(x, f.alpha) * math.Log(1.0/x)
}

type f3 struct {
	alpha float64
}

func (f f3) Evaluate(x float64) float64 {
	return math.Cos(math.Pow(2.0, f.alpha) * math.Sin(x))
}

type f11 struct {
	alpha float64
}

func (f f11) Evaluate(x float64) float64 {
	return math.Pow(math.Log(1/x), f.alpha-1)
}

type f15 struct {
	alpha float64
}

func (f f15) Evaluate(x float64) float64 {
	return x * x * math.Exp(-math.Pow(2.0, -f.alpha)*x)
}

type f16 struct {
	alpha float64
}

func (f f16) Evaluate(x float64) float64 {
	if x == 0 && f.alpha == 1 {
		return 1
	}
	if x == 0 && f.alpha > 1 {
		return 0
	}
	return math.Pow(x, f.alpha-1) / math.Pow(1+10.0*x, 2.0)
}

type f454 struct {
}

func (f f454) Evaluate(x float64) float64 {
	x2 := x * x
	x3 := x * x2
	return x3 * math.Log(math.Abs((x2-1.0)*(x2-2.0)))
}

type f455 struct {
}

func (f f455) Evaluate(x float64) float64 {
	return math.Log(x) / (1.0 + 100.0*x*x)
}

type f456 struct {
}

func (f f456) Evaluate(x float64) float64 {
	if x == 0.0 {
		return 0
	}
	return math.Log(x)
}

type f457 struct {
}

func (f f457) Evaluate(x float64) float64 {
	if x == 0.0 {
		return 0
	}
	return 1 / math.Sqrt(x)
}

type f458 struct {
}

func (f f458) Evaluate(x float64) float64 {

	if x == 0.0 {
		return 0
	}

	u := math.Log(x)
	v := 1 + u*u

	return 1.0 / (v * v)
}

type f459 struct {
}

func (f f459) Evaluate(x float64) float64 {
	return 1.0 / (5.0*x*x*x + 6.0)
}

type myfn1 struct {
}

func (f myfn1) Evaluate(x float64) float64 {
	return math.Exp(-x - x*x)
}

type myfn2 struct {
	alpha float64
}

func (f myfn2) Evaluate(x float64) float64 {
	return math.Exp(f.alpha * x)
}

type f_sin struct {
}

func (f f_sin) Evaluate(x float64) float64 {
	return math.Sin(x)
}

/* integ(f_sin,x,a,b) */
func integ_f_sin(a, b float64) float64 {
	return -math.Cos(b) + math.Cos(a)
}

type cqf11 struct {
}

func (f cqf11) Evaluate(x float64) float64 {
	return 1.0 / (1 + math.Exp(x))
}

type monomial_params struct {
	degree   int
	constant float64
}

type f_monomial struct {
	params *monomial_params
}

func (f f_monomial) Evaluate(x float64) float64 {
	return f.params.constant * sys.PowInt(x, f.params.degree)
}

func integ_f_monomial(a, b float64, p *monomial_params) float64 {
	degreep1 := p.degree + 1
	bnp1 := sys.PowInt(b, degreep1)
	anp1 := sys.PowInt(a, degreep1)
	return (p.constant / float64(degreep1)) * (bnp1 - anp1)
}

/* The test functions. */
type cqf1 struct{}

func (f cqf1) Evaluate(x float64) float64 {
	return math.Exp(x)
}

type cqf2 struct{}

func (f cqf2) Evaluate(x float64) float64 {
	if x >= 0.3 {
		return 1
	}

	return 0
}

type cqf3 struct{}

func (f cqf3) Evaluate(x float64) float64 {
	return math.Sqrt(x)
}

type cqf4 struct{}

func (f cqf4) Evaluate(x float64) float64 {
	return (23.0/25)*math.Cosh(x) - math.Cos(x)
}

type cqf5 struct{}

func (f cqf5) Evaluate(x float64) float64 {
	x2 := x * x
	return 1.0 / (x2*(x2+1) + 0.9)
}

type cqf6 struct{}

func (f cqf6) Evaluate(x float64) float64 {
	return x * math.Sqrt(x)
}

type cqf7 struct{}

func (f cqf7) Evaluate(x float64) float64 {
	x2 := x * x
	return 1.0 / (1 + x2*x2)
}

type cqf8 struct{}

func (f cqf8) Evaluate(x float64) float64 {
	x2 := x * x
	return 1.0 / (1 + x2*x2)
}

type cqf9 struct{}

func (f cqf9) Evaluate(x float64) float64 {
	return 2.0 / (2 + math.Sin(10*gsl.Pi*x))
}

type cqf10 struct{}

func (f cqf10) Evaluate(x float64) float64 {
	return 1.0 / (1 + x)
}

// type cqf11 struct{}

// func (f cqf11) Evaluate(x float64) float64 {
// 	return 1.0 / (1 + math.Exp(x))
// }

type cqf12 struct{}

func (f cqf12) Evaluate(x float64) float64 {
	return x / (math.Exp(x) - 1.0)
}

type cqf13 struct{}

func (f cqf13) Evaluate(x float64) float64 {
	return math.Sin(100*gsl.Pi*x) / (gsl.Pi * x)
}

type cqf14 struct{}

func (f cqf14) Evaluate(x float64) float64 {
	return math.Sqrt(50.0) * math.Exp(-50*gsl.Pi*x*x)
}

type cqf15 struct{}

func (f cqf15) Evaluate(x float64) float64 {
	return 25.0 * math.Exp(-25*x)
}

type cqf16 struct{}

func (f cqf16) Evaluate(x float64) float64 {
	return 50 / gsl.Pi * (2500*x*x + 1)
}

type cqf17 struct{}

func (f cqf17) Evaluate(x float64) float64 {
	t1 := 50 * gsl.Pi * x
	t2 := math.Sin(t1) / t1
	return 50 * t2 * t2
}

type cqf18 struct{}

func (f cqf18) Evaluate(x float64) float64 {
	return math.Cos(math.Cos(x) + 3*math.Sin(x) + 2*math.Cos(2*x) + 3*math.Sin(2*x) + 3*math.Cos(3*x))
}

type cqf19 struct{}

func (f cqf19) Evaluate(x float64) float64 {
	return math.Log(x)
}

type cqf20 struct{}

func (f cqf20) Evaluate(x float64) float64 {
	return 1 / (x*x + 1.005)
}

type cqf21 struct{}

func (f cqf21) Evaluate(x float64) float64 {
	return 1/math.Cosh(10*(x-0.2)*2) + 1/math.Cosh(100*(x-0.4)*4) + 1/math.Cosh(1000*(x-0.6)*8)
}

type cqf22 struct{}

func (f cqf22) Evaluate(x float64) float64 {
	return 4 * gsl.Pi * gsl.Pi * x * math.Sin(20*gsl.Pi*x) * math.Cos(2*gsl.Pi*x)
}

type cqf23 struct{}

func (f cqf23) Evaluate(x float64) float64 {
	t := 230*x - 30
	return 1 / (1 + t*t)
}

type cqf24 struct{}

func (f cqf24) Evaluate(x float64) float64 {
	return math.Floor(math.Exp(x))
}

type cqf25 struct{}

func (f cqf25) Evaluate(x float64) float64 {
	var s1, s2, s3 float64
	if x < 1 {
		s1 = 1
	}
	if 1 <= x && x <= 3 {
		s2 = 1
	}
	if x > 3 {
		s3 = 1
	}
	return s1*(x+1) + s2*(3-x) + s3*2
}
