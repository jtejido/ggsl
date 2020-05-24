package integration

import (
	"fmt"
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/test"
	"math"
	"testing"
)

func TestGLFixed(t *testing.T) {
	initTesting()
	/* Test the fixed-order Gauss-Legendre rules with a monomial. */
	{

		a := 0.0
		b := 1.2
		params := new(monomial_params)
		params.constant = 1.0
		f := f_monomial{params}
		for n := 1; n < 1025; n++ {
			var expected, result float64

			tbl, _ := NewGLFixedTable(n)
			f.params.degree = 2*n - 1 /* n point rule exact for 2n-1 degree poly */
			expected = integ_f_monomial(a, b, params)
			result = GLFixed(&f, a, b, tbl)

			if tbl.precomputed == 1 {
				test.TestRel(t, result, expected, 1.0e-12, fmt.Sprintf("glfixed %d-point: Integrating (%g*x^%d) over [%g,%g]", n, params.constant, params.degree, a, b))
			} else {
				test.TestRel(t, result, expected, 1.0e-7, fmt.Sprintf("glfixed %d-point: Integrating (%g*x^%d) over [%g,%g]", n, params.constant, params.degree, a, b))
			}

		}
	}

	/* Sanity check sin(x) test function for fixed Gauss-Legendre rules */
	{
		f := &f_sin{}

		test.TestAbs(t, f.Evaluate(2.0), math.Sin(2.0), 0.0, "f_sin sanity check 1")
		test.TestAbs(t, f.Evaluate(7.0), math.Sin(7.0), 0.0, "f_sin sanity check 2")
		test.TestAbs(t, integ_f_sin(0.0, gsl.Pi), 2.0, gsl.Float64Eps, "integ_f_sin sanity check")
	}

	/* Test the fixed-order Gauss-Legendre rules against sin(x) on [0, pi] */
	{
		n_max := 1024
		f := &f_sin{}
		a := 0.0
		b := gsl.Pi
		expected := integ_f_sin(a, b)
		var result, abserr, prev_abserr float64

		for n := 1; n <= n_max; n++ {
			tbl, _ := NewGLFixedTable(n)

			result = GLFixed(f, a, b, tbl)
			abserr = math.Abs(expected - result)

			if n == 1 {
				test.TestAbs(t, result, f.Evaluate((b+a)/2)*(b-a), 0.0, fmt.Sprintf("glfixed %d-point: behavior for n == 1", n))
			} else if n < 9 {
				s := 0
				if !(abserr < prev_abserr) {
					s = 1
				}

				test.Test(t, s, fmt.Sprintf("glfixed %d-point: observed drop in absolute error versus %d-points", n, n-1))
			} else if tbl.precomputed == 1 {
				test.TestAbs(t, result, expected, 2.0*float64(n)*gsl.Float64Eps, fmt.Sprintf("glfixed %d-point: very low absolute error for high precision coefficients", n))
			} else {
				test.TestAbs(t, result, expected, 1.0e6*gsl.Float64Eps, fmt.Sprintf("glfixed %d-point: acceptable absolute error for on-the-fly coefficients", n))
			}

			prev_abserr = abserr
		}
	}

	/* Test some fixed-order Gauss-Legendre rule points and weights on [-1, 1] */
	/* This verifies the (point, weight) retrieval API behaves sanely */
	{
		eps := gsl.Float64Eps
		var n, i int
		var xi, wi float64

		/* Analytical results for points and weights on [-1, 1]
		   Pulled from http://en.wikipedia.org/wiki/Gaussian_quadrature
		   Sorted in increasing order of Gauss points */
		e1 := [][]float64{
			{0, 2},
		}
		e2 := [][]float64{
			{-1.0 / gsl.Sqrt3, 1},
			{1.0 / gsl.Sqrt3, 1},
		}
		e3 := [][]float64{
			{-sqrt15 / 5, 5. / 9},
			{0, 8. / 9},
			{sqrt15 / 5, 5. / 9},
		}
		e4 := [][]float64{
			{-const1, (18 - sqrt30) / 36},
			{-const2, (18 + sqrt30) / 36},
			{const2, (18 + sqrt30) / 36},
			{const1, (18 - sqrt30) / 36},
		}
		e5 := [][]float64{
			{-const3, (322 - 13*sqrt70) / 900},
			{-const4, (322 + 13*sqrt70) / 900},
			{0, 128. / 225},
			{const4, (322 + 13*sqrt70) / 900},
			{const3, (322 - 13*sqrt70) / 900},
		}

		n = 1
		tbl, _ := NewGLFixedTable(n)
		for i = 0; i < n; i++ {
			GLFixedPoint(-1, 1, i, &xi, &wi, tbl)
			test.TestAbs(t, xi, e1[i][0], eps, fmt.Sprintf("glfixed %d-point lookup: x(%d)", n, i))
			test.TestAbs(t, wi, e1[i][1], eps, fmt.Sprintf("glfixed %d-point lookup: w(%d)", n, i))
		}
		n = 2
		tbl2, _ := NewGLFixedTable(n)
		for i = 0; i < n; i++ {
			GLFixedPoint(-1, 1, i, &xi, &wi, tbl2)
			test.TestAbs(t, xi, e2[i][0], eps, fmt.Sprintf("glfixed %d-point lookup: x(%d)", n, i))
			test.TestAbs(t, wi, e2[i][1], eps, fmt.Sprintf("glfixed %d-point lookup: w(%d)", n, i))
		}

		n = 3
		tbl3, _ := NewGLFixedTable(n)
		for i = 0; i < n; i++ {
			GLFixedPoint(-1, 1, i, &xi, &wi, tbl3)
			test.TestAbs(t, xi, e3[i][0], eps, fmt.Sprintf("glfixed %d-point lookup: x(%d)", n, i))
			test.TestAbs(t, wi, e3[i][1], eps, fmt.Sprintf("glfixed %d-point lookup: w(%d)", n, i))
		}

		n = 4
		tbl4, _ := NewGLFixedTable(n)
		for i = 0; i < n; i++ {
			GLFixedPoint(-1, 1, i, &xi, &wi, tbl4)
			test.TestAbs(t, xi, e4[i][0], eps, fmt.Sprintf("glfixed %d-point lookup: x(%d)", n, i))
			test.TestAbs(t, wi, e4[i][1], eps, fmt.Sprintf("glfixed %d-point lookup: w(%d)", n, i))
		}

		n = 5
		tbl5, _ := NewGLFixedTable(n)
		for i = 0; i < n; i++ {
			GLFixedPoint(-1, 1, i, &xi, &wi, tbl5)
			test.TestAbs(t, xi, e5[i][0], eps, fmt.Sprintf("glfixed %d-point lookup: x(%d)", n, i))
			test.TestAbs(t, wi, e5[i][1], eps, fmt.Sprintf("glfixed %d-point lookup: w(%d)", n, i))
		}
	}

}
