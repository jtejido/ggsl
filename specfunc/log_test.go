package specfunc

import (
	"fmt"
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestLog(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{0.1, -2.3025850929940456840, TEST_TOL0},
		{1.1, 0.09531017980432486004, TEST_TOL1},
		{1000.0, 6.907755278982137052, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Log_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Log_e(%v)", c.x))

		})
	}
}

func TestComplexLog(t *testing.T) {
	r1, r2 := new(Result), new(Result)
	cases := []struct {
		zr, zi, expected1, tol1, expected2, tol2 float64
	}{
		{1.0, 1.0,
			0.3465735902799726547, TEST_TOL0,
			0.7853981633974483096, TEST_TOL0},

		{1.0, -1.0,
			0.3465735902799726547, TEST_TOL0,
			-0.7853981633974483096, TEST_TOL0},

		{1.0, 100.0,
			4.605220183488258022, TEST_TOL0,
			1.560796660108231381, TEST_TOL0},

		{-1000.0, -1.0,
			6.907755778981887052, TEST_TOL0,
			-3.1405926539231263718, TEST_TOL0},

		{-1.0, 0.0,
			0.0, TEST_TOL0,
			3.1415926535897932385, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Complex_log_e(c.zr, c.zi, r1, r2)
			run_test_sf_2(t, stat, r1, c.expected1, c.tol1, r2, c.expected2, c.tol2, err.SUCCESS, fmt.Sprintf("Complex_log_e(%v,%v)", c.zr, c.zi))
		})
	}
}

func TestLogAbs(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{-0.1, -2.3025850929940456840, TEST_TOL0},
		{-1.1, 0.09531017980432486004, TEST_TOL1},
		{-1000.0, 6.907755278982137052, TEST_TOL0},
		{0.1, -2.3025850929940456840, TEST_TOL0},
		{1.1, 0.09531017980432486004, TEST_TOL1},
		{1000.0, 6.907755278982137052, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Log_abs_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Log_abs_e(%v)", c.x))

		})
	}
}

func TestLog1plusx(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{1.0e-10, 9.999999999500000000e-11, TEST_TOL0},
		{1.0e-8, 9.999999950000000333e-09, TEST_TOL0},
		{1.0e-4, 0.00009999500033330833533, TEST_TOL0},
		{0.1, 0.09531017980432486004, TEST_TOL0},
		{0.49, 0.3987761199573677730, TEST_TOL0},

		{-0.49, -0.6733445532637655964, TEST_TOL0},
		{1.0, gsl.Ln2, TEST_TOL0},
		{-0.99, -4.605170185988091368, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Log_1plusx_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Log_1plusx_e(%v)", c.x))

		})
	}
}

func TestLog1plusxMx(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{1.0e-10, -4.999999999666666667e-21, TEST_TOL0},
		{1.0e-8, -4.999999966666666917e-17, TEST_TOL0},
		{1.0e-4, -4.999666691664666833e-09, TEST_TOL0},
		{0.1, -0.004689820195675139956, TEST_TOL0},
		{0.49, -0.09122388004263222704, TEST_TOL0},

		{-0.49, -0.18334455326376559639, TEST_TOL0},
		{1.0, gsl.Ln2 - 1.0, TEST_TOL0},
		{-0.99, -3.615170185988091368, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Log_1plusx_mx_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Log_1plusx_mx_e(%v)", c.x))

		})
	}
}
