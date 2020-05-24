package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestGegenpoly1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		lambda, x, expected float64
		tol                 float64
	}{
		{-0.2, 1.0, -0.4, TEST_TOL0},
		{0.0, 1.0, 2.0, TEST_TOL0},
		{1.0, 1.0, 2.0, TEST_TOL0},
		{1.0, 0.5, 1.0, TEST_TOL0},
		{5.0, 1.0, 10.0, TEST_TOL0},
		{100.0, 0.5, 100.0, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gegenpoly_1_e(c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gegenpoly_1_e(%v,%v)", c.lambda, c.x))
		})
	}
}

func TestGegenpoly2(t *testing.T) {
	r := new(Result)
	cases := []struct {
		lambda, x, expected float64
		tol                 float64
	}{
		{-0.2, 0.5, 0.12, TEST_TOL0},
		{0.0, 1.0, 1.00, TEST_TOL0},
		{1.0, 1.0, 3.00, TEST_TOL0},
		{1.0, 0.1, -0.96, TEST_TOL0},
		{5.0, 1.0, 55.0, TEST_TOL0},
		{100.0, 0.5, 4950.0, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gegenpoly_2_e(c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gegenpoly_2_e(%v,%v)", c.lambda, c.x))
		})
	}
}

func TestGegenpoly3(t *testing.T) {
	r := new(Result)
	cases := []struct {
		lambda, x, expected float64
		tol                 float64
	}{
		{-0.2, 0.5, 0.112, TEST_TOL0},
		{0.0, 1.0, -2.0 / 3.0, TEST_TOL0},
		{1.0, 1.0, 4.000, TEST_TOL0},
		{1.0, 0.1, -0.392, TEST_TOL0},
		{5.0, 1.0, 220.000, TEST_TOL0},
		{100.0, 0.5, 161600.000, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gegenpoly_3_e(c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gegenpoly_3_e(%v,%v)", c.lambda, c.x))
		})
	}
}

func TestGegenpolyN(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n                   int
		lambda, x, expected float64
		tol                 float64
	}{
		{1, 1.0, 1.0, 2.000, TEST_TOL0},
		{10, 1.0, 1.0, 11.000, TEST_TOL0},
		{10, 1.0, 0.1, -0.4542309376, TEST_TOL0},
		{10, 5.0, 1.0, 9.23780e+4, TEST_TOL0},
		{10, 100.0, 0.5, 1.5729338392690000e+13, TEST_TOL0},
		{1000, 100.0, 1.0, 3.3353666135627322e+232, TEST_TOL1},
		{100, 2000.0, 1.0, 5.8753432034937579e+202, TEST_TOL0},
		{103, 207.0, 2.0, 1.4210272202235983e+145, TEST_TOL0},
		{103, -0.4, 0.3, -1.64527498094522e-04, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gegenpoly_n_e(c.n, c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gegenpoly_n_e(%v,%v,%v)", c.n, c.lambda, c.x))
		})
	}
}
