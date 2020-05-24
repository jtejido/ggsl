package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselInScaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{-4, 0.1, 2.3575258620054605307e-07, TEST_TOL0},
		{4, 0.1, 2.3575258620054605307e-07, TEST_TOL0},
		{5, 2.0, 0.0013297610941881578142, TEST_TOL0},
		{100, 100.0, 1.7266862628167695785e-22, TEST_TOL0},
		{2, 1e7, 1.261566024466416433e-4, TEST_TOL2},
		{2, 1e8, 3.989422729212649531e-5, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_In_scaled_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_In_scaled_e(%v,%v)", c.n, c.x))
		})
	}
}

func TestBesselIn(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{4, 0.1, 2.6054690212996573677e-07, TEST_TOL0},
		{5, 2.0, 0.009825679323131702321, TEST_TOL0},
		{100, 100.0, 4.641534941616199114e+21, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_In_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_In_e(%v,%v)", c.n, c.x))
		})
	}
}
