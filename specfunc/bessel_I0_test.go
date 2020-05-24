package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselI0Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1e-10, 0.99999999990000000001, TEST_TOL0},
		{0.1, 0.90710092578230109640, TEST_TOL0},
		{2, 0.30850832255367103953, TEST_TOL0},
		{100.0, 0.03994437929909668265, TEST_TOL0},
		{65536.0, 0.0015583712551952223537, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_I0_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_I0_scaled_e(%v)", c.x))
		})
	}
}

func TestBesselI0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 1.0025015629340956014, TEST_TOL0},
		{2.0, 2.2795853023360672674, TEST_TOL0},
		{100.0, 1.0737517071310738235e+42, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_I0_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_I0_e(%v)", c.x))
		})
	}
}
