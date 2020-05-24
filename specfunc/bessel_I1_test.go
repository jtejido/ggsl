package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselI1Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 0.04529844680880932501, TEST_TOL0},
		{2, 0.21526928924893765916, TEST_TOL0},
		{100.0, 0.03974415302513025267, TEST_TOL0},
		{65536.0, 0.0015583593657207350452, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_I1_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_I1_scaled_e(%v)", c.x))
		})
	}
}

func TestBesselI1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 0.05006252604709269211, TEST_TOL0},
		{2.0, 1.59063685463732906340, TEST_TOL0},
		{100.0, 1.0683693903381624812e+42, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_I1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_I1_e(%v)", c.x))
		})
	}
}
