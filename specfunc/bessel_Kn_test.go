package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselKnScaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{4, 0.1, 530040.2483725626207, TEST_TOL1},
		{5, 2.0, 69.68655087607675118, TEST_TOL0},
		{100, 100.0, 2.0475736731166756813e+19, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Kn_scaled_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Kn_scaled_e(%v,%v)", c.n, c.x))
		})
	}
}

func TestBesselKn(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{4, 0.1, 479600.2497925682849, TEST_TOL1},
		{5, 2.0, 9.431049100596467443, TEST_TOL0},
		{100, 100.0, 7.617129630494085416e-25, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Kn_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Kn_e(%v,%v)", c.n, c.x))
		})
	}
}
