package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselY0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, -1.5342386513503668441, TEST_TOL0},
		{2, 0.5103756726497451196, TEST_TOL0},
		{256.0, -0.03381290171792454909, TEST_TOL0},
		{4294967296.0, 3.657903190017678681e-06, TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Y0_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Y0_e(%v)", c.x))
		})
	}
}
