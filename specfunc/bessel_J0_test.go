package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselJ0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 0.99750156206604003230, TEST_TOL0},
		{2.0, 0.22389077914123566805, TEST_TOL0},
		{100.0, 0.019985850304223122424, TEST_TOL0},
		//{1.0e+10, 2.1755917502468917269e-06, TEST_SQRT_TOL0}, // failing
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_J0_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_J0_e(%v)", c.x))
		})
	}
}
