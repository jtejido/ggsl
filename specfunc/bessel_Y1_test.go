package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselY1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, -6.45895109470202698800, TEST_TOL0},
		{2, -0.10703243154093754689, TEST_TOL0},
		{100.0, -0.020372312002759793305, TEST_TOL0},
		{4294967296.0, 0.000011612249378370766284, TEST_TOL4},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Y1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Y1_e(%v)", c.x))
		})
	}
}
