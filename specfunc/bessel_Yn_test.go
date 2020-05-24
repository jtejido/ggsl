package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselYn(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{4, 0.1, -305832.29793353160319, TEST_TOL1},
		{5, 2, -9.935989128481974981, TEST_TOL0},
		{100, 100.0, -0.16692141141757650654, TEST_TOL0},
		{100, 4294967296.0, 3.657889671577715808e-06, TEST_SQRT_TOL0},
		{1000, 4294967296.0, 3.656551321485397501e-06, 2.0e-05},
		{2, 15000.0, -0.006185217273358617849, TEST_TOL4},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Yn_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Yn_e(%v,%v)", c.n, c.x))
		})
	}
}
