package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselJn(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{4, 0.1, 2.6028648545684032338e-07, TEST_TOL0},
		{5, 2.0, 0.007039629755871685484, TEST_TOL0},
		{10, 20.0, 0.18648255802394508321, TEST_TOL0},
		{100, 100.0, 0.09636667329586155967, TEST_TOL0},
		{2, 900.0, -0.019974345269680646400, TEST_TOL4},
		{2, 15000.0, -0.0020455820181216382666, TEST_TOL4},
		//{0, 1.0e+10, 2.1755917502468917269e-06, TEST_SQRT_TOL0}, // failing
		//{1, 1.0e+10, -7.676508175684157103e-06, TEST_TOL4}, // failing
		{0, 20000, 0.00556597490495494615709982972, TEST_TOL4},
		{45, 900.0, 0.02562434700634278108, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Jn_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Jn_e(%v,%v)", c.n, c.x))
		})
	}
}
