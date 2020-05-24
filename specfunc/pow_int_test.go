package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestPowInt(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x             float64
		n             int
		expected, tol float64
	}{
		{2.0, 3, 8.0, TEST_TOL0},
		{-2.0, 3, -8.0, TEST_TOL0},
		{2.0, -3, 1.0 / 8.0, TEST_TOL0},
		{-2.0, -3, -1.0 / 8.0, TEST_TOL0},

		{10.0, 4, 1.0e+4, TEST_TOL0},
		{10.0, -4, 1.0e-4, TEST_TOL0},
		{-10.0, 4, 1.0e+4, TEST_TOL0},
		{-10.0, -4, 1.0e-4, TEST_TOL0},

		{10.0, 40, 1.0e+40, TEST_TOL0},
		{8.0, -40, 7.523163845262640051e-37, TEST_TOL0},
		{-10.0, 40, 1.0e+40, TEST_TOL0},
		{-8.0, -40, 7.523163845262640051e-37, TEST_TOL0},

		{10.0, 41, 1.0e+41, TEST_TOL0},
		{8.0, -41, 9.403954806578300064e-38, TEST_TOL0},
		{-10.0, 41, -1.0e+41, TEST_TOL0},
		{-8.0, -41, -9.403954806578300064e-38, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Pow_int_e(c.x, c.n, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Pow_int_e(%v, %v)", c.x, c.n))

		})
	}
}
