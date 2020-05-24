package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestAtanint(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0e-10, 1.0e-10, TEST_TOL0},
		{1.0e-05, 9.99999999988888888889e-06, TEST_TOL0},
		{0.1, 0.09988928686033618404, TEST_TOL0},
		{1.0, 0.91596559417721901505, TEST_TOL0},
		{2.0, 1.57601540344632342236, TEST_TOL0},
		{10.0, 3.71678149306806859029, TEST_TOL0},
		{50.0, 6.16499047850274874222, TEST_TOL0},
		{300.0, 8.96281388924518959990, TEST_TOL0},
		{1.0e+5, 18.084471031038661920, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Atanint_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Atanint_e(%v)", c.x))
		})
	}
}
