package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestSi(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-1.0, -0.9460830703671830149, TEST_TOL0},
		{1.0e-10, 1.0e-10, TEST_TOL0},
		{1.0e-05, 9.999999999944444444e-06, TEST_TOL0},
		{0.1, 0.09994446110827695016, TEST_TOL0},
		{1.0, 0.9460830703671830149, TEST_TOL0},
		{10.0, 1.6583475942188740493, TEST_TOL0},
		{50.0, 1.5516170724859358947, TEST_TOL0},
		{300.0, 1.5708810882137495193, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Si_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Si_e(%v)", c.x))
		})
	}
}

func TestCi(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0 / 4294967296.0, -21.603494113016717041, TEST_TOL0},
		{1.0 / 65536.0, -10.513139224115799751, TEST_TOL0},
		{1.0 / 8.0, -1.5061295845296396649, TEST_TOL0},
		{1.0, 0.3374039229009681347, TEST_TOL0},
		{10.0, -0.04545643300445537263, TEST_TOL0},
		{50.0, -0.005628386324116305440, TEST_TOL0},
		{300.0, -0.003332199918592111780, TEST_TOL0},
		{65536.0, 0.000010560248837656279453, TEST_TOL0},
		{4294967296.0, -1.0756463261957757485e-10, TEST_SQRT_TOL0},
		{1099511627776.0, -3.689865584710764214e-13, 1024.0 * TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ci_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Ci_e(%v)", c.x))
		})
	}
}
