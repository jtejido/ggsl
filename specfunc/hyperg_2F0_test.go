package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestHyperg2f0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, b, x, expected float64
		tol               float64
	}{
		{0.01, 1.0, -0.02, .99980388665511730901180717, TEST_TOL0},
		{0.1, 0.5, -0.02, .99901595171179281891589794, TEST_TOL0},
		{1, 1, -0.02, .98075549650574351826538049000, TEST_TOL0},
		{8, 8, -0.02, .32990592849626965538692141, TEST_TOL0},
		{50, 50, -0.02, .2688995263772964415245902e-12, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Hyperg_2F0_e(c.a, c.b, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Hyperg_2F0_e(%v,%v,%v)", c.a, c.b, c.x))
		})
	}
}
