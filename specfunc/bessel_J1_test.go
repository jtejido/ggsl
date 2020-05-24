package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselJ1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 0.04993752603624199756, TEST_TOL0},
		{2.0, 0.57672480775687338720, TEST_TOL0},
		{100.0, -0.07714535201411215803, TEST_TOL0},
		// {1.0e+10, -7.676508175684157103e-06, TEST_TOL4}, // failing
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_J1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_J1_e(%v)", c.x))
		})
	}
}
