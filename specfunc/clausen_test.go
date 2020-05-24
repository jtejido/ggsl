package specfunc

import (
	"fmt"
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestClausen(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{gsl.Pi / 20.0, 0.4478882448133546, TEST_TOL0},
		{gsl.Pi / 6.0, 0.8643791310538927, TEST_TOL0},
		{gsl.Pi / 3.0, 1.0149416064096535, TEST_TOL0},
		{2.0*gsl.Pi + gsl.Pi/3.0, 1.0149416064096535, TEST_TOL0},
		{100.0*gsl.Pi + gsl.Pi/3.0, 1.0149416064096535, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Clausen_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Clausen_e(%v)", c.x))
		})
	}
}
