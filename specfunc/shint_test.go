package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestShi(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-1.0, -1.0572508753757285146, TEST_TOL0},
		{1.0 / 4294967296.0, 2.3283064365386962891e-10, TEST_TOL0},
		{1.0 / 65536.0, 0.00001525878906269737298, TEST_TOL0},
		{0.1, 0.1000555722250569955, TEST_TOL0},
		{1.0, 1.0572508753757285146, TEST_TOL0},
		{10.0, 1246.1144901994233444, TEST_TOL1},
		{50.0, 5.292818448565845482e+19, TEST_TOL2},
		{300.0, 3.248241254044332895e+127, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Shi_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Shi_e(%v)", c.x))
		})
	}
}

func TestChi(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-1.0, 0.8378669409802082409, TEST_TOL0},
		{1.0 / 4294967296.0, -21.603494113016717041, TEST_TOL0},
		{1.0 / 65536.0, -10.513139223999384429, TEST_TOL0},
		{1.0 / 8.0, -1.4983170827635760646, TEST_TOL0},
		{1.0, 0.8378669409802082409, TEST_TOL0},
		{10.0, 1246.1144860424544147, TEST_TOL1},
		{50.0, 5.292818448565845482e+19, TEST_TOL2},
		{300.0, 3.248241254044332895e+127, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Chi_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Chi_e(%v)", c.x))
		})
	}
}
