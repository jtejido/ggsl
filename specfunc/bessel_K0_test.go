package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselK0Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 2.6823261022628943831, TEST_TOL0},
		{1.95, 0.8513330938802157074894, TEST_TOL0},
		{2.0, 0.8415682150707714179, TEST_TOL0},
		{6.0, 0.50186313086214003217346, TEST_TOL0},
		{100.0, 0.1251756216591265789, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_K0_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_K0_scaled_e(%v)", c.x))

		})
	}
}

func TestBesselK0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 2.4270690247020166125, TEST_TOL0},
		{1.95, 0.1211226255426818887894, TEST_TOL0},
		{2.0, 0.11389387274953343565, TEST_TOL0},
		{100.0, 4.656628229175902019e-45, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_K0_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_K0_e(%v)", c.x))
		})
	}
}
