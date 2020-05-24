package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesseli0Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.0, 1.0, TEST_TOL0},
		{0.1, 0.9063462346100907067, TEST_TOL0},
		{2.0, 0.24542109027781645493, TEST_TOL0},
		{100.0, 0.005000000000000000000, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_i0_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_i0_scaled_e(%v)", c.x))
		})
	}
}

func TestBesseli1Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.0, 0.0, TEST_TOL0},
		{0.1, 0.030191419289002226846, TEST_TOL0},
		{2.0, 0.131868364583275317610, TEST_TOL0},
		{100.0, 0.004950000000000000000, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_i1_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_i1_scaled_e(%v)", c.x))
		})
	}
}

func TestBesseli2Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.0, 0.0, TEST_TOL0},
		{0.1, 0.0006036559400239012567, TEST_TOL0},
		{2.0, 0.0476185434029034785100, TEST_TOL0},
		{100.0, 0.0048515000000000000000, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_i2_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_i2_scaled_e(%v)", c.x))
		})
	}
}

func TestBesselilScaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l           int
		x, expected float64
		tol         float64
	}{
		{0, 0.0, 1.0, TEST_TOL0},
		{1, 0.0, 0.0, TEST_TOL0},
		{4, 0.001, 1.0571434341190365013e-15, TEST_TOL0},
		{4, 0.1, 9.579352242057134927e-08, TEST_TOL1},
		{5, 2.0, 0.0004851564602127540059, TEST_TOL0},
		{5, 100.0, 0.004300446777500000000, TEST_TOL1},
		{100, 100.0, 1.3898161964299132789e-23, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_il_scaled_e(c.l, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_il_scaled_e(%v, %v)", c.l, c.x))
		})
	}
}
