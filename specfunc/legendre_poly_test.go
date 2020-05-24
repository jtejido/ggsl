package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestLegendreP1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{-0.5, -0.5, TEST_TOL0},
		{0.5, 0.5, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_P1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_P1_e(%v)", c.x))

		})
	}
}

func TestLegendreP2(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{0.0, -0.5, TEST_TOL0},
		{0.5, -0.125, TEST_TOL0},
		{1.0, 1.0, TEST_TOL0},
		{100.0, 14999.5, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_P2_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_P2_e(%v)", c.x))

		})
	}
}

func TestLegendreP3(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{-0.5, 0.4375, TEST_TOL0},
		{0.5, -0.4375, TEST_TOL0},
		{1.0, 1.0, TEST_TOL0},
		{100.0, 2.49985e+06, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_P3_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_P3_e(%v)", c.x))

		})
	}
}

func TestLegendrePl(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l                int
		x, expected, tol float64
	}{
		{1, -0.5, -0.5, TEST_TOL0},
		{1, 1.0e-8, 1.0e-08, TEST_TOL0},
		{1, 0.5, 0.5, TEST_TOL0},
		{1, 1.0, 1.0, TEST_TOL0},

		{10, -0.5, -0.1882286071777345, TEST_TOL0},
		{10, 1.0e-8, -0.24609374999999864648, TEST_TOL0},
		{10, 0.5, -0.18822860717773437500, TEST_TOL0},
		{10, 1.0, 1.0, TEST_TOL0},

		{99, -0.5, 0.08300778172138770477, TEST_TOL0},
		{99, 1.0e-8, -7.958923738716563193e-08, TEST_TOL0},
		{99, 0.5, -0.08300778172138770477, TEST_TOL0},
		{99, 0.999, -0.3317727359254778874, TEST_TOL2},
		{99, 1.0, 1.0, TEST_TOL0},

		{1000, -0.5, -0.019168251091650277878, TEST_TOL2},
		{1000, 1.0e-8, 0.0252250181770982897470252620, TEST_TOL2},
		{1000, 0.5, -0.019168251091650277878, TEST_TOL2},
		{1000, 1.0, 1.0, TEST_TOL0},

		{4000, -0.5, -0.009585404456573080972, TEST_TOL2},
		{4000, 0.5, -0.009585404456573080972, TEST_TOL2},
		{4000, 1.0, 1.0, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_Pl_e(c.l, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_Pl_e(%v, %v)", c.l, c.x))

		})
	}
}

func TestLegendrePlm(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l, m             int
		x, expected, tol float64
	}{
		{10, 0, -0.5, -0.18822860717773437500, TEST_TOL0},
		{10, 0, 1.0e-08, -0.24609374999999864648, TEST_TOL0},
		{10, 0, 0.5, -0.18822860717773437500, TEST_TOL0},

		{10, 1, -0.5, -2.0066877394361256516, TEST_TOL0},
		{10, 1, 1.0e-08, -2.7070312499999951725e-07, TEST_TOL0},
		{10, 1, 0.5, 2.0066877394361256516, TEST_TOL0},

		{10, 5, -0.5, -30086.169706116174977, TEST_TOL0},
		{10, 5, 1.0e-08, -0.0025337812499999964949, TEST_TOL0},
		{10, 5, 0.5, 30086.169706116174977, TEST_TOL0},
		{10, 5, 0.999, -0.5036411489013270406, TEST_TOL1},

		{100, 5, -0.5, -6.617107444248382171e+08, TEST_TOL0},
		{100, 5, 1.0e-08, 817.8987598063712851, TEST_TOL0},
		{100, 5, 0.5, 6.617107444248382171e+08, TEST_TOL0},
		{100, 5, 0.999, -1.9831610803806212189e+09, TEST_TOL2},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_Plm_e(c.l, c.m, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_Plm_e(%v, %v, %v)", c.l, c.m, c.x))

		})
	}
}

func TestLegendreSphPlm(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l, m             int
		x, expected, tol float64
	}{
		{10, 0, -0.5, -0.24332702369300133776, TEST_TOL0},
		{10, 0, 0.5, -0.24332702369300133776, TEST_TOL0},
		{10, 0, 0.999, 1.2225754122797385990, TEST_TOL1},

		{10, 5, -0.5, -0.3725739049803293972, TEST_TOL0},
		{10, 5, 1.0e-08, -3.1377233589376792243e-08, TEST_TOL0},
		{10, 5, 0.5, 0.3725739049803293972, TEST_TOL0},
		{10, 5, 0.999, -6.236870674727370094e-06, TEST_TOL2},

		{10, 10, -0.5, 0.12876871185785724117, TEST_TOL1},
		{10, 10, 0.5, 0.12876871185785724117, TEST_TOL1},
		{10, 10, 0.999, 1.7320802307583118647e-14, TEST_TOL2},

		{200, 1, -0.5, 0.3302975570099492931, TEST_TOL1},
		{200, 1, 0.5, -0.3302975570099492931, TEST_TOL1},
		{200, 1, 0.999, -1.4069792055546256912, TEST_TOL2},

		/* Test case from alberto@physik.fu-berlin.de */

		{3, 1, 0.0, 0.323180184114150653007, TEST_TOL2},

		/* Other test cases */

		{200, 1, -0.5, 0.3302975570099492931418227583, TEST_TOL2},
		{140, 135, 1, 0.0, TEST_TOL2},
		// #ifdef EXTENDED
		{140, 135, 0.99998689456491752, -6.54265253269093276310395668335e-305, TEST_TOL6},
		// #endif
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_sphPlm_e(c.l, c.m, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_sphPlm_e(%v, %v, %v)", c.l, c.m, c.x))

		})
	}
}
