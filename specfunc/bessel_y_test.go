package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestBessely0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.001, -999.99950000004166670, TEST_TOL0},
		{1.0, -0.5403023058681397174, TEST_TOL0},
		{10.0, 0.08390715290764524523, TEST_TOL0},
		{100.0, -0.008623188722876839341, TEST_TOL0},
		{65536.0, 0.000011014324202158573930, TEST_TOL0},
		{4294967296.0, 2.0649445131370357007e-10, TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_y0_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_y0_e(%v)", c.x))
		})
	}
}

func TestBessely1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.01, -10000.499987500069444, TEST_TOL0},
		{1.0, -1.3817732906760362241, TEST_TOL0},
		{10.0, 0.06279282637970150586, TEST_TOL0},
		{100.0, 0.004977424523868819543, TEST_TOL0},
		{4294967296.0, 1.0756463271573404688e-10, TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_y1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_y1_e(%v)", c.x))
		})
	}
}

func TestBessely2(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.01, -3.0000500012499791668e+06, TEST_TOL0},
		{1.0, -3.605017566159968955, TEST_TOL0},
		{10.0, -0.06506930499373479347, TEST_TOL0},
		{100.0, 0.008772511458592903927, TEST_TOL0},
		{4294967296.0, -2.0649445123857054207e-10, TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_y2_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_y2_e(%v)", c.x))
		})
	}
}

func TestBesselyl(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l           int
		x, expected float64
		tol         float64
	}{
		{0, 0.01, -99.995000041666528, TEST_TOL0},
		{0, 1.0, -0.54030230586813972, TEST_TOL0},
		{1, 10.0, 0.062792826379701506, TEST_TOL0},
		{5, 1.0, -999.44034339223641, TEST_TOL0},
		{10, 0.01, -6.5473079797378378e+30, TEST_TOL0},
		{10, 10.0, -0.172453672088057849, TEST_TOL0},
		{100, 1.0, -6.6830794632586775e+186, TEST_TOL1},
		{100, 100.0, -0.0229838504915622811, TEST_TOL1},
		{2000, 1048576.0, 5.9545201447146155e-07, TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_yl_e(c.l, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_yl_e(%v, %v)", c.l, c.x))
		})
	}
}
