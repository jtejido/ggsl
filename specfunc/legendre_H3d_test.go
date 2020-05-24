package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestLegendreH3d0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		lambda, eta, expected, tol float64
	}{
		{1.0e-06, 1.0e-06, 0.9999999999998333333, TEST_TOL0},
		{1.0, 0.0, 1.0, TEST_TOL0},
		{1.0, 1.0, 0.7160229153604338713, TEST_TOL0},
		{1.0, 100.0, -3.767437313149604566e-44, TEST_TOL2},
		{1.0, 500.0, -6.665351935878582205e-218, TEST_TOL2},
		{100.0, 1.0, -0.004308757035378200029, TEST_TOL0},
		{100.0, 10.0, 7.508054627912986427e-07, TEST_TOL0},
		{1000.0, 1.0, 0.0007036067909088818319, TEST_TOL0},
		{1.0e+08, 1.0, 7.927485371429105968e-09, TEST_TOL3},
		{1.0e+08, 100.0, -3.627118904186918957e-52, 32.0 * TEST_SQRT_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_H3d_0_e(c.lambda, c.eta, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_H3d_0_e(%v,%v)", c.lambda, c.eta))

		})
	}
}

func TestLegendreH3d1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		lambda, eta, expected, tol float64
	}{
		{1.0e-06, 1.0e-06, 3.333333333334222222e-07, TEST_TOL0},
		{1.0, 1.0e-10, 4.714045207910316829e-11, TEST_TOL0},
		{1.0, 1.0, 0.3397013994799344639, TEST_TOL0},
		{1.0, 100.0, -7.200624449531811272e-44, TEST_TOL2},
		{1.0, 500.0, 4.192260336821728677e-218, TEST_TOL2},
		{100.0, 0.01, 0.30117664944267412324, TEST_TOL1},
		{100.0, 1.0, -0.007393833425336299309, TEST_TOL0},
		{100.0, 10.0, -5.031062029821254982e-07, TEST_TOL0},
		{1000.0, 0.001, 0.30116875865090396421, TEST_TOL0},
		{1000.0, 1.0, -0.0004776144516074971885, TEST_TOL0},
		{1.0e+08, 1.0e-08, 0.30116867893975679722, TEST_TOL1},
		{1.0e+08, 1.0, 3.0921097047369081582e-09, TEST_TOL4},
		{1.0e+08, 100.0, -6.496142701296286936e-52, 32.0 * TEST_SQRT_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_H3d_1_e(c.lambda, c.eta, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_H3d_1_e(%v,%v)", c.lambda, c.eta))

		})
	}
}

func TestLegendreH3d(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l                          int
		lambda, eta, expected, tol float64
	}{
		{5, 1.0e-06, 1.0e-06, 1.1544011544013627977e-32, TEST_TOL2},
		{5, 1.0, 1.0e-10, 2.0224912016958766992e-52, TEST_TOL2},
		{5, 1.0, 1.0, 0.011498635037491577728, TEST_TOL1},
		{5, 1.0, 5.0, 0.0020696945662545205776, TEST_TOL4},
		{5, 1.0, 7.0, -0.0017555303787488993676, TEST_TOL4},
		{5, 1.0, 10.0, 0.00008999979724504887101, TEST_TOL2},
		{5, 1.0, 100.0, -4.185397793298567945e-44, TEST_TOL2},
		{5, 1.0, 500.0, 1.4235113901091961263e-217, TEST_TOL3},
		{5, 100.0, 0.001, 9.642762597222417946e-10, TEST_TOL2},
		{5, 100.0, 0.002, 3.0821201254308036109e-08, TEST_TOL2},
		{5, 100.0, 0.01, 0.00009281069019005840532, TEST_TOL1},
		{5, 100.0, 1.0, -0.008043100696178624653, TEST_TOL2},
		{5, 100.0, 10.0, -3.927678432813974207e-07, TEST_TOL3},
		{5, 1000.0, 0.001, 0.00009256365284253254503, TEST_TOL1},
		{5, 1000.0, 0.01, -0.05553733815473079983, TEST_TOL0},
		{5, 1.0e+08, 1.0e-08, 0.00009256115861125841299, TEST_TOL2},
		{5, 1.0e+08, 100.0, -6.496143209092860765e-52, 128.0 * TEST_SQRT_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_H3d_e(c.l, c.lambda, c.eta, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_H3d_e(%v,%v,%v)", c.l, c.lambda, c.eta))

		})
	}
}
