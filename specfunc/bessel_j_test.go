package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselj0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10.0, -0.05440211108893698134, TEST_TOL0},
		{0.001, 0.9999998333333416667, TEST_TOL0},
		{1.0, 0.84147098480789650670, TEST_TOL0},
		{10.0, -0.05440211108893698134, TEST_TOL0},
		{100.0, -0.005063656411097587937, TEST_TOL1},
		// #ifdef FIXME
		//   {1048576.0,  3.1518281938718287624e-07, TEST_TOL2},
		// #endif

		//   /* these values are from Mathematica */
		// #ifdef FIXME
		//   {1.0e18,  -9.9296932074040507620955e-19, TEST_TOL0},
		//   {1.0e20,  -6.4525128526578084420581e-21, TEST_TOL0},
		// #endif
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_j0_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_j0_e(%v)", c.x))
		})
	}
}

func TestBesselj1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10.0, -0.07846694179875154709, TEST_TOL0},
		{0.01, 0.003333300000119047399, TEST_TOL0},
		{1.0, 0.30116867893975678925, TEST_TOL0},
		{10.0, 0.07846694179875154709, TEST_TOL0},
		{100.0, -0.008673825286987815220, TEST_TOL0},
		{1048576.0, -9.000855242905546158e-07, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_j1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_j1_e(%v)", c.x))
		})
	}
}

func TestBesselj2(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10.0, 0.07794219362856244547, TEST_TOL0},
		{0.01, 6.666619047751322551e-06, TEST_TOL0},
		{1.0, 0.06203505201137386110, TEST_TOL0},
		{10.0, 0.07794219362856244547, TEST_TOL0},
		{100.0, 0.004803441652487953480, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_j2_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_j2_e(%v)", c.x))
		})
	}
}

func TestBesseljl(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l           int
		x, expected float64
		tol         float64
	}{
		{0, 0.0, 1.0, TEST_TOL0},

		{1, 10.0, 0.07846694179875154709000, TEST_TOL0},
		{5, 1.0, 0.00009256115861125816357, TEST_TOL0},
		{10, 10.0, 0.06460515449256426427, TEST_TOL0},
		{100, 100.0, 0.010880477011438336539, TEST_TOL1},
		{2000, 1048576.0, 7.449384239168568534e-07, TEST_SQRT_TOL0},

		/* related to BUG#3 problem */
		{2, 900.0, -0.0011089115568832940086, TEST_TOL4},
		{2, 15000.0, -0.00005955592033075750554, TEST_TOL4},

		/* Bug report by Mario Santos, value computed from AS 10.1.8 */
		{100, 1000.0, -0.00025326311230945818285, TEST_TOL4},

		/* Bug reported by Koichi Takahashi <ktakahashi@molsci.org>,
		   computed from Pari besseljh(n,x) and AS 10.1.1 */

		{30, 3878.62, -0.00023285567034330878410434732790, TEST_TOL4},
		{49, 9912.63, 5.2043354544842669214485107019e-5, TEST_TOL4},
		{49, 9950.35, 5.0077368819565969286578715503e-5, TEST_TOL4},
		{52, 9930.51, -7.4838588266727718650124475651e-6, TEST_TOL4},

		/* bug report #37209 */
		{364, 36.62, 1.118907148986954e-318, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_jl_e(c.l, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_jl_e(%v, %v)", c.l, c.x))
		})
	}
}
