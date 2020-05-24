package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestLegendreQ0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		/* x = -1 + 2^-16 */
		{-0.9999847412109375, -5.8917472200477175158028143531855, TEST_TOL4},
		{-0.5, -0.5493061443340548457, TEST_TOL0},
		{-1e-10, -1.000000000000000000e-10, TEST_TOL0},
		{0.0, 0.0, TEST_TOL0},
		{1e-10, 1.000000000000000000e-10, TEST_TOL0},
		/* x = 1 - 2^-16 */
		{0.9999847412109375, 5.8917472200477175158028143531855, TEST_TOL4},
		/* x = 1 + 2^-16 */
		{1.0000152587890625, 5.8917548494422489138325509750429, TEST_TOL4},
		{1.5, 0.8047189562170501873, TEST_TOL0},
		{9.99, 0.1004364599660005447, TEST_TOL0},
		{10.0, 0.1003353477310755806, TEST_TOL0},
		{10.01, 0.1002344395571710243, TEST_TOL0},
		{100, 0.010000333353334762015, TEST_TOL0},
		{1e10, 1.000000000000000000e-10, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_Q0_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_Q0_e(%v)", c.x))

		})
	}
}

func TestLegendreQ1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{-0.9999847412109375, 4.8916573191196772369, TEST_TOL4},
		{-0.5, -0.7253469278329725772, TEST_TOL0},
		{-0.01, -0.9998999966664666524, TEST_TOL0},
		{-1e-10, -0.999999999999999999, TEST_TOL0},
		{0.0, -1.0, TEST_TOL0},
		{1e-10, -0.999999999999999999, TEST_TOL0},
		{0.0001, -0.9999999899999999667, TEST_TOL0},
		{0.01, -0.9998999966664666524, TEST_TOL0},
		{0.5, -0.7253469278329725772, TEST_TOL0},
		{0.9999847412109375, 4.8916573191196772369, TEST_TOL4},
		{1.0000152587890625, 4.8918447504867045145, TEST_TOL4},
		{1.5, 0.20707843432557528095, TEST_TOL0},
		{9.99, 3.360235060345441639e-3, TEST_TOL0},
		{10.0, 3.353477310755806357e-3, TEST_TOL0},
		{10.01, 3.346739967281953346e-3, TEST_TOL0},
		{100.0, 3.333533347620158821e-5, TEST_TOL0},
		{1e10, 3.333333333333333333e-21, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_Q1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_Q1_e(%v)", c.x))

		})
	}
}

func TestLegendreQl(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l                int
		x, expected, tol float64
	}{
		{10, -0.5, -0.29165813966586752393, TEST_TOL0},
		{10, 0.5, 0.29165813966586752393, TEST_TOL0},
		{10, 1.5, 0.000014714232718207477406, TEST_TOL0},

		{100, -0.5, -0.09492507395207282096, TEST_TOL1},
		{100, 0.5, 0.09492507395207282096, TEST_TOL1},
		{100, 1.5, 1.1628163435044121988e-43, TEST_TOL2},

		{1000, -0.5, -0.030105074974005303500, TEST_TOL1},
		{1000, 0.5, 0.030105074974005303500, TEST_TOL1},
		{1000, 1.1, 1.0757258447825356443e-194, TEST_TOL3},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Legendre_Ql_e(c.l, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Legendre_Ql_e(%v,%v)", c.l, c.x))

		})
	}
}
