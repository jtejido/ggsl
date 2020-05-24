package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"math"
	"strconv"
	"testing"
)

func TestErfc(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{-10.0, 2.0, TEST_TOL0},
		{-5.0000002, 1.9999999999984625433, TEST_TOL0},
		{-5.0, 1.9999999999984625402, TEST_TOL0},
		{-1.0, 1.8427007929497148693, TEST_TOL0},
		{-0.5, 1.5204998778130465377, TEST_TOL0},
		{1.0, 0.15729920705028513066, TEST_TOL0},
		{3.0, 0.000022090496998585441373, TEST_TOL1},
		{7.0, 4.183825607779414399e-23, TEST_TOL2},
		{10.0, 2.0884875837625447570e-45, TEST_TOL2},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Erfc_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Erfc_e(%v)", c.x))

		})
	}
}

func TestLogErfc(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{-1.0, math.Log(1.842700792949714869), TEST_TOL0},
		{-0.1, 0.106576400586522485015, TEST_TOL0},
		{-1e-10, 1.1283791670318505967e-10, TEST_TOL0},
		{0.0, math.Log(1.0), TEST_TOL0},
		{1e-10, -1.128379167159174551e-10, TEST_TOL0},
		{0.001, -0.0011290158896213548027, TEST_TOL0},
		{0.1, -0.119304973737395598329, TEST_TOL0},
		{1.0, math.Log(0.15729920705028513066), TEST_TOL0},
		{10.0, math.Log(2.0884875837625447570e-45), TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Log_erfc_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Log_erfc_e(%v)", c.x))

		})
	}
}

func TestErf(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{-10.0, -1.0000000000000000000, TEST_TOL0},
		{0.5, 0.5204998778130465377, TEST_TOL0},
		{1.0, 0.8427007929497148693, TEST_TOL0},
		{10.0, 1.0000000000000000000, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Erf_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Erf_e(%v)", c.x))

		})
	}
}

func TestErfZ(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{1.0, 0.24197072451914334980, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Erf_Z_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Erf_Z_e(%v)", c.x))

		})
	}
}

func TestErfQ(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{10.0, 7.619853024160526066e-24, TEST_TOL2},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Erf_Q_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Erf_Q_e(%v)", c.x))

		})
	}
}

func TestHazard(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{-20.0, 5.5209483621597631896e-88, TEST_TOL2},
		{-10.0, 7.6945986267064193463e-23, TEST_TOL2},
		{-1.0, 0.28759997093917836123, TEST_TOL0},
		{0.0, 0.79788456080286535588, TEST_TOL0},
		{1.0, 1.5251352761609812091, TEST_TOL0},
		{10.0, 10.098093233962511963, TEST_TOL2},
		{20.0, 20.049753068527850542, TEST_TOL2},
		{30.0, 30.033259667433677037, TEST_TOL2},
		{50.0, 50.019984031905639809, TEST_TOL0},
		{80.0, 80.012496096798234468, TEST_TOL0},
		{150.0, 150.00666607420571802, TEST_TOL0},
		{300.0, 300.00333325926337415, TEST_TOL0},
		{900.0, 900.00111110836764382, TEST_TOL0},
		{1001.0, 1001.0009989990049990, TEST_TOL0},
		{2000.0, 2000.0004999997500003, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Hazard_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Hazard_e(%v)", c.x))

		})
	}
}
