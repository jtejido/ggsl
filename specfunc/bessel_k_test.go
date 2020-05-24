package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselk0Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 15.707963267948966192, TEST_TOL0},
		{2.0, 0.7853981633974483096, TEST_TOL0},
		{100.0, 0.015707963267948966192, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_k0_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_k0_scaled_e(%v)", c.x))
		})
	}
}

func TestBesselk1Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 172.78759594743862812, TEST_TOL0},
		{2.0, 1.1780972450961724644, TEST_TOL0},
		{100.0, 0.015865042900628455854, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_k1_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_k1_scaled_e(%v)", c.x))
		})
	}
}

func TestBesselk2Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 5199.335841691107810, TEST_TOL0},
		{2.0, 2.5525440310417070063, TEST_TOL0},
		{100.0, 0.016183914554967819868, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_k2_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_k2_scaled_e(%v)", c.x))
		})
	}
}

func TestBesselklScaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l           int
		x, expected float64
		tol         float64
	}{
		{4, 1.0 / 256.0, 1.8205599816961954439e+14, TEST_TOL0},
		{4, 1.0 / 8.0, 6.1173217814406597530e+06, TEST_TOL0},
		{5, 2.0, 138.10735829492005119, TEST_TOL0},
		{100, 100.0, 3.985930768060258219e+18, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_kl_scaled_e(c.l, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_kl_scaled_e(%v, %v)", c.l, c.x))
		})
	}
}
