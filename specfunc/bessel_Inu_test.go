package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselInuScaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		nu, x, expected float64
		tol             float64
	}{
		{0.0001, 10.0, 0.12783333709581669672, TEST_TOL0},
		{1.0, 0.001, 0.0004995003123542213370, TEST_TOL0},
		{1.0, 1.0, 0.20791041534970844887, TEST_TOL0},
		{30.0, 1.0, 1.3021094983785914437e-42, TEST_TOL0},
		{30.0, 100.0, 0.0004486987756920986146, TEST_TOL3},
		{10.0, 1.0, 1.0127529864692066036e-10, TEST_TOL0},
		{10.0, 100.0, 0.024176682718258828365, TEST_TOL3},
		{10.2, 100.0, 0.023691628843913810043, TEST_TOL3},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Inu_scaled_e(c.nu, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Inu_scaled_e(%v,%v)", c.nu, c.x))

		})
	}
}

func TestBesselInu(t *testing.T) {
	r := new(Result)
	cases := []struct {
		nu, x, expected float64
		tol             float64
	}{
		{0.0001, 10.0, 2815.7166269770030352, TEST_TOL0},
		{1.0, 0.001, 0.0005000000625000026042, TEST_TOL0},
		{1.0, 1.0, 0.5651591039924850272, TEST_TOL0},
		{30.0, 1.0, 3.539500588106447747e-42, TEST_TOL0},
		{30.0, 100.0, 1.2061548704498434006e+40, TEST_TOL2},
		{10.0, 1.0, 2.7529480398368736252e-10, TEST_TOL0},
		{10.0, 100.0, 6.498975524720147799e+41, TEST_TOL2},
		{10.2, 100.0, 6.368587361287030443e+41, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Inu_e(c.nu, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Inu_e(%v,%v)", c.nu, c.x))
		})
	}
}
