package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselK1Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 10.890182683049696574, TEST_TOL0},
		{1.95, 1.050086915104152747182, TEST_TOL0},
		{2.0, 1.0334768470686885732, TEST_TOL0},
		{6.0, 0.5421759102771335382849, TEST_TOL0},
		{100.0, 0.1257999504795785293, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_K1_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_K1_scaled_e(%v)", c.x))
		})
	}
}

func TestBesselK1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, 9.853844780870606135, TEST_TOL0},
		{1.95, 0.1494001409315894276793, TEST_TOL0},
		{2.0, 0.13986588181652242728, TEST_TOL0},
		{100.0, 4.679853735636909287e-45, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_K1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_K1_e(%v)", c.x))
		})
	}
}
