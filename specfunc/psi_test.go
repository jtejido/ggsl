package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestPsiInt(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n        int
		expected float64
		tol      float64
	}{
		{1, -0.57721566490153286060, TEST_TOL0},
		{2, 0.42278433509846713939, TEST_TOL0},
		{3, 0.92278433509846713939, TEST_TOL0},
		{4, 1.2561176684318004727, TEST_TOL0},
		{5, 1.5061176684318004727, TEST_TOL0},
		{100, 4.600161852738087400, TEST_TOL0},
		{110, 4.695928024251535633, TEST_TOL0},
		{5000, 8.517093188082904107, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Psi_int_e(c.n, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Psi_int_e(%v)", c.n))
		})
	}
}

func TestPsi(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{5000.0, 8.517093188082904107, TEST_TOL0},
		{5.0, 1.5061176684318004727, TEST_TOL0},
		{-10.5, 2.3982391295357816134, TEST_TOL0},
		{-100.5, 4.615124601338064117, TEST_TOL2},
		{-1.0e+5 - 0.5, 11.512935464924395337, 4.0 * TEST_TOL4},
		{-262144.0 - 0.5, 12.476653064769611581, 4.0 * TEST_TOL4},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Psi_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Psi_e(%v)", c.x))
		})
	}
}

func TestPsi1piy(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.8, -0.07088340212750589223, TEST_TOL1},
		{1.0, 0.09465032062247697727, TEST_TOL0},
		{5.0, 1.6127848446157465854, TEST_TOL2},
		{100.0, 4.605178519404762003, TEST_TOL0},
		{2000.0, 7.600902480375416216, TEST_TOL0},
		{-0.8, -0.07088340212750589223, TEST_TOL1},
		{-1.0, 0.09465032062247697727, TEST_TOL0},
		{-5.0, 1.6127848446157465854, TEST_TOL2},
		{-100.0, 4.605178519404762003, TEST_TOL0},
		{-2000.0, 7.600902480375416216, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Psi_1piy_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Psi_1piy_e(%v)", c.x))
		})
	}
}

func TestPsi1int(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n        int
		expected float64
		tol      float64
	}{
		{1, 1.6449340668482264364, TEST_TOL0},
		{2, 0.64493406684822643647, TEST_TOL0},
		{3, 0.39493406684822643647, TEST_TOL0},
		{4, 0.28382295573711532536, TEST_TOL0},
		{1, 1.6449340668482264365, TEST_TOL0},
		{5, 0.22132295573711532536, TEST_TOL0},
		{100, 0.010050166663333571395, TEST_TOL0},
		{110, 0.009132356622022545705, TEST_TOL0},
		{500, 0.0020020013333322666697, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Psi_1_int_e(c.n, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Psi_1_int_e(%v)", c.n))
		})
	}
}

func TestPsi1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0 / 32.0, 1025.5728544782377089, TEST_TOL0},
		{1.0, 1.6449340668482264365, TEST_TOL0},
		{5.0, 0.22132295573711532536, TEST_TOL0},
		{100.0, 0.010050166663333571395, TEST_TOL0},
		{110.0, 0.009132356622022545705, TEST_TOL0},
		{500.0, 0.0020020013333322666697, TEST_TOL0},
		{-1.0 - 1.0/128.0, 16386.648472598746587, TEST_TOL0},
		{-1.50, 9.3792466449891237539, TEST_TOL0},
		{-10.5, 9.7787577398148123845, TEST_TOL0},
		{-15.5, 9.8071247184113896201, TEST_TOL0},
		{-50.5, 9.8499971860824842274, TEST_TOL0},
		{-1000.5, 9.8686054001734414233, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Psi_1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Psi_1_e(%v)", c.x))
		})
	}
}

func TestPsiN(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{1, 1, 1.6449340668482264364, TEST_TOL0},
		{1, 2, 0.64493406684822643647, TEST_TOL0},
		{1, 3, 0.39493406684822643647, TEST_TOL0},
		{1, 4, 0.28382295573711532536, TEST_TOL0},
		{1, 5, 0.22132295573711532536, TEST_TOL0},
		{1, 100, 0.010050166663333571395, TEST_TOL0},
		{1, 110, 0.009132356622022545705, TEST_TOL0},
		{1, 500, 0.0020020013333322666697, TEST_TOL0},
		{3, 5.0, 0.021427828192755075022, TEST_TOL0},
		{3, 500.0, 1.6048063999872000683e-08, TEST_TOL0},
		{10, 5.0, -0.08675107579196581317, TEST_TOL1},
		{10, 50.0, -4.101091112731268288e-12, TEST_TOL0},
		{0, -1.5, 0.70315664064524318723, TEST_TOL0},
		{1, -1.5, 9.3792466449891237539, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Psi_n_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Psi_n_e(%v, %v)", c.n, c.x))

		})
	}
}
