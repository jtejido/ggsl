package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselKnuScaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		nu, x, expected float64
		tol             float64
	}{
		// {0.0001, 10.0, 0.3916319346235421817, TEST_TOL0},
		// {1.0, 0.001, 1000.9967345590684524, TEST_TOL0},
		// {1.0, 1.0, 1.6361534862632582465, TEST_TOL0},
		// {30.0, 1.0, 1.2792629867539753925e+40, TEST_TOL0},
		// {30.0, 100.0, 10.673443449954850040, TEST_TOL0},
		// {10.0, 1.0, 4.912296520990198599e+08, TEST_TOL0},
		// {10.0, 100.0, 0.20578687173955779807, TEST_TOL0},
		// {10.0, 1000.0, 0.04165905142800565788, TEST_TOL0},
		// {10.0, 1.0e+8, 0.00012533147624060789938, TEST_TOL0},
		// {10.2, 100.0, 0.20995808355244385075, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Knu_scaled_e(c.nu, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Knu_scaled_e(%v,%v)", c.nu, c.x))
		})
	}
}

func TestBesselKnu(t *testing.T) {
	r := new(Result)
	cases := []struct {
		nu, x, expected float64
		tol             float64
	}{
		// {0.0001, 0.001, 7.023689431812884141, TEST_TOL0},
		// {0.0001, 10.0, 0.000017780062324654874306, TEST_TOL0},
		// {1.0, 0.001, 999.9962381560855743, TEST_TOL0},
		// {1.0, 1.0, 0.6019072301972345747, TEST_TOL0},
		// {10.0, 0.001, 1.8579455483904008064e+38, TEST_TOL0},
		// {10.0, 1.0, 1.8071328990102945469e+08, TEST_TOL0},
		// {10.0, 100.0, 7.655427977388100611e-45, TEST_TOL2},
		// {10.2, 100.0, 7.810600225948217841e-45, TEST_TOL2},
		// {30.0, 1.0, 4.706145526783626883e+39, TEST_TOL1},
		// {30.0, 100.0, 3.970602055959398739e-43, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_Knu_e(c.nu, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_Knu_e(%v,%v)", c.nu, c.x))
		})
	}
}

func TestBesselLnKnu(t *testing.T) {
	r := new(Result)
	cases := []struct {
		nu, x, expected float64
		tol             float64
	}{
		// {0.0001, 1.0e-100, 5.439794449319847, TEST_TOL0},
		// {0.0001, 0.0001, 2.232835507214331, TEST_TOL0},
		// {0.0001, 10.0, -10.93743282256098, TEST_TOL0},
		// {1.0, 1.0e-100, 230.2585092994045, TEST_TOL0},
		// {1.0, 1.0e-10, 23.025850929940456840, TEST_TOL0},
		// {1.0, 0.001, 6.907751517131146, TEST_TOL0},
		// {1.0, 1.0, -0.5076519482107523309, TEST_TOL0},
		// {30.0, 1.0e-100, 6999.113586185543475, TEST_TOL0},
		// {30.0, 1.0, 91.34968784026325464, TEST_TOL0},
		// {30.0, 100.0, -97.63224126416760932, TEST_TOL0},
		// {100.0, 1.0e-100, 23453.606706185466825, TEST_TOL0},
		// {100.0, 1.0, 427.7532510250188083, TEST_TOL0},
		// {100.0, 100.0, -55.53422771502921431, TEST_TOL0},
		// {1000.0, 1.0e-100, 236856.183755993135, TEST_TOL0},
		// {10000.0, 1.0e-100, 2.39161558914890695e+06, TEST_TOL0},
		// /* [bug #31528] gsl_sf_bessel_lnKnu overflows for large nu */
		// {180.0, 2.2, 735.1994170369583930752590258, TEST_TOL1},
		// {3500.5, 1500.0, 1731.220077116482710070986699, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_lnKnu_e(c.nu, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Bessel_lnKnu_e(%v,%v)", c.nu, c.x))
		})
	}
}
