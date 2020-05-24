package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestBesselZeroJ0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		s        uint
		expected float64
		tol      float64
		status   int
	}{
		{0, 0, TEST_TOL0, err.EINVAL},
		{1, 2.404825557695771, TEST_TOL1, err.SUCCESS},
		{2, 5.520078110286304, TEST_TOL1, err.SUCCESS},
		{20, 62.048469190227081, TEST_TOL1, err.SUCCESS},
		{25, 77.756025630388058, TEST_TOL1, err.SUCCESS},
		{100, 313.37426607752784, TEST_TOL1, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_zero_J0_e(c.s, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Bessel_zero_J0_e(%v)", c.s))
		})
	}
}

func TestBesselZeroJ1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		s        uint
		expected float64
		tol      float64
		status   int
	}{
		{0, 0, TEST_TOL0, err.SUCCESS},
		{1, 3.831705970207512, TEST_TOL2, err.SUCCESS},
		{2, 7.015586669815619, TEST_TOL2, err.SUCCESS},
		{20, 63.61135669848124, TEST_TOL2, err.SUCCESS},
		{25, 79.32048717547630, TEST_TOL2, err.SUCCESS},
		{100, 314.9434728377672, TEST_TOL2, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_zero_J1_e(c.s, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Bessel_zero_J1_e(%v)", c.s))
		})
	}
}

func TestBesselZeroJnu(t *testing.T) {
	r := new(Result)
	cases := []struct {
		nu       float64
		s        uint
		expected float64
		tol      float64
		status   int
	}{
		{0.0, 0, 0, TEST_TOL0, err.EINVAL},
		{0.0, 1, 2.404825557695771, TEST_TOL1, err.SUCCESS},
		{0.0, 2, 5.520078110286304, TEST_TOL1, err.SUCCESS},
		{0.0, 20, 62.048469190227081, TEST_TOL1, err.SUCCESS},
		{0.0, 25, 77.756025630388058, TEST_TOL1, err.SUCCESS},
		{0.0, 100, 313.37426607752784, TEST_TOL1, err.SUCCESS},
		{1.0, 0, 0, TEST_TOL0, err.SUCCESS},
		{1.5, 1, 4.4934094579090641, TEST_TOL1, err.SUCCESS},
		{5.0, 1, 8.7714838159599540, TEST_TOL1, err.SUCCESS},
		{1.5, 2, 7.7252518369377072, TEST_TOL1, err.SUCCESS},
		{5.0, 2, 12.338604197466944, TEST_TOL1, err.SUCCESS},
		{1.5, 3, 10.904121659428900, TEST_TOL1, err.SUCCESS},
		{5.0, 3, 15.700174079711671, TEST_TOL1, err.SUCCESS},
		{1.5, 4, 14.066193912831473, TEST_TOL1, err.SUCCESS},
		{5.0, 4, 18.980133875179921, TEST_TOL1, err.SUCCESS},
		{1.5, 5, 17.220755271930768, TEST_TOL1, err.SUCCESS},
		/* Something wrong with the tolerances on these */
		{5.0, 5, 22.217799896561268, TEST_SQRT_TOL0, err.SUCCESS},
		{8.0, 5, 26.266814641176644, TEST_SQRT_TOL0, err.SUCCESS},
		{20.0, 5, 41.413065513892636, TEST_SQRT_TOL0, err.SUCCESS},
		{1.5, 6, 20.371302959287563, TEST_TOL1, err.SUCCESS},
		{5.0, 6, 25.430341154222704, TEST_TOL1, err.SUCCESS},
		{8.0, 6, 29.545659670998550, TEST_TOL1, err.SUCCESS},
		{1.5, 7, 23.519452498689007, TEST_TOL1, err.SUCCESS},
		{5.0, 7, 28.626618307291138, TEST_TOL1, err.SUCCESS},
		{8.0, 7, 32.795800037341462, TEST_TOL1, err.SUCCESS},
		{1.5, 8, 26.666054258812674, TEST_TOL1, err.SUCCESS},
		{5.0, 8, 31.811716724047763, TEST_TOL1, err.SUCCESS},
		{10.0, 8, 38.761807017881651, TEST_TOL1, err.SUCCESS},
		{1.5, 9, 29.811598790892959, TEST_TOL1, err.SUCCESS},
		{5.0, 9, 34.988781294559295, TEST_TOL1, err.SUCCESS},
		{10.0, 9, 42.004190236671805, TEST_TOL1, err.SUCCESS},
		{1.5, 10, 32.956389039822477, TEST_TOL1, err.SUCCESS},
		{5.0, 10, 38.159868561967132, TEST_TOL1, err.SUCCESS},
		{15.0, 10, 52.017241278881633, TEST_TOL1, err.SUCCESS},
		{5.0, 11, 41.326383254047406, TEST_TOL1, err.SUCCESS},
		{15.0, 11, 55.289204146560061, TEST_TOL1, err.SUCCESS},
		{5.0, 12, 44.4893191232197314, TEST_TOL1, err.SUCCESS},
		{15.0, 12, 58.5458289043850856, TEST_TOL1, err.SUCCESS},
		{5.0, 13, 47.6493998066970948, TEST_TOL1, err.SUCCESS},
		{15.0, 13, 61.7897598959450550, TEST_TOL1, err.SUCCESS},
		{5.0, 14, 50.8071652030063595, TEST_TOL1, err.SUCCESS},
		{15.0, 14, 65.0230502510422545, TEST_TOL1, err.SUCCESS},
		{5.0, 15, 53.9630265583781707, TEST_TOL1, err.SUCCESS},
		{15.0, 15, 68.2473219964207837, TEST_TOL1, err.SUCCESS},
		{5.0, 16, 57.1173027815042647, TEST_TOL1, err.SUCCESS},
		{15.0, 16, 71.4638758850226630, TEST_TOL1, err.SUCCESS},
		{5.0, 17, 60.2702450729428077, TEST_TOL1, err.SUCCESS},
		{15.0, 17, 74.6737687121404241, TEST_TOL1, err.SUCCESS},
		{5.0, 18, 63.4220540458757799, TEST_TOL1, err.SUCCESS},
		{15.0, 18, 77.8778689734863729, TEST_TOL1, err.SUCCESS},
		{5.0, 19, 66.5728918871182703, TEST_TOL1, err.SUCCESS},
		{15.0, 19, 81.0768977206328326, TEST_TOL1, err.SUCCESS},
		{5.0, 20, 69.722891161716742, TEST_TOL1, err.SUCCESS},
		{15.0, 20, 84.271459069716442, TEST_TOL1, err.SUCCESS},
		{23.0, 11, 65.843393469524653, TEST_TOL6, err.SUCCESS},
		{30.0, 11, 74.797306585175426, TEST_TOL6, err.SUCCESS},
		{32.0, 15, 90.913637691861741, TEST_TOL6, err.SUCCESS},
		{50.0, 15, 113.69747988073942, TEST_TOL6, err.SUCCESS},
		{5.0, 22, 76.020793430591605, TEST_TOL2, err.SUCCESS},
		{10.0, 22, 83.439189796105756, TEST_TOL3, err.SUCCESS},
		{12.0, 22, 86.345496520534055, TEST_TOL6, err.SUCCESS},
		{100.0, 22, 199.82150220122519, TEST_TOL4, err.SUCCESS},
		{500.0, 22, 649.34132440891735, TEST_TOL2, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Bessel_zero_Jnu_e(c.nu, c.s, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Bessel_zero_Jnu_e(%v, %v)", c.nu, c.s))
		})
	}
}
