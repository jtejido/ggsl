package specfunc

import (
	"fmt"
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"math"
	"strconv"
	"testing"
)

func TestLambertW0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
		status           int
	}{
		{0.0, 0.0, TEST_TOL0, err.SUCCESS},
		{1.0, 0.567143290409783872999969, TEST_TOL0, err.SUCCESS},
		{2.0, 0.852605502013725491346472, TEST_TOL0, err.SUCCESS},
		{20.0, 2.205003278024059970493066, TEST_TOL0, err.SUCCESS},
		{1000.0, 5.24960285240159622712606, TEST_TOL0, err.SUCCESS},
		{1.0e+6, 11.38335808614005262200016, TEST_TOL0, err.SUCCESS},
		{1.0e+12, 24.43500440493491313826305, TEST_TOL0, err.SUCCESS},
		{1.0e+308, 702.641362034106812081125, TEST_TOL0, err.SUCCESS},

		/* Test case from Katrin Wolff <katrin_wolff@gmx.de> fails under
		   double-precision */

		{1.6849341956993852953416990, 0.775706963944252869680440, TEST_TOL0, err.SUCCESS},
		{-1.0/math.E - gsl.Float64Eps, -1.0, TEST_TOL0, err.EDOM},
		{-1.0/math.E + 1.0/(1024.0*1024.0*1024.0), -0.999928845560308370714970, TEST_TOL0, err.SUCCESS},
		{-1.0/math.E + 1.0/(1024.0*1024.0), -0.997724730359774141620354, TEST_TOL0, err.SUCCESS},
		{-1.0/math.E + 1.0/512.0, -0.900335676696088773044678, TEST_TOL0, err.SUCCESS},
		{-1.0/math.E + 0.25, -0.1349044682661213545487599, TEST_TOL0, err.SUCCESS},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Lambert_W0_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Lambert_W0_e(%v)", c.x))

		})
	}
}

func TestLambertWm1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
		status           int
	}{
		{0.0, 0.0, TEST_TOL0, err.SUCCESS},
		{1.0, 0.567143290409783872999969, TEST_TOL0, err.SUCCESS},
		{2.0, 0.852605502013725491346472, TEST_TOL0, err.SUCCESS},
		{20.0, 2.205003278024059970493066, TEST_TOL0, err.SUCCESS},

		{-1.0/math.E - gsl.Float64Eps, -1.0, TEST_TOL0, err.EDOM},
		{-1.0/math.E + 1.0/(1024.0*1024.0*1024.0), -1.000071157815154608049055, TEST_TOL1, err.SUCCESS},
		{-1.0/math.E + 1.0/(1024.0*1024.0), -1.002278726118593023934693, TEST_TOL1, err.SUCCESS},
		{-1.0/math.E + 1.0/512.0, -1.106761200865743124599130, TEST_TOL1, err.SUCCESS},
		{-1.0/math.E + 1.0/64.0, -1.324240940341812125489772, TEST_TOL1, err.SUCCESS},
		{-1.0/math.E + 0.25, -3.345798131120112, TEST_TOL1, err.SUCCESS},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Lambert_Wm1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Lambert_Wm1_e(%v)", c.x))

		})
	}
}
