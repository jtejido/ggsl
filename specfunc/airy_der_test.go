package specfunc

import (
	"fmt"
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestAiryAiDeriv(t *testing.T) {
	mode := gsl.MODE_DEFAULT
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-5.0, 0.3271928185544435, TEST_TOL1},
		{-0.5500000000000094, -0.1914604987143629, TEST_TOL0},
		{0.4999999999999906, -0.2249105326646850, TEST_TOL0},
		{1.899999999999992, -0.06043678178575718, TEST_TOL0},
		{3.249999999999988, -0.007792687926790889, TEST_TOL0},
		{5.199999999999981, -0.0001589434526459543, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_Ai_deriv_e(c.x, mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_Ai_deriv_e(%v)", c.x))
		})
	}
}

func TestAiryAiDerivScaled(t *testing.T) {
	mode := gsl.MODE_DEFAULT
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-5.0, 0.3271928185544435, TEST_TOL1},
		{0.5499999999999906, -0.2874057279170166, TEST_TOL0},
		{1.499999999999991, -0.3314199796863637, TEST_TOL0},
		{2.49999999999999, -0.3661089384751620, TEST_TOL0},
		{3.649999999999986, -0.3974033831453963, TEST_TOL0},
		{6.299999999999977, -0.4508799189585947, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_Ai_deriv_scaled_e(c.x, mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_Ai_deriv_scaled_e(%v)", c.x))

		})
	}
}

func TestAiryBiDeriv(t *testing.T) {
	mode := gsl.MODE_DEFAULT
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-5.0, 0.778411773001899, TEST_TOL0},
		{-0.5500000000000094, 0.5155785358765014, TEST_TOL0},
		{0.4999999999999906, 0.5445725641405883, TEST_TOL0},
		{1.899999999999992, 3.495165862891568, TEST_TOL0},
		{3.249999999999988, 36.55485149250338, TEST_TOL0},
		{5.199999999999981, 2279.748293583233, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_Bi_deriv_e(c.x, mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_Bi_deriv_e(%v)", c.x))

		})
	}
}

func TestAiryBiDerivScaled(t *testing.T) {
	mode := gsl.MODE_DEFAULT
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-5.0, 0.778411773001899, TEST_TOL0},
		{0.5499999999999906, 0.4322811281817566, TEST_TOL0},
		{1.499999999999991, 0.5542307563918037, TEST_TOL0},
		{2.49999999999999, 0.6755384441644985, TEST_TOL0},
		{3.649999999999986, 0.7613959373000228, TEST_TOL0},
		{6.299999999999977, 0.8852064139737571, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_Bi_deriv_scaled_e(c.x, mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_Bi_deriv_scaled_e(%v)", c.x))

		})
	}
}
