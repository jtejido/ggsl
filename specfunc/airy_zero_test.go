package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestAiry0Ai(t *testing.T) {
	r := new(Result)
	cases := []struct {
		s        int
		expected float64
		tol      float64
	}{
		{2, -4.087949444130970617, TEST_TOL0},
		{50, -38.02100867725525443, TEST_TOL0},
		{100, -60.45555727411669871, TEST_TOL0},
		{110, -64.43135670991324811, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_zero_Ai_e(c.s, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_zero_Ai_e(%v)", c.s))
		})
	}
}

func TestAiry0Bi(t *testing.T) {
	r := new(Result)
	cases := []struct {
		s        int
		expected float64
		tol      float64
	}{
		{2, -3.271093302836352716, TEST_TOL0},
		{50, -37.76583438165180116, TEST_TOL0},
		{100, -60.25336482580837088, TEST_TOL0},
		{110, -64.2355167606561537, TEST_TOL0},
		{111, -64.6268994819519378, TEST_TOL0},
		{200, -95.88699147356682665, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_zero_Bi_e(c.s, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_zero_Bi_e(%v)", c.s))
		})
	}
}

func TestAiry0AiDeriv(t *testing.T) {
	r := new(Result)
	cases := []struct {
		s        int
		expected float64
		tol      float64
	}{
		{2, -3.248197582179836561, TEST_TOL0},
		{50, -37.76565910053887108, TEST_TOL0},
		{100, -60.25329596442479317, TEST_TOL0},
		{110, -64.23545617243546956, TEST_TOL0},
		{1000, -280.9378080358935071, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_zero_Ai_deriv_e(c.s, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_zero_Ai_deriv_e(%v)", c.s))
		})
	}
}

func TestAiry0BiDeriv(t *testing.T) {
	r := new(Result)
	cases := []struct {
		s        int
		expected float64
		tol      float64
	}{
		{2, -4.073155089071828216, TEST_TOL0},
		{50, -38.02083574095788210, TEST_TOL0},
		{100, -60.45548887257140819, TEST_TOL0},
		{110, -64.43129648944845060, TEST_TOL0},
		{111, -64.82208737584206093, TEST_TOL0},
		{200, -96.04731050310324450, TEST_TOL0},
		{1000, -281.0315164471118527, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_zero_Bi_deriv_e(c.s, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_zero_Bi_deriv_e(%v)", c.s))
		})
	}
}
