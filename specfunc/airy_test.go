package specfunc

import (
	"fmt"
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestAiryAi(t *testing.T) {
	mode := gsl.MODE_DEFAULT
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-500.0, 0.0725901201040411396, TEST_TOL4},
		{-5.0, 0.3507610090241142, TEST_TOL0},
		{-0.3000000000000094, 0.4309030952855831, TEST_TOL0},
		{0.6999999999999907, 0.1891624003981519, TEST_TOL0},
		{1.649999999999991, 0.0583105861872088521, TEST_TOL0},
		{2.54999999999999, 0.01446149513295428, TEST_TOL0},
		{3.499999999999987, 0.002584098786989702, TEST_TOL1},
		{5.39999999999998, 4.272986169411866e-05, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_Ai_e(c.x, mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_Ai_e(%v)", c.x))

		})
	}
}

func TestAiryAiScaled(t *testing.T) {
	mode := gsl.MODE_DEFAULT
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-5.0, 0.3507610090241142, TEST_TOL0},
		{0.6999999999999907, 0.2795125667681217, TEST_TOL0},
		{1.649999999999991, 0.2395493001442741, TEST_TOL0},
		{2.54999999999999, 0.2183658595899388, TEST_TOL0},
		{3.499999999999987, 0.2032920808163519, TEST_TOL0},
		{5.39999999999998, 0.1836050093282229, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_Ai_scaled_e(c.x, mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_Ai_scaled_e(%v)", c.x))
		})
	}
}

func TestAiryBi(t *testing.T) {
	mode := gsl.MODE_DEFAULT
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-500.0, -0.094688570132991028, TEST_TOL4},
		{-5.0, -0.1383691349016005, TEST_TOL1},
		{0.6999999999999907, 0.9733286558781599, TEST_TOL0},
		{1.649999999999991, 2.196407956850028, TEST_TOL0},
		{2.54999999999999, 6.973628612493443, TEST_TOL0},
		{3.499999999999987, 33.05550675461069, TEST_TOL1},
		{5.39999999999998, 1604.476078241272, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_Bi_e(c.x, mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_Bi_e(%v)", c.x))
		})
	}
}

func TestAiryBiScaled(t *testing.T) {
	mode := gsl.MODE_DEFAULT
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-5.0, -0.1383691349016005, TEST_TOL1},
		{0.6999999999999907, 0.6587080754582302, TEST_TOL0},
		{1.649999999999991, 0.5346449995597539, TEST_TOL0},
		{2.54999999999999, 0.461835455542297, TEST_TOL0},
		{3.499999999999987, 0.4201771882353061, TEST_TOL1},
		{5.39999999999998, 0.3734050675720473, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Airy_Bi_scaled_e(c.x, mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Airy_Bi_scaled_e(%v)", c.x))
		})
	}
}
