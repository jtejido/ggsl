package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestZeta(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-151, 8.195215221831378294e+143, TEST_TOL2},
		{-51, 9.68995788746359406565e+24, TEST_TOL1},
		{-5, -0.003968253968253968253968, TEST_TOL1},
		{-8, 0.0, TEST_TOL1},
		{-6, 0.0, TEST_TOL1},
		{-4, 0.0, TEST_TOL1},
		{-3, 1.0 / 120.0, TEST_TOL1},
		{-2, 0.0, TEST_TOL1},
		{-1, -1.0 / 12.0, TEST_TOL1},
		{-0.5, -0.207886224977354566017307, TEST_TOL1},
		{-1e-10, -0.49999999990810614668948, TEST_TOL1},
		{0.0, -0.5, TEST_TOL0},
		{1e-10, -0.50000000009189385333058, TEST_TOL0},
		{0.5, -1.460354508809586812889499, TEST_TOL0},
		{1.0 - 1.0/1024.0, -1023.4228554489429787, TEST_TOL0},
		{1.0 + 1.0/1048576, 1.0485765772157343441e+06, TEST_TOL0},
		{5.0, 1.036927755143369926331365, TEST_TOL0},
		{25.5, 1.000000021074106110269959, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Zeta_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Zeta_e(%v)", c.x))
		})
	}
}

func TestZetaInt(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x        int
		expected float64
		tol      float64
	}{
		{-61.0, -3.30660898765775767257e+34, TEST_TOL0},
		{-8, 0.0, TEST_TOL0},
		{-6, 0.0, TEST_TOL0},
		{-5.0, -0.003968253968253968253968, TEST_TOL0},
		{-4, 0.0, TEST_TOL0},
		{-3, 1.0 / 120.0, TEST_TOL0},
		{-2, 0.0, TEST_TOL0},
		{-1, -1.0 / 12.0, TEST_TOL0},
		{5.0, 1.0369277551433699263313655, TEST_TOL0},
		{31.0, 1.0000000004656629065033784, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Zeta_int_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Zeta_int_e(%v)", c.x))
		})
	}
}

func TestZetaM1Int(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x        int
		expected float64
		tol      float64
	}{
		{-61.0, -3.30660898765775767257e+34, TEST_TOL0},
		{-5.0, -1.003968253968253968253968, TEST_TOL0},
		{-8, -1.0, TEST_TOL0},
		{-6, -1.0, TEST_TOL0},
		{-4, -1.0, TEST_TOL0},
		{-3, -119.0 / 120.0, TEST_TOL0},
		{-2, -1.0, TEST_TOL0},
		{-1, -13.0 / 12.0, TEST_TOL0},
		{5.0, 0.0369277551433699263313655, TEST_TOL0},
		{31.0, 0.0000000004656629065033784, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Zetam1_int_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Zetam1_int_e(%v)", c.x))
		})
	}
}

func TestZetaM1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-8, -1.0, TEST_TOL1},
		{-6, -1.0, TEST_TOL1},
		{-4, -1.0, TEST_TOL1},
		{-3, -119.0 / 120.0, TEST_TOL1},
		{-2, -1.0, TEST_TOL1},
		{-1, -13.0 / 12.0, TEST_TOL1},
		{-0.5, -1.207886224977354566017307, TEST_TOL1},
		{-1e-10, -1.49999999990810614668948, TEST_TOL1},
		{0.0, -1.5, TEST_TOL0},
		{1e-10, -1.50000000009189385333058, TEST_TOL0},
		{0.5, -2.460354508809586812889499, TEST_TOL0},
		{2.0, 0.64493406684822643647, TEST_TOL1},
		{3.0, 0.20205690315959428540, TEST_TOL1},
		{5.0, 0.0369277551433699263314, TEST_TOL1},
		{9.5, 0.0014125906121736622712, TEST_TOL1},
		{10.5, 0.000700842641736155219500, TEST_TOL1},
		{12.5, 0.000173751733643178193390, TEST_TOL1},
		{13.5, 0.000086686727462338155188, TEST_TOL1},
		{15.5, 0.000021619904246069108133, TEST_TOL1},
		{16.5, 0.000010803124900178547671, TEST_TOL0},
		{25.5, 0.000000021074106110269959, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Zetam1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Zetam1_e(%v)", c.x))
		})
	}
}

func TestEtaInt(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x        int
		expected float64
		tol      float64
	}{
		{-91, -4.945598888750002040e+94, TEST_TOL0},
		{-51, -4.363969073121683116e+40, TEST_TOL0},
		{-5, 0.25, TEST_TOL0},
		{-1, 0.25, TEST_TOL0},
		{0, 0.5, TEST_TOL0},
		{5, 0.9721197704469093059, TEST_TOL0},
		{6, 0.9855510912974351041, TEST_TOL0},
		{20, 0.9999990466115815221, TEST_TOL0},
		{1000, 1.0, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Eta_int_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Eta_int_e(%v)", c.x))
		})
	}
}

func TestEta(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-51.5, -1.2524184036924703656e+41, TEST_TOL2},
		{-5, 0.25, TEST_TOL0},
		{0.5, 0.6048986434216303702, TEST_TOL0},
		{0.999, 0.6929872789683383574, TEST_TOL0},
		{1.0, 0.6931471805599453094, TEST_TOL0},
		{1.0 + 1.0e-10, 0.6931471805759321998, TEST_TOL0},
		{5, 0.9721197704469093059, TEST_TOL0},
		{5.2, 0.9755278712546684682, TEST_TOL0},
		{6, 0.9855510912974351041, TEST_TOL0},
		{20, 0.9999990466115815221, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Eta_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Eta_e(%v)", c.x))

		})
	}
}
