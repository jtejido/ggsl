package specfunc

import (
	"fmt"
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"math"
	"strconv"
	"testing"
)

func TestSin(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10.0, 0.5440211108893698134, TEST_TOL0},
		{1.0, 0.8414709848078965067, TEST_TOL0},
		{1000.0, 0.8268795405320025603, TEST_TOL0},
		{1048576.75, 0.8851545351115651914, TEST_TOL1},
		{62831853.75, 0.6273955953485000827, TEST_TOL3},
		{1073741822.5, -0.8284043541754465988, TEST_SQRT_TOL0},
		{1073741824.0, -0.6173264150460421708, TEST_SQRT_TOL0},
		{1073741825.5, 0.7410684679436226926, TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Sin_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Sin_e(%v)", c.x))
		})
	}
}

func TestCos(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10.0, -0.8390715290764524523, TEST_TOL0},
		{1.0, 0.5403023058681397174, TEST_TOL0},
		{1000.0, 0.5623790762907029911, TEST_TOL1},
		{1048576.75, 0.4652971620066351799, TEST_TOL2},
		{62831853.75, 0.7787006914966116436, TEST_TOL2},
		{1073741822.5, -0.5601305436977716102, TEST_SQRT_TOL0},
		{1073741824.0, 0.7867071229411881196, TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Cos_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Cos_e(%v)", c.x))
		})
	}
}

func TestSinc(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0 / 1024.0, 0.9999984312693665404, TEST_TOL0},
		{1.0 / 2.0, 2.0 / gsl.Pi, TEST_TOL0},
		{80.5, 0.0039541600768172754, TEST_TOL0},
		{100.5, 0.0031672625490924445, TEST_TOL0},
		{1.0e+06 + 0.5, 3.18309727028927157e-07, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Sinc_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Sinc_e(%v)", c.x))
		})
	}
}

func TestComplexSin(t *testing.T) {
	r1, r2 := new(Result), new(Result)
	cases := []struct {
		zr, zi, expected1, tol1, expected2, tol2 float64
	}{
		{1.0, 5.0, 62.44551846769653403, TEST_TOL0, 40.09216577799840254, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Complex_sin_e(c.zr, c.zi, r1, r2)
			run_test_sf_2(t, stat, r1, c.expected1, c.tol1, r2, c.expected2, c.tol2, err.SUCCESS, fmt.Sprintf("Complex_sin_e(%v,%v)", c.zr, c.zi))
		})
	}
}

func TestComplexCos(t *testing.T) {
	r1, r2 := new(Result), new(Result)
	cases := []struct {
		zr, zi, expected1, tol1, expected2, tol2 float64
	}{
		{1.0, 5.0, 40.09580630629882573, TEST_TOL0, -62.43984868079963017, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Complex_cos_e(c.zr, c.zi, r1, r2)
			run_test_sf_2(t, stat, r1, c.expected1, c.tol1, r2, c.expected2, c.tol2, err.SUCCESS, fmt.Sprintf("Complex_cos_e(%v,%v)", c.zr, c.zi))
		})
	}
}

func TestComplexLogSin(t *testing.T) {
	r1, r2 := new(Result), new(Result)
	cases := []struct {
		zr, zi, expected1, tol1, expected2, tol2 float64
	}{
		{1.0, 100.0, 99.3068528194400546900, TEST_TOL0, 0.5707963267948966192, TEST_TOL0},
		{1.0, -100.0, 99.3068528194400546900, TEST_TOL1, -0.5707963267948966192, TEST_TOL1},
		{5.0, 5.0, 4.3068909128079757420, TEST_TOL0, 2.8540063315538773952, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Complex_logsin_e(c.zr, c.zi, r1, r2)
			run_test_sf_2(t, stat, r1, c.expected1, c.tol1, r2, c.expected2, c.tol2, err.SUCCESS, fmt.Sprintf("Complex_logsin_e(%v,%v)", c.zr, c.zi))
		})
	}
}

func TestLnSinh(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.1, -2.3009189815304652235, TEST_TOL0},
		{1.0, 0.16143936157119563361, TEST_TOL0},
		{5.0, 4.306807418479684201, TEST_TOL0},
		{100.0, 99.30685281944005469, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Lnsinh_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Lnsinh_e(%v)", c.x))
		})
	}
}

func TestLnCosh(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.125, 0.007792239318898252791, TEST_TOL0},
		{1.0, 0.4337808304830271870, TEST_TOL0},
		{5.0, 4.306898218339271555, TEST_TOL0},
		{100.0, 99.30685281944005469, TEST_TOL0},
		{1000.0, 999.30685281944005469, TEST_TOL0},

		{-0.125, 0.007792239318898252791, TEST_TOL0},
		{-1.0, 0.4337808304830271870, TEST_TOL0},
		{-5.0, 4.306898218339271555, TEST_TOL0},
		{-100.0, 99.30685281944005469, TEST_TOL0},
		{-1000.0, 999.30685281944005469, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Lncosh_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Lncosh_e(%v)", c.x))
		})
	}
}

func TestPolarToRect(t *testing.T) {
	r1, r2 := new(Result), new(Result)
	cases := []struct {
		r, theta, expected1, tol1, expected2, tol2 float64
	}{
		{10.0, gsl.Pi / 6.0, (10.0 * math.Sqrt(3) / 2.0), TEST_TOL0, (10.0 * 0.5), TEST_TOL0},
		{10.0, -2.0 / 3.0 * gsl.Pi, (10.0 * (-0.5)), TEST_TOL1, (10.0 * (-math.Sqrt(3.0) / 2.0)), TEST_TOL1},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Polar_to_rect(c.r, c.theta, r1, r2)
			run_test_sf_2(t, stat, r1, c.expected1, c.tol1, r2, c.expected2, c.tol2, err.SUCCESS, fmt.Sprintf("Polar_to_rect(%v,%v)", c.r, c.theta))
		})
	}
}

func TestAngleRestrictPos(t *testing.T) {
	DELTA := 1.2246467991473531772e-16

	cases := []struct {
		theta, expected float64
		tol             float64
	}{
		{2.0 * gsl.Pi, 2 * gsl.Pi, TEST_TOL1},
		{-2.0 * gsl.Pi, 2 * DELTA, TEST_TOL1},
		{2.0*gsl.Pi + 4*gsl.Float64Eps, 4*gsl.Float64Eps - 2*DELTA, TEST_TOL1},
		{-2.0*gsl.Pi - 4*gsl.Float64Eps, 2*gsl.Pi - 4*gsl.Float64Eps + 2*DELTA, TEST_TOL1},
		{4.0*gsl.Pi + 8*gsl.Float64Eps, 8*gsl.Float64Eps - 4*DELTA, TEST_TOL1},
		{-4.0*gsl.Pi - 8*gsl.Float64Eps, 2*gsl.Pi - 8*gsl.Float64Eps + 4*DELTA, TEST_TOL1},
		{1e9, 0.5773954235013851694, TEST_TOL1},
		{1e12, 5.625560548042800009446, TEST_SNGL},
		{-1e9, 5.7057898836782013075, TEST_TOL1},
		{-1e12, 0.6576247591367864674792517289, 100 * TEST_SNGL},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			v := c.theta
			Angle_restrict_pos_e(&v)
			run_test_theta(t, v, c.expected, c.tol, fmt.Sprintf("Angle_restrict_pos_e(%v)", c.theta))
		})
	}
}

func TestAngleRestrictPosErr(t *testing.T) {
	DELTA := 1.2246467991473531772e-16
	r := new(Result)
	cases := []struct {
		theta, expected float64
		tol             float64
		status          int
	}{
		{2.0 * gsl.Pi, 2 * gsl.Pi, TEST_TOL1, err.SUCCESS},
		{-2.0 * gsl.Pi, 2 * DELTA, TEST_TOL1, err.SUCCESS},
		{1e9, 0.5773954235013851694, TEST_TOL1, err.SUCCESS},
		{1e12, 5.625560548042800009446, TEST_SNGL, err.SUCCESS},
		{-1e9, 5.7057898836782013075, TEST_TOL1, err.SUCCESS},
		{-1e12, 0.6576247591367864674792517289, 100 * TEST_SNGL, err.SUCCESS},
		{1e15, math.NaN(), TEST_TOL1, err.ELOSS},
		{-1e15, math.NaN(), TEST_TOL1, err.ELOSS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Angle_restrict_pos_err_e(c.theta, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Angle_restrict_pos_err_e(%v)", c.theta))
		})
	}
}

func TestAngleRestrictSymm(t *testing.T) {
	DELTA := 1.2246467991473531772e-16

	cases := []struct {
		theta, expected float64
		tol             float64
	}{
		{2.0 * gsl.Pi, -2 * DELTA, TEST_TOL1},
		{-2.0 * gsl.Pi, 2 * DELTA, TEST_TOL1},
		{gsl.Pi, gsl.Pi, TEST_TOL1},
		{-gsl.Pi, -gsl.Pi, TEST_TOL1},
		{gsl.Pi + 2*gsl.Float64Eps, -gsl.Pi + 2*(gsl.Float64Eps-DELTA), TEST_TOL1},
		{-gsl.Pi - 2*gsl.Float64Eps, gsl.Pi - 2*(gsl.Float64Eps-DELTA), TEST_TOL1},
		{3*gsl.Pi + 6*gsl.Float64Eps, -gsl.Pi + 6*gsl.Float64Eps - 4*DELTA, TEST_TOL1},
		{-3*gsl.Pi - 6*gsl.Float64Eps, gsl.Pi - 6*gsl.Float64Eps + 4*DELTA, TEST_TOL1},
		{1e9, 0.5773954235013851694, TEST_TOL1},
		{1e12, -0.6576247591367864674792517289, 100 * TEST_SNGL},
		{-1e9, -0.5773954235013851694, TEST_TOL1},
		{-1e12, 0.6576247591367864674792517289, 100 * TEST_SNGL},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			v := c.theta
			Angle_restrict_symm_e(&v)
			run_test_theta(t, v, c.expected, c.tol, fmt.Sprintf("Angle_restrict_symm_e(%v)", c.theta))
		})
	}
}

func TestAngleRestrictSymmErr(t *testing.T) {
	DELTA := 1.2246467991473531772e-16
	r := new(Result)
	cases := []struct {
		theta, expected float64
		tol             float64
		status          int
	}{
		{2.0 * gsl.Pi, -2 * DELTA, TEST_TOL1, err.SUCCESS},
		{-2.0 * gsl.Pi, 2 * DELTA, TEST_TOL1, err.SUCCESS},
		{1e9, 0.5773954235013851694, TEST_TOL1, err.SUCCESS},
		{1e12, -0.6576247591367864674792517289, 100 * TEST_SNGL, err.SUCCESS},
		{-1e9, -0.5773954235013851694, TEST_TOL1, err.SUCCESS},
		{-1e12, 0.6576247591367864674792517289, 100 * TEST_SNGL, err.SUCCESS},
		{1e15, math.NaN(), TEST_TOL1, err.ELOSS},
		{-1e15, math.NaN(), TEST_TOL1, err.ELOSS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Angle_restrict_symm_err_e(c.theta, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Angle_restrict_symm_err_e(%v)", c.theta))
		})
	}
}
