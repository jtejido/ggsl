package specfunc

import (
	"fmt"
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"github.com/lucky-se7en/ggsl/test"
	"math"
	"strconv"
	"testing"
)

func TestExp(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10.0, math.Exp(-10.0), TEST_TOL0},
		{10.0, math.Exp(10.0), TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exp_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Exp_e(%v)", c.x))
		})
	}
}

func TestExpE10(t *testing.T) {
	re := new(Result_e10)
	sa := 0
	stat := Exp_e10_e(1.0, re)
	if stat != nil {
		sa += stat.(err.GSLError).Status()
	}
	if test_sf_frac_diff(re.val, math.E) > TEST_TOL0 {
		sa++
	}
	if re.err > TEST_TOL1 {
		sa++
	}
	if re.e10 != 0 {
		sa++
	}

	test.Test(t, sa, "Exp_e10_e(1.0)")

	sa = 0
	stat = Exp_e10_e(2000.0, re)
	if stat != nil {
		sa += stat.(err.GSLError).Status()
	}
	if test_sf_frac_diff(re.val, 3.88118019428363725) > TEST_TOL3 {
		sa++
	}
	if re.err > TEST_TOL5 {
		sa++
	}
	if re.e10 != 868 {
		sa++
	}
	test.Test(t, sa, "Exp_e10_e(2000.0)")
}

func TestExpErr(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, dx, expected float64
		tol             float64
	}{
		{-10.0, TEST_TOL1, math.Exp(-10.0), TEST_TOL1},
		{10.0, TEST_TOL1, math.Exp(10.0), TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exp_err_e(c.x, c.dx, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Exp_err_e(%v, %v)", c.x, c.dx))
		})
	}
}

func TestExpErrE10(t *testing.T) {
	re := new(Result_e10)
	sa := 0
	stat := Exp_err_e10_e(1.0, TEST_SQRT_TOL0, re)
	if stat != nil {
		sa += stat.(err.GSLError).Status()
	}

	if test_sf_frac_diff(re.val, math.E) > TEST_TOL1 {
		sa++
	}
	if re.err > 32.0*TEST_SQRT_TOL0 {
		sa++
	}
	if re.e10 != 0 {
		sa++
	}

	test.Test(t, sa, "Exp_err_e10_e(1.0, TEST_SQRT_TOL0)")

	sa = 0
	stat = Exp_err_e10_e(2000.0, 1.0e-10, re)
	if stat != nil {
		sa += stat.(err.GSLError).Status()
	}

	if test_sf_frac_diff(re.val, 3.88118019428363725) > TEST_TOL3 {
		sa++
	}
	if re.err > 1.0e-07 {
		sa++
	}
	if re.e10 != 868 {
		sa++
	}
	test.Test(t, sa, "Exp_err_e10_e(2000.0, 1.0e-10)")
}

func TestExpmult(t *testing.T) {
	r := new(Result)
	x := 0.8 * gsl.LnMaxFloat64
	cases := []struct {
		x, y, expected float64
		tol            float64
	}{
		{-10.0, 1.0e-06, 1.0e-06 * math.Exp(-10.0), TEST_TOL0},
		{-10.0, 2.0, 2.0 * math.Exp(-10.0), TEST_TOL0},
		{-10.0, -2.0, -2.0 * math.Exp(-10.0), TEST_TOL0},
		{10.0, 1.0e-06, 1.0e-06 * math.Exp(10.0), TEST_TOL0},
		{10.0, -2.0, -2.0 * math.Exp(10.0), TEST_TOL0},
		{x, 1.00001, 1.00001 * math.Exp(x), TEST_TOL3},
		{x, 1.000001, 1.000001 * math.Exp(x), TEST_TOL3},
		{x, 1.000000001, 1.000000001 * math.Exp(x), TEST_TOL3},
		{x, 100.0, 100.0 * math.Exp(x), TEST_TOL3},
		{x, 1.0e+20, 1.0e+20 * math.Exp(x), TEST_TOL3},
		{x, math.Exp(-x) * math.Exp(math.Ln2), 2.0, TEST_TOL4},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exp_mult_e(c.x, c.y, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Exp_mult_e(%v, %v)", c.x, c.y))
		})
	}
}

func TestExpmultErr(t *testing.T) {
	r := new(Result)
	x := 0.8 * gsl.LnMaxFloat64
	cases := []struct {
		x, dx, y, dy, expected float64
		tol                    float64
	}{
		{-10.0, TEST_SQRT_TOL0, 2.0, TEST_SQRT_TOL0, 2.0 * math.Exp(-10.0), TEST_SQRT_TOL0},
		{x, TEST_SQRT_TOL0 * x, math.Exp(-x) * math.Exp(math.Ln2), TEST_SQRT_TOL0 * math.Exp(-x) * math.Exp(math.Ln2), 2.0, TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exp_mult_err_e(c.x, c.dx, c.y, c.dy, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Exp_mult_err_e(%v, %v, %v, %v)", c.x, c.dx, c.y, c.dy))
		})
	}
}

func TestExpMultE10(t *testing.T) {
	re := new(Result_e10)
	sa := 0
	stat := Exp_mult_e10_e(1.0, 1.0, re)
	if stat != nil {
		sa += stat.(err.GSLError).Status()
	}
	if test_sf_frac_diff(re.val, math.E) > TEST_TOL0 {
		sa++
	}
	if re.err > TEST_TOL2 {
		sa++
	}
	if re.e10 != 0 {
		sa++
	}

	test.Test(t, sa, "Exp_mult_e10_e(1.0, 1.0)")
	run_test_sf_e10(t, stat, re, math.E, 0, TEST_TOL0, err.SUCCESS, "Exp_mult_e10_e(1.0, 1.0)")

	sa = 0
	stat = Exp_mult_e10_e(1000.0, 1.0e+200, re)
	if stat != nil {
		sa += stat.(err.GSLError).Status()
	}

	if test_sf_frac_diff(re.val, 1.970071114017046993888879352) > TEST_TOL3 {
		sa++
	}
	if re.err > 1.0e-11 {
		sa++
	}
	if re.e10 != 634 {
		sa++
	}

	test.Test(t, sa, "Exp_mult_e10_e(1000.0, 1.0e+200)")
	run_test_sf_e10(t, stat, re, 1.970071114017046993888879352, 634, TEST_TOL3, err.SUCCESS, "Exp_mult_e10_e(1000.0, 1.0e+200)")

}

func TestExpMultErrE10(t *testing.T) {
	re := new(Result_e10)

	sa := 0
	stat := Exp_mult_err_e10_e(1.0, TEST_TOL0, 1.0, TEST_TOL0, re)
	if stat != nil {
		sa += stat.(err.GSLError).Status()
	}

	if test_sf_frac_diff(re.val, math.E) > TEST_TOL0 {
		sa++
	}
	if re.err > TEST_TOL2 {
		sa++
	}
	if re.e10 != 0 {
		sa++
	}

	test.Test(t, sa, "Exp_mult_err_e10_e(1.0, TEST_TOL0, 1.0, TEST_TOL0)")
	run_test_sf_e10(t, stat, re, math.E, 0, TEST_TOL0, err.SUCCESS, "Exp_mult_err_e10_e(1.0, TEST_TOL0, 1.0, TEST_TOL0)")

	sa = 0
	stat = Exp_mult_err_e10_e(1000.0, 1.0e-12, 1.0e+200, 1.0e+190, re)
	if stat != nil {
		sa += stat.(err.GSLError).Status()
	}

	if test_sf_frac_diff(re.val, 1.9700711140165661) > TEST_TOL3 {
		sa++
	}
	if re.err > 1.0e-09 {
		sa++
	}
	if re.e10 != 634 {
		sa++
	}

	test.Test(t, sa, "Exp_mult_err_e10_e(1000.0, 1.0e-12, 1.0e+200, 1.0e+190)")
	run_test_sf_e10(t, stat, re, 1.9700711140165661, 634, TEST_TOL3, err.SUCCESS, "Exp_mult_err_e10_e(1000.0, 1.0e-12, 1.0e+200, 1.0e+190)")

}

func TestExpMultE10_2(t *testing.T) {
	re := new(Result_e10)
	cases := []struct {
		x, y, expected float64
		e10            int
		tol            float64
	}{
		{10000.0, 1.0, 8.806818225662921587261496007, 4342, TEST_TOL5},
		{100.0, 1.0, 2.688117141816135448412625551e43, 0, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exp_mult_e10_e(c.x, c.y, re)
			run_test_sf_e10(t, stat, re, c.expected, c.e10, c.tol, err.SUCCESS, fmt.Sprintf("Exp_mult_e10_e(%v, %v)", c.x, c.y))
		})
	}
}

func TestExpE10_2(t *testing.T) {
	re := new(Result_e10)
	cases := []struct {
		x, expected float64
		e10         int
		tol         float64
	}{
		{100.0, 2.688117141816135448412625551e43, 0, TEST_TOL2},
		{1000.0, 1.970071114017046993888879352, 434, TEST_TOL3},
		{-100.0, 3.720075976020835962959695803e-44, 0, TEST_TOL2},
		{-1000.0, 5.075958897549456765291809479, -435, TEST_TOL3},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exp_e10_e(c.x, re)
			run_test_sf_e10(t, stat, re, c.expected, c.e10, c.tol, err.SUCCESS, fmt.Sprintf("Exp_e10_e(%v)", c.x))
		})
	}
}

func TestExpm1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10.0, math.Exp(-10.0) - 1.0, TEST_TOL0},
		{-0.001, -0.00099950016662500845, TEST_TOL0},
		{-1.0e-8, -1.0e-08 + 0.5e-16, TEST_TOL0},
		{1.0e-8, 1.0e-08 + 0.5e-16, TEST_TOL0},
		{0.001, 0.0010005001667083417, TEST_TOL0},
		{10.0, math.Exp(10.0) - 1.0, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expm1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expm1_e(%v)", c.x))
		})
	}
}

func TestExprel(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10.0, 0.0999954600070237515, TEST_TOL0},
		{-0.001, 0.9995001666250084, TEST_TOL0},
		{-1.0e-8, 1.0 - 0.5e-08, TEST_TOL0},
		{1.0e-8, 1.0 + 0.5e-08, TEST_TOL0},
		{0.001, 1.0005001667083417, TEST_TOL0},
		{10.0, 2202.5465794806716517, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exprel_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Exprel_e(%v)", c.x))
		})
	}
}

func TestExprel2(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10.0, 0.18000090799859524970, TEST_TOL0},
		{-0.001, 0.9996667499833361107, TEST_TOL0},
		{-1.0e-8, 0.9999999966666666750, TEST_TOL0},
		{1.0e-8, 1.0000000033333333417, TEST_TOL0},
		{0.001, 1.0003334166833361115, TEST_TOL0},
		{10.0, 440.3093158961343303, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exprel_2_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Exprel_2_e(%v)", c.x))
		})
	}
}

func TestExprelN(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{3, -1000.0, 0.00299400600000000000, TEST_TOL0},
		{3, -100.0, 0.02940600000000000000, TEST_TOL0},
		{3, -10.0, 0.24599972760042142509, TEST_TOL0},
		{3, -3.0, 0.5444917625849191238, TEST_TOL0},
		{3, -0.001, 0.9997500499916678570, TEST_TOL0},
		{3, -1.0e-8, 0.9999999975000000050, TEST_TOL0},
		{3, 1.0e-8, 1.0000000025000000050, TEST_TOL0},
		{3, 0.001, 1.0002500500083345240, TEST_TOL0},
		{3, 3.0, 2.5745637607083706091, TEST_TOL0},
		{3, 3.1, 2.6772417068460206247, TEST_TOL0},
		{3, 10.0, 131.79279476884029910, TEST_TOL1},
		{3, 100.0, 1.6128702850896812690e+38, TEST_TOL2},
		{50, -1000.0, 0.04766231609253975959, TEST_TOL0},
		{50, -100.0, 0.3348247572345889317, TEST_TOL0},
		{50, -10.0, 0.8356287051853286482, TEST_TOL0},
		{50, -3.0, 0.9443881609152163615, TEST_TOL0},
		{50, -1.0, 0.980762245565660617, TEST_TOL0},
		{50, -1.0e-8, 1.0 - 1.0e-8/51.0, TEST_TOL0},
		{50, 1.0e-8, 1.0 + 1.0e-8/51.0, TEST_TOL0},
		{50, 1.0, 1.01999216583666790, TEST_TOL0},
		{50, 3.0, 1.0624205757460368307, TEST_TOL0},
		{50, 48.0, 7.499573876877194416, TEST_TOL0},
		{50, 50.1, 9.311803306230992272, TEST_TOL4},
		{50, 100.0, 8.175664432485807634e+07, TEST_TOL4},
		{50, 500.0, 4.806352370663185330e+146, TEST_TOL3},
		{500, -1000.0, 0.3334815803127619256, TEST_TOL0},
		{500, -100.0, 0.8335646217536183909, TEST_TOL0},
		{500, -10.0, 0.9804297803131823066, TEST_TOL0},
		{500, -3.0, 0.9940475488850672997, TEST_TOL0},
		{500, -1.0, 0.9980079602383488808, TEST_TOL0},
		{500, -1.0e-8, 1.0 - 1.0e-8/501.0, TEST_TOL0},
		{500, 1.0e-8, 1.0 + 1.0e-8/501.0, TEST_TOL0},
		{500, 1.0, 1.0019999920160634252, TEST_TOL0},
		{500, 3.0, 1.0060240236632444934, TEST_TOL0},
		{500, 48.0, 1.1059355517981272174, TEST_TOL0},
		{500, 100.0, 1.2492221464878287204, TEST_TOL1},
		{500, 500.0, 28.363019877927630858, TEST_TOL2},
		{500, 1000.0, 2.4037563160335300322e+68, TEST_TOL4},
		{500, 1600.0, 7.899293535320607403e+226, TEST_TOL4},
		{1263131.0, 1261282.3637, 545.0113107238425900305428360, TEST_TOL4},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exprel_n_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Exprel_n_e(%v, %v)", c.n, c.x))
		})
	}
}

func TestExprelNCF(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n, x, expected float64
		tol            float64
	}{
		{6.315655e+05, 6.302583168053568806e+05, 385.425369029433473098652465720, TEST_TOL4},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Exprel_n_CF_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Exprel_n_CF_e(%v, %v)", c.n, c.x))
		})
	}
}
