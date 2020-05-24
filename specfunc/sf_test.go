package specfunc

import (
	"bytes"
	"fmt"
	gsl "github.com/jtejido/ggsl"
	. "github.com/jtejido/ggsl/err"
	"github.com/jtejido/ggsl/sys"
	"github.com/jtejido/ggsl/test"
	"math"
	"os"
	"testing"
)

func initTesting() {
	// original test just disables handler but if you'd like to know where it catches the flu.
	// SetErrorHandler(myErrorHandler)
	SetErrorHandlerOff()
	os.Setenv("GSL_TEST_VERBOSE", "1")
}

func myErrorHandler(reason, file string, line, gsl_errno int) {
	fmt.Printf("(error expected [%s:%d: %s (%d)])\n", file, line, reason, gsl_errno)
}

const (
	TEST_TOL0      = 2.0 * gsl.Float64Eps
	TEST_TOL1      = 16.0 * gsl.Float64Eps
	TEST_TOL2      = 256.0 * gsl.Float64Eps
	TEST_TOL3      = 2048.0 * gsl.Float64Eps
	TEST_TOL4      = 16384.0 * gsl.Float64Eps
	TEST_TOL5      = 131072.0 * gsl.Float64Eps
	TEST_TOL6      = 1048576.0 * gsl.Float64Eps
	TEST_SNGL      = 1.0e-06
	TEST_SQRT_TOL0 = 2.0 * gsl.SqrtFloat64Eps
)

const (
	test_sf_incons = 1
	test_sf_errneg = 2
	test_sf_tolbad = 4
	test_sf_retbad = 8
	test_sf_errbad = 16
	test_sf_errbig = 32
	test_sf_expbad = 64
	test_sigma     = 1.5
	test_factor    = 100.0
)

const (
	incons = "  value/expected not consistent within reported error\n"
	errneg = "  reported error negative\n"
	errbad = "  reported error is bad\n"
	errbig = "  reported error is much too big\n"
	tolbad = "  value not within tolerance of expected value\n"
	expbad = "  exponent is incorrect\n"
)

func test_sf_frac_diff(x1, x2 float64) float64 {
	if x1 == 0.0 && x2 == 0.0 {
		return 0.0
	} else if x1 == 0.0 {
		return math.Abs(x2)
	} else if x1 <= math.MaxFloat64 && x2 <= math.MaxFloat64 && (x1+x2 != 0.0) {
		return math.Abs((x1 - x2) / (x1 + x2))
	}

	return 1.0
}

func isAlmostEqual(a, b, tol float64) bool {
	return test_sf_frac_diff(a, b) <= tol
}

func run_test_sf_return(t *testing.T, sstat GSLError, expect_return int, desc string) {
	var status int
	if sstat != nil {
		status = sstat.Status()
	}

	initTesting()
	test_sf_return(t, status, expect_return, desc)
}

func run_test_sf(t *testing.T, sstat GSLError, r *Result, val_in, tol float64, expect_return int, desc string) {
	var status int
	if sstat != nil {
		status = sstat.Status()
	}
	initTesting()
	test_sf(t, r, val_in, tol, status, expect_return, desc)
}

func run_test_sf_2(t *testing.T, sstat GSLError, r1 *Result, val1, tol1 float64, r2 *Result, val2, tol2 float64, expect_return int, desc string) {
	var status int
	if sstat != nil {
		status = sstat.Status()
	}
	initTesting()
	test_sf_2(t, r1, val1, tol1, r2, val2, tol2, status, expect_return, desc)
}

func run_test_sf_e10(t *testing.T, sstat GSLError, re *Result_e10, val_in float64, e10_in int, tol float64, expect_return int, desc string) {
	var status int
	if sstat != nil {
		status = sstat.Status()
	}
	initTesting()
	test_sf_e10(t, re, val_in, e10_in, tol, status, expect_return, desc)
}

func run_test_theta(t *testing.T, val, val_in, tol float64, desc string) {
	initTesting()
	test_sf_val(t, val, val_in, tol, desc)
}

func run_test_val(t *testing.T, val, val_in, tol float64, desc string) {
	initTesting()
	test_sf_val(t, val, val_in, tol, desc)
}

func test_sf_val(t *testing.T, val, val_in, tol float64, desc string) int {
	var message_buff bytes.Buffer
	local_s := 0

	message_buff.WriteByte(byte('\000'))

	local_s |= test_sf_check_val(&message_buff, val, val_in, tol)

	test.Test(t, local_s, desc)

	if local_s != 0 {
		t.Errorf("\n%s\n  %22.18e\n", message_buff.String(), val)
	}

	return local_s
}

func test_sf_return(t *testing.T, status, expect_return int, desc string) int {
	var message_buff bytes.Buffer
	local_s := 0

	message_buff.WriteByte(byte('\000'))

	local_s |= test_sf_check_return(&message_buff, status, expect_return)

	test.Test(t, local_s, desc)
	if local_s != 0 {
		t.Errorf("\n%s\n", message_buff.String())
	}
	return local_s
}

func test_sf(t *testing.T, r *Result, val_in, tol float64, status, expect_return int, desc string) int {
	var message_buff bytes.Buffer
	local_s := 0

	message_buff.WriteByte(byte('\000'))
	local_s |= test_sf_check_result(&message_buff, r, val_in, tol)
	local_s |= test_sf_check_return(&message_buff, status, expect_return)

	test.Test(t, local_s, desc)

	if local_s != 0 {
		t.Errorf("\n%s  result: %22.18e  error: %22.18e\n", message_buff.String(), r.val, r.err)
	}

	return local_s
}

func test_sf_e10(t *testing.T, re *Result_e10, val_in float64, e10_in int, tol float64, status, expect_return int, desc string) int {
	var message_buff bytes.Buffer
	local_s := 0
	r := new(Result)
	r.val = re.val
	r.err = re.err

	message_buff.WriteByte(byte('\000'))
	local_s |= test_sf_check_result(&message_buff, r, val_in, tol)
	local_s |= test_sf_check_e10(&message_buff, re.e10, e10_in)
	local_s |= test_sf_check_return(&message_buff, status, expect_return)

	test.Test(t, local_s, desc)
	if local_s != 0 {
		t.Errorf("\n%s  result: %22.18e  error: %22.18e  e10: 10^%d\n", message_buff.String(), re.val, re.err, re.e10)
	}
	return local_s
}

func test_sf_rlx(t *testing.T, r *Result, val_in, tol float64, status, expect_return int, desc string) int {
	var message_buff bytes.Buffer
	local_s := 0

	message_buff.WriteByte(byte('\000'))

	local_s |= test_sf_check_result_relax(&message_buff, r, val_in, tol)
	local_s |= test_sf_check_return(&message_buff, status, expect_return)

	test.Test(t, local_s, desc)
	if local_s != 0 {
		t.Errorf("\n%s  result: %22.18e  error: %22.18e\n", message_buff.String(), r.val, r.err)
	}
	return local_s
}

func test_sf_2(t *testing.T, r1 *Result, val1, tol1 float64, r2 *Result, val2, tol2 float64, status, expect_return int, desc string) int {
	var message_buff bytes.Buffer
	local_s := 0

	message_buff.WriteByte(byte('\000'))

	local_s |= test_sf_check_result(&message_buff, r1, val1, tol1)
	local_s |= test_sf_check_result(&message_buff, r2, val2, tol2)
	local_s |= test_sf_check_return(&message_buff, status, expect_return)

	test.Test(t, local_s, desc)
	if local_s != 0 {
		t.Errorf("\n%s  result1: %22.18e  error1: %22.18e\n   result2: %22.18e  error2: %22.18e\n", message_buff.String(), r1.val, r1.err, r2.val, r2.err)
	}
	return local_s
}

func test_sf_check_result(message_buff *bytes.Buffer, r *Result, val, tol float64) int {
	s := 0
	f := 0.0
	d := 0.0

	if sys.IsNaN(r.val) == 1 || sys.IsNaN(val) == 1 {
		if sys.IsNaN(r.val) != sys.IsNaN(val) {
			s = test_sf_incons
		}
	} else if sys.IsInf(r.val) == 1 || sys.IsInf(val) == 1 {
		if sys.IsInf(r.val) != sys.IsInf(val) {
			s = test_sf_incons
		}
	} else {
		f = test_sf_frac_diff(val, r.val)
		d = math.Abs(val - r.val)

		if d > 2.0*test_sigma*r.err {
			s |= test_sf_incons
		}
		if r.err < 0.0 {
			s |= test_sf_errneg
		}
		if sys.IsInf(r.err) == 1 {
			s |= test_sf_errbad
		}
		if f > test_factor*tol {
			s |= test_sf_tolbad
		}
	}

	if s != 0 {
		message_buff.Write([]byte(fmt.Sprintf("  expected: %20.16e\n", val)))
		message_buff.Write([]byte(fmt.Sprintf("  obtained: %20.16e +/- %.16e (rel=%g)\n", r.val, r.err, r.err/(math.Abs(r.val)+r.err))))
		message_buff.Write([]byte(fmt.Sprintf("  fracdiff: %20.16e\n", f)))
		message_buff.Write([]byte(fmt.Sprintf("  tolerance: %20.16e\n", tol)))
	}

	if s&test_sf_incons == 1 {
		message_buff.Write([]byte(incons))
	}
	if s&test_sf_errneg == 1 {
		message_buff.Write([]byte(errneg))
	}
	if s&test_sf_errbad == 1 {
		message_buff.Write([]byte(errbad))
	}
	if s&test_sf_errbig == 1 {
		message_buff.Write([]byte(errbig))
	}
	if s&test_sf_tolbad == 1 {
		message_buff.Write([]byte(tolbad))
	}

	return s
}

func test_sf_check_e10(message_buff *bytes.Buffer, e10, e10_in int) int {
	s := 0
	if e10 != e10_in {
		s = test_sf_expbad
	}

	if s != 0 {
		message_buff.Write([]byte(fmt.Sprintf("  expected exponent: 10^%d\n", e10_in)))
		message_buff.Write([]byte(fmt.Sprintf("  obtained exponent: 10^%d", e10)))
	}

	if s&test_sf_expbad == 1 {
		message_buff.Write([]byte(expbad))
	}

	return s
}

func test_sf_check_val(message_buff *bytes.Buffer, rval, val, tol float64) int {
	s := 0
	f := test_sf_frac_diff(val, rval)

	if f > test_factor*tol {
		s |= test_sf_tolbad
	}

	if s != 0 {
		message_buff.Write([]byte(fmt.Sprintf("  expected: %20.16e\n", val)))
		message_buff.Write([]byte(fmt.Sprintf("  obtained: %20.16e\n", rval)))
		message_buff.Write([]byte(fmt.Sprintf("  fracdiff: %20.16e\n", f)))
	}

	if s&test_sf_tolbad == 1 {
		message_buff.Write([]byte(fmt.Sprintf(tolbad)))
	}

	return s
}

func test_sf_check_result_relax(message_buff *bytes.Buffer, r *Result, val, tol float64) int {
	s := 0
	f := test_sf_frac_diff(val, r.val)

	if f > gsl.Max(TEST_SNGL, test_factor*tol) {
		s |= test_sf_incons
	}
	if r.err < 0.0 {
		s |= test_sf_errneg
	}
	if sys.IsInf(r.err) == 1 {
		s |= test_sf_errbad
	}
	if f > test_factor*tol {
		s |= test_sf_tolbad
	}

	if s != 0 {
		message_buff.Write([]byte(fmt.Sprintf("  expected: %20.16e\n", val)))
		message_buff.Write([]byte(fmt.Sprintf("  obtained: %20.16e +/- %.16e  (rel=%g)\n", r.val, r.err, r.err/(math.Abs(r.val)+r.err))))
		message_buff.Write([]byte(fmt.Sprintf("  fracdiff: %20.16e\n", f)))
	}

	if s&test_sf_incons == 1 {
		message_buff.Write([]byte(fmt.Sprintf("  value/expected not consistent MAX(tol,SINGLE_PREC)\n")))
	}
	if s&test_sf_errneg == 1 {
		message_buff.Write([]byte(fmt.Sprintf(errneg)))
	}
	if s&test_sf_errbad == 1 {
		message_buff.Write([]byte(fmt.Sprintf(errbad)))
	}
	if s&test_sf_tolbad == 1 {
		message_buff.Write([]byte(fmt.Sprintf(tolbad)))
	}

	return s
}

func test_sf_check_return(message_buff *bytes.Buffer, val_return, expected_return int) int {
	if val_return != expected_return {
		message_buff.Write([]byte(fmt.Sprintf("  unexpected return code: %d\n", val_return)))
		return test_sf_retbad
	} else {
		return 0
	}
}
