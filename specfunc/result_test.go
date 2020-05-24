package specfunc

import (
	"testing"
)

func TestResult(t *testing.T) {
	var s int
	var re Result_e10
	var r Result
	re.val = -1.0
	re.err = 0.5
	re.e10 = 0
	Result_smash_e(&re, &r)
	if test_sf_frac_diff(r.val, -1.0) > TEST_TOL0 {
		s++
	}
	if test_sf_frac_diff(r.err, 0.5) > TEST_TOL0 {
		s++
	}

	re.val = -1.0
	re.err = 0.5
	re.e10 = 10
	Result_smash_e(&re, &r)
	if test_sf_frac_diff(r.val, -1.0e+10) > TEST_TOL1 {
		s++
	}
	if test_sf_frac_diff(r.err, 0.5e+10) > TEST_TOL1 {
		s++
	}

	if s != 0 {
		t.Errorf("Failed")
	}
}
