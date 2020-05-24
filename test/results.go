package test

import (
	"fmt"
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/sys"
	"math"
	"os"
	"strconv"
	"testing"
)

var (
	tests   = 0
	passed  = 0
	failed  = 0
	verbose uint64
)

func initialise() {
	p := os.Getenv("GSL_TEST_VERBOSE")
	var err error

	/* 0 = show failures only (we always want to see these) */
	/* 1 = show passes and failures */

	if p == "" { /* environment variable is not set */
		return
	}

	verbose, err = strconv.ParseUint(p, 10, 64)

	if err != nil {
		verbose = 0
	}

	return
}

func update(s int) {
	tests++
	if s == 0 {
		passed++
	} else {
		failed++
	}
}

func Test(t *testing.T, status int, test_description string) {
	if tests == 0 {
		initialise()
	}

	update(status)

	if status != 0 || verbose != 0 {
		if status == 0 {
			fmt.Printf("PASS: %s", test_description)
		} else {
			s := fmt.Sprintf("FAIL: %s", test_description)
			if verbose == 0 {
				s += fmt.Sprintf(" [%v]", tests)
			}

			t.Errorf("%s", s)
		}

		fmt.Printf("\n")
	}
}

func TestRel(t *testing.T, result, expected, relative_error float64, test_description string) {
	var status int

	if tests == 0 {
		initialise()
	}

	/* Check for NaN vs inf vs number */
	if sys.IsNaN(result) == 1 || sys.IsNaN(expected) == 1 {
		if sys.IsNaN(result) != sys.IsNaN(expected) {
			status = 1
		}
	} else if sys.IsInf(result) == 1 || sys.IsInf(expected) == 1 {
		if sys.IsInf(result) != sys.IsInf(expected) {
			status = 1
		}
	} else if (expected > 0 && expected < gsl.MinFloat64) || (expected < 0 && expected > -(gsl.MinFloat64)) {
		status = -1
	} else if expected != 0 {
		if math.Abs(result-expected)/math.Abs(expected) > relative_error {
			status = 1
		}
	} else {
		if math.Abs(result) > relative_error {
			status = 1
		}
	}

	update(status)
	if status != 0 || verbose != 0 {
		if status == 0 {
			fmt.Printf("PASS: %s (%g observed vs %g expected)", test_description, result, expected)
		} else {
			s := fmt.Sprintf("FAIL: %s (%.18g observed vs %.18g expected)", test_description, result, expected)
			if status == -1 {
				s += " [test uses subnormal value]"
			}

			if verbose == 0 {
				s += fmt.Sprintf(" [%v]", tests)
			}

			t.Errorf("%s", s)
		}

		fmt.Printf("\n")
	}
}

func TestAbs(t *testing.T, result, expected, absolute_error float64, test_description string) {
	var status int

	if tests == 0 {
		initialise()
	}

	/* Check for NaN vs inf vs number */
	if sys.IsNaN(result) == 1 || sys.IsNaN(expected) == 1 {
		if sys.IsNaN(result) != sys.IsNaN(expected) {
			status = 1
		}
	} else if sys.IsInf(result) == 1 || sys.IsInf(expected) == 1 {
		if sys.IsInf(result) != sys.IsInf(expected) {
			status = 1
		}
	} else if (expected > 0 && expected < gsl.MinFloat64) || (expected < 0 && expected > -(gsl.MinFloat64)) {
		status = -1
	} else {
		if math.Abs(result-expected) > absolute_error {
			status = 1
		}
	}

	update(status)

	if status != 0 || verbose != 0 {
		if status == 0 {
			fmt.Printf("PASS: %s (%g observed vs %g expected)", test_description, result, expected)
		} else {
			s := fmt.Sprintf("FAIL: %s (%.18g observed vs %.18g expected)", test_description, result, expected)
			if status == -1 {
				s += " [test uses subnormal value]"
			}

			if verbose == 0 {
				s += fmt.Sprintf(" [%v]", tests)
			}

			t.Errorf("%s", s)
		}

		fmt.Printf("\n")
	}
}

func TestInt(t *testing.T, result, expected int, test_description string) {
	var status int
	if result != expected {
		status = 1
	}

	if tests == 0 {
		initialise()
	}

	update(status)

	if status != 0 || verbose != 0 {
		if status == 0 {
			fmt.Printf("PASS: %s (%d observed vs %d expected)", test_description, result, expected)
		} else {
			s := fmt.Sprintf("FAIL: %s (%d observed vs %d expected)", test_description, result, expected)
			if verbose == 0 {
				s += fmt.Sprintf(" [%v]", tests)
			}

			t.Errorf("%s", s)
		}

		fmt.Printf("\n")
	}
}
