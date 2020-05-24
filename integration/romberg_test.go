package integration

import (
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/test"
	"testing"
)

func TestRomberg(t *testing.T) {
	initTesting()
	/* test Romberg integration */
	status := 0
	var result, exp_result float64
	var neval int
	exp_ier := 0

	w, _ := NewRombergWorkspace(20)

	f := f_sin{}
	exp_result = 1.0

	errqg := Romberg(&f, 0.0, gsl.PiOver2, 0.0, 1e-10, &result, &neval, w)
	if errqg != nil {
		status = errqg.Status()
	}
	test.TestInt(t, status, exp_ier, "romberg(f_sin) status")
	test.TestRel(t, result, exp_result, 1e-15, "romberg(f_sin) result")

	errqg = Romberg(&f, gsl.PiOver2, 0.0, 0.0, 1e-10, &result, &neval, w)
	if errqg != nil {
		status = errqg.Status()
	}
	test.TestInt(t, status, exp_ier, "romberg(f_sin) reverse status")
	test.TestRel(t, result, -exp_result, 1e-15, "romberg(f_sin) reverse result")

	f2 := cqf11{}
	exp_result = 5.0

	errqg = Romberg(&f2, -5.0, 5.0, 0.0, 1e-10, &result, &neval, w)
	if errqg != nil {
		status = errqg.Status()
	}
	test.TestInt(t, status, exp_ier, "romberg(cqf11) status")
	test.TestRel(t, result, exp_result, 1e-15, "romberg(cqf11) result")

	errqg = Romberg(&f2, 5.0, -5.0, 0.0, 1e-10, &result, &neval, w)
	if errqg != nil {
		status = errqg.Status()
	}
	test.TestInt(t, status, exp_ier, "romberg(cqf11) reverse status")
	test.TestRel(t, result, -exp_result, 1e-15, "romberg(cqf11) reverse result")
}
