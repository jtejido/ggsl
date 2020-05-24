package integration

import (
	"github.com/lucky-se7en/ggsl/err"
	"github.com/lucky-se7en/ggsl/test"
	"testing"
)

func TestQng(t *testing.T) {
	initTesting()
	{
		var result, abserr float64
		var status, neval int
		exp_result := 7.716049379303083211e-02
		exp_abserr := 9.424302199601294244e-08
		exp_neval := 21
		exp_ier := 0

		alpha := 2.6
		f := &f1{alpha}

		errqg := Qng(f, 0.0, 1.0, 1e-1, 0.0, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, exp_result, 1e-15, "qng(f1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qng(f1) smooth abserr")
		test.TestInt(t, neval, exp_neval, "qng(f1) smooth neval")
		test.TestInt(t, status, exp_ier, "qng(f1) smooth status")

		errqg = Qng(f, 1.0, 0.0, 1e-1, 0.0, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, -exp_result, 1e-15, "qng(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qng(f1) reverse abserr")
		test.TestInt(t, neval, exp_neval, "qng(f1) reverse neval")
		test.TestInt(t, status, exp_ier, "qng(f1) reverse status")
	}

	{
		var result, abserr float64
		var status, neval int

		exp_result := 7.716049382706505200e-02
		exp_abserr := 2.666893044866214501e-12
		exp_neval := 43
		exp_ier := 0

		alpha := 2.6
		f := &f1{alpha}

		errqg := Qng(f, 0.0, 1.0, 0.0, 1e-9, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, exp_result, 1e-15, "qng(f1) smooth 43pt result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qng(f1) smooth 43pt abserr")
		test.TestInt(t, neval, exp_neval, "qng(f1) smooth 43pt neval")
		test.TestInt(t, status, exp_ier, "qng(f1) smooth 43pt status")

		errqg = Qng(f, 1.0, 0.0, 0.0, 1e-9, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, -exp_result, 1e-15, "qng(f1) reverse 43pt result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qng(f1) reverse 43pt abserr")
		test.TestInt(t, neval, exp_neval, "qng(f1) reverse 43pt neval")
		test.TestInt(t, status, exp_ier, "qng(f1) reverse 43pt status")
	}

	{
		var result, abserr float64
		var status, neval int
		exp_result := -7.238969575482961938e-01
		exp_abserr := 1.277676889520056369e-14
		exp_neval := 43
		exp_ier := 0

		alpha := 1.3
		f := &f3{alpha}

		errqg := Qng(f, 0.3, 2.71, 0.0, 1e-12, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, exp_result, 1e-15, "qnq(f3) oscill result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qng(f3) oscill abserr")
		test.TestInt(t, neval, exp_neval, "qng(f3) oscill neval")
		test.TestInt(t, status, exp_ier, "qng(f3) oscill status")

		errqg = Qng(f, 2.71, 0.3, 0.0, 1e-12, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, -exp_result, 1e-15, "qnq(f3) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qng(f3) reverse abserr")
		test.TestInt(t, neval, exp_neval, "qng(f3) reverse neval")
		test.TestInt(t, status, exp_ier, "qng(f3) reverse status")
	}

	{
		var result, abserr float64
		var status, neval int

		exp_result := 7.716049382716029525e-02
		exp_abserr := 8.566535680046930668e-16
		exp_neval := 87
		exp_ier := 0

		alpha := 2.6
		f := &f1{alpha}

		errqg := Qng(f, 0.0, 1.0, 0.0, 1e-13, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, exp_result, 1e-15, "qng(f1) 87pt smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qng(f1) 87pt smooth abserr")
		test.TestInt(t, neval, exp_neval, "qng(f1) 87pt smooth neval")
		test.TestInt(t, status, exp_ier, "qng(f1) 87pt smooth status")

		errqg = Qng(f, 1.0, 0.0, 0.0, 1e-13, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, -exp_result, 1e-15, "qng(f1) 87pt reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qng(f1) 87pt reverse abserr")
		test.TestInt(t, neval, exp_neval, "qng(f1) 87pt reverse neval")
		test.TestInt(t, status, exp_ier, "qng(f1) 87pt reverse status")
	}

	{
		var result, abserr float64
		var status, neval int

		exp_result := 3.222948711817264211e+01
		exp_abserr := 2.782360287710622870e+01
		exp_neval := 87
		exp_ier := err.ETOL

		alpha := -0.9
		f := f1{alpha}

		errqg := Qng(f, 0.0, 1.0, 0.0, 1e-3, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, exp_result, 1e-15, "qng(f1) sing beyond 87pt result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qng(f1) sing beyond 87pt abserr")
		test.TestInt(t, neval, exp_neval, "qng(f1) sing beyond 87pt neval")
		test.TestInt(t, status, exp_ier, "qng(f1) sing beyond 87pt status")

		errqg = Qng(f, 1.0, 0.0, 0.0, 1e-3, &result, &abserr, &neval)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, -exp_result, 1e-15, "qng(f1) reverse beyond 87pt result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qng(f1) rev beyond 87pt abserr")
		test.TestInt(t, neval, exp_neval, "qng(f1) rev beyond 87pt neval")
		test.TestInt(t, status, exp_ier, "qng(f1) rev beyond 87pt status")
	}
}
