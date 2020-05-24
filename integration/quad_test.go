package integration

import (
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"github.com/jtejido/ggsl/test"
	"math"
	"testing"
)

func TestQuad(t *testing.T) {
	/*	QAG  */
	{
		initTesting()
		w, _ := NewWorkspace(1000)
		exp_result := 7.716049382715854665e-02
		exp_abserr := 6.679384885865053037e-12
		exp_last := 6
		exp_neval := 165
		var result, abserr float64
		var status, exp_ier int

		a := []float64{0, 0.5, 0.25, 0.125, 0.0625, 0.03125}
		b := []float64{0.03125, 1, 0.5, 0.25, 0.125, 0.0625}
		r := []float64{3.966769831709074375e-06, 5.491842501998222409e-02,
			1.909827770934243926e-02, 2.776531175604360531e-03,
			3.280661030752063693e-04, 3.522704932261797744e-05}
		e := []float64{6.678528276336181873e-12, 6.097169993333454062e-16,
			2.120334764359736934e-16, 3.082568839745514608e-17,
			3.642265412331439511e-18, 3.910988124757650942e-19}
		order := []int{1, 2, 3, 4, 5, 6}

		alpha := 2.6
		f := f1{alpha: alpha}
		fc := &counter{f, 0}
		errqg := Qag(fc, 0.0, 1.0, 0.0, 1e-10, w.limit, w, &result, &abserr, Qk15)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-15, "qag(f1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qag(f1) smooth abserr")
		test.TestInt(t, fc.neval, exp_neval, "qag(f1) smooth neval")
		test.TestInt(t, w.size, exp_last, "qag(f1) smooth last")
		test.TestInt(t, status, exp_ier, "qag(f1) smooth status")

		for i := 0; i < 6; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qag(f1) smooth alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qag(f1) smooth blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-15, "qag(f1) smooth rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-6, "qag(f1) smooth elist")
			test.TestInt(t, int(w.order[i]), order[i]-1, "qag(f1) smooth order")
		}

		fc.neval = 0
		status = 0
		errqg = Qag(fc, 1.0, 0.0, 0.0, 1e-10, w.limit, w, &result, &abserr, Qk15)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, -exp_result, 1e-15, "qag(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qag(f1) reverse abserr")
		test.TestInt(t, fc.neval, exp_neval, "qag(f1) smooth neval")
		test.TestInt(t, int(w.size), exp_last, "qag(f1) reverse last")
		test.TestInt(t, status, exp_ier, "qag(f1) reverse status")
	}

	/* Test the same function using an absolute error bound and the
	   21-point rule */

	{
		w, _ := NewWorkspace(1000)
		exp_result := 7.716049382716050342e-02
		exp_abserr := 2.227969521869139532e-15
		exp_neval := 315
		exp_last := 8
		var result, abserr float64
		var status, exp_ier int
		var i int

		a := []float64{0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125}
		b := []float64{0.0078125, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625}
		r := []float64{3.696942726831556522e-08, 5.491842501998223103e-02,
			1.909827770934243579e-02, 2.776531175604360097e-03,
			3.280661030752062609e-04, 3.522704932261797744e-05,
			3.579060884684503576e-06, 3.507395216921808047e-07}
		e := []float64{1.371316364034059572e-15, 6.097169993333454062e-16,
			2.120334764359736441e-16, 3.082568839745514608e-17,
			3.642265412331439511e-18, 3.910988124757650460e-19,
			3.973555800712018091e-20, 3.893990926286736620e-21}
		order := []int{1, 2, 3, 4, 5, 6, 7, 8}

		alpha := 2.6
		f := f1{alpha: alpha}
		fc := &counter{f, 0}
		errqg := Qag(fc, 0.0, 1.0, 1e-14, 0.0, w.limit, w, &result, &abserr, Qk21)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, exp_result, 1e-15, "qag(f1,21pt) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qag(f1,21pt) smooth abserr")
		test.TestInt(t, fc.neval, exp_neval, "qag(f1,21pt) smooth neval")
		test.TestInt(t, int(w.size), exp_last, "qag(f1,21pt) smooth last")
		test.TestInt(t, status, exp_ier, "qag(f1,21pt) smooth status")

		for i = 0; i < 8; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qag(f1,21pt) smooth alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qag(f1,21pt) smooth blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-15, "qag(f1,21pt) smooth rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-6, "qag(f1,21pt) smooth elist")
			test.TestInt(t, int(w.order[i]), order[i]-1, "qag(f1,21pt) smooth order")
		}

		fc.neval = 0
		errqg = Qag(fc, 1.0, 0.0, 1e-14, 0.0, w.limit, w, &result, &abserr, Qk21)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, -exp_result, 1e-15, "qag(f1,21pt) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qag(f1,21pt) reverse abserr")
		test.TestInt(t, fc.neval, exp_neval, "qag(f1,21pt) reverse neval")
		test.TestInt(t, int(w.size), exp_last, "qag(f1,21pt) reverse last")
		test.TestInt(t, status, exp_ier, "qag(f1,21pt) reverse status")
	}

	/* Adaptive integration of an oscillatory function which terminates because
	   of roundoff error, uses the 31-pt rule */
	{
		w, _ := NewWorkspace(1000)
		exp_result := -7.238969575482959717e-01
		exp_abserr := 1.285805464427459261e-14
		exp_last := 1
		exp_neval := 31
		exp_ier := err.EROUND
		var result, abserr float64
		var status int

		alpha := 1.3
		f := f3{alpha: alpha}
		fc := &counter{f, 0}
		errqg := Qag(fc, 0.3, 2.71, 1e-14, 0.0, w.limit, w, &result, &abserr, Qk31)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, exp_result, 1e-15, "qag(f3,31pt) oscill result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qag(f3,31pt) oscill abserr")
		test.TestInt(t, fc.neval, exp_neval, "qag(f3,31pt) oscill neval")
		test.TestInt(t, int(w.size), exp_last, "qag(f3,31pt) oscill last")
		test.TestInt(t, status, exp_ier, "qag(f3,31pt) oscill status")
		fc.neval = 0

		errqg = Qag(fc, 2.71, 0.3, 1e-14, 0.0, w.limit, w, &result, &abserr, Qk31)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, -exp_result, 1e-15, "qag(f3,31pt) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qag(f3,31pt) reverse abserr")
		test.TestInt(t, fc.neval, exp_neval, "qag(f3,31pt) reverse neval")
		test.TestInt(t, int(w.size), exp_last, "qag(f3,31pt) reverse last")
		test.TestInt(t, status, exp_ier, "qag(f3,31pt) reverse status")
	}

	/* Check the singularity detection (singularity at x=-0.1 in this example) */
	{
		w, _ := NewWorkspace(1000)
		exp_neval := 5151
		exp_ier := err.ESING
		exp_last := 51
		var result, abserr float64
		var status int

		alpha := 2.0
		f := f16{alpha: alpha}
		fc := &counter{f, 0}
		errqg := Qag(fc, -1.0, 1.0, 1e-14, 0.0, w.limit, w, &result, &abserr, Qk51)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestInt(t, fc.neval, exp_neval, "qag(f16,51pt) sing neval")
		test.TestInt(t, w.size, exp_last, "qag(f16,51pt) sing last")
		test.TestInt(t, status, exp_ier, "qag(f16,51pt) sing status")

		fc.neval = 0
		errqg = Qag(fc, 1.0, -1.0, 1e-14, 0.0, w.limit, w, &result, &abserr, Qk51)

		if errqg != nil {
			status = errqg.Status()
		}

		test.TestInt(t, fc.neval, exp_neval, "qag(f16,51pt) rev neval")
		test.TestInt(t, w.size, exp_last, "qag(f16,51pt) rev last")
		test.TestInt(t, status, exp_ier, "qag(f16,51pt) rev status")
	}

	/* Check for hitting the iteration limit */
	{
		w, _ := NewWorkspace(3)
		exp_result := 9.565151449233894709
		exp_abserr := 1.570369823891028460e+01
		exp_neval := 305
		exp_ier := err.EMAXITER
		exp_last := 3
		var result, abserr float64
		var status int

		a := []float64{-5.000000000000000000e-01,
			0.000000000000000000,
			-1.000000000000000000}
		b := []float64{0.000000000000000000,
			1.000000000000000000,
			-5.000000000000000000e-01}
		r := []float64{9.460353469435913709,
			9.090909090909091161e-02,
			1.388888888888888812e-02}
		e := []float64{1.570369823891028460e+01,
			1.009293658750142399e-15,
			1.541976423090495140e-16}
		order := []int{1, 2, 3}

		alpha := 1.0
		f := f16{alpha: alpha}
		fc := &counter{f, 0}
		errqg := Qag(fc, -1.0, 1.0, 1e-14, 0.0, w.limit, w, &result, &abserr, Qk61)
		if errqg != nil {
			status = errqg.Status()
		}
		test.TestRel(t, result, exp_result, 1e-15, "qag(f16,61pt) limit result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qag(f16,61pt) limit abserr")
		test.TestInt(t, fc.neval, exp_neval, "qag(f16,61pt) limit neval")
		test.TestInt(t, w.size, exp_last, "qag(f16,61pt) limit last")
		test.TestInt(t, status, exp_ier, "qag(f16,61pt) limit status")

		for i := 0; i < 3; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qag(f16,61pt) limit alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qag(f16,61pt) limit blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-15, "qag(f16,61pt) limit rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-6, "qag(f16,61pt) limit elist")
			test.TestInt(t, w.order[i], order[i]-1, "qag(f16,61pt) limit order")
		}

		fc.neval = 0
		errqg = Qag(fc, 1.0, -1.0, 1e-14, 0.0, w.limit, w, &result, &abserr, Qk61)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, -exp_result, 1e-15, "qag(f16,61pt) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qag(f16,61pt) reverse abserr")
		test.TestInt(t, fc.neval, exp_neval, "qag(f16,61pt) reverse neval")
		test.TestInt(t, w.size, exp_last, "qag(f16,61pt) reverse last")
		test.TestInt(t, status, exp_ier, "qag(f16,61pt) reverse status")
	}

	/*	QAGS  */

	/* Test the adaptive integrator with extrapolation QAGS */
	{
		var result, abserr float64
		var status int

		w, _ := NewWorkspace(1000)
		exp_result := 7.716049382715789440e-02
		exp_abserr := 2.216394961010438404e-12
		exp_neval := 189
		exp_ier := 0
		exp_last := 5

		a := []float64{0, 0.5, 0.25, 0.125, 0.0625}
		b := []float64{0.0625, 1, 0.5, 0.25, 0.125}
		r := []float64{3.919381915366914693e-05,
			5.491842501998223103e-02,
			1.909827770934243579e-02,
			2.776531175604360097e-03,
			3.280661030752062609e-04}
		e := []float64{2.215538742580964735e-12,
			6.097169993333454062e-16,
			2.120334764359736441e-16,
			3.082568839745514608e-17,
			3.642265412331439511e-18}
		order := []int{1, 2, 3, 4, 5}

		alpha := 2.6
		f := f1{alpha: alpha}
		fc := &counter{f, 0}

		errqg := Qags(fc, 0.0, 1.0, 0.0, 1e-10, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-15, "qags(f1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qags(f1) smooth abserr")
		test.TestInt(t, fc.neval, exp_neval, "qags(f1) smooth neval")
		test.TestInt(t, w.size, exp_last, "qags(f1) smooth last")
		test.TestInt(t, status, exp_ier, "qags(f1) smooth status")

		for i := 0; i < 5; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qags(f1) smooth alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qags(f1) smooth blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-15, "qags(f1) smooth rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-6, "qags(f1) smooth elist")
			test.TestInt(t, w.order[i], order[i]-1, "qags(f1) smooth order")
		}

		fc.neval = 0
		errqg = Qags(fc, 1.0, 0.0, 0.0, 1e-10, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, -exp_result, 1e-15, "qags(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qags(f1) reverse abserr")
		test.TestInt(t, fc.neval, exp_neval, "qags(f1) reverse neval")
		test.TestInt(t, w.size, exp_last, "qags(f1) reverse last")
		test.TestInt(t, status, exp_ier, "qags(f1) reverse status")
	}

	/* Test f11 using an absolute error bound */
	{
		var result, abserr float64
		var status int

		w, _ := NewWorkspace(1000)

		exp_result := -5.908755278982136588e+03
		exp_abserr := 1.299646281053874554e-10
		exp_neval := 357
		exp_ier := 0
		exp_last := 9

		a := []float64{1.000000000000000000e+00,
			5.005000000000000000e+02,
			2.507500000000000000e+02,
			1.258750000000000000e+02,
			6.343750000000000000e+01,
			3.221875000000000000e+01,
			1.660937500000000000e+01,
			8.804687500000000000e+00,
			4.902343750000000000e+00}
		b := []float64{4.902343750000000000e+00,
			1.000000000000000000e+03,
			5.005000000000000000e+02,
			2.507500000000000000e+02,
			1.258750000000000000e+02,
			6.343750000000000000e+01,
			3.221875000000000000e+01,
			1.660937500000000000e+01,
			8.804687500000000000e+00}
		r := []float64{-3.890977835520834649e+00,
			-3.297343675805121620e+03,
			-1.475904154146372775e+03,
			-6.517404019686431411e+02,
			-2.829354222635842007e+02,
			-1.201692001973227519e+02,
			-4.959999906099650246e+01,
			-1.971441499411640308e+01,
			-7.457032710459004399e+00}
		e := []float64{6.448276035006137169e-11,
			3.660786868980994028e-11,
			1.638582774073219226e-11,
			7.235772003440423011e-12,
			3.141214202790722909e-12,
			1.334146129098576244e-12,
			5.506706097890446534e-13,
			2.188739744348345039e-13,
			8.278969410534525339e-14}
		order := []int{1, 2, 3, 4, 5, 6, 7, 8, 9}
		alpha := 2.0
		f := f11{alpha: alpha}
		fc := &counter{f, 0}

		errqg := Qags(fc, 1.0, 1000.0, 1e-7, 0.0, w.limit, w, &result, &abserr)

		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-15, "qags(f11) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-3, "qags(f11) smooth abserr")
		test.TestInt(t, fc.neval, exp_neval, "qags(f11) smooth neval")
		test.TestInt(t, w.size, exp_last, "qags(f11) smooth last")
		test.TestInt(t, status, exp_ier, "qags(f11) smooth status")

		for i := 0; i < 9; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qags(f11) smooth alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qags(f11) smooth blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-15, "qags(f11) smooth rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-5, "qags(f11) smooth elist")
			test.TestInt(t, w.order[i], order[i]-1, "qags(f11) smooth order")
		}

		fc.neval = 0
		errqg = Qags(fc, 1000.0, 1.0, 1e-7, 0.0, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, -exp_result, 1e-15, "qags(f11) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-3, "qags(f11) reverse abserr")
		test.TestInt(t, fc.neval, exp_neval, "qags(f11) reverse neval")
		test.TestInt(t, w.size, exp_last, "qags(f11) reverse last")
		test.TestInt(t, status, exp_ier, "qags(f11) reverse status")
	}

	/* Test infinite range integral f455 using a relative error bound */

	{
		var result, abserr float64
		var status int
		w, _ := NewWorkspace(1000)
		exp_result := -3.616892186127022568e-01
		exp_abserr := 3.016716913328831851e-06
		exp_neval := 285
		exp_ier := 0
		exp_last := 10

		a := []float64{9.687500000000000000e-01,
			0.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			7.500000000000000000e-01,
			1.250000000000000000e-01,
			8.750000000000000000e-01,
			6.250000000000000000e-02,
			9.375000000000000000e-01,
			3.125000000000000000e-02}
		b := []float64{1.000000000000000000e+00,
			3.125000000000000000e-02,
			7.500000000000000000e-01,
			5.000000000000000000e-01,
			8.750000000000000000e-01,
			2.500000000000000000e-01,
			9.375000000000000000e-01,
			1.250000000000000000e-01,
			9.687500000000000000e-01,
			6.250000000000000000e-02}
		r := []float64{-1.390003415539725340e-01,
			1.429785306003466313e-03,
			-1.229943369113085765e-02,
			2.995321156568048898e-03,
			-4.980050133751051655e-02,
			2.785385934678596704e-03,
			-8.653752279614615461e-02,
			1.736218164975512294e-03,
			-8.398745675010892142e-02,
			1.041689192004495576e-03}
		e := []float64{2.395037249893453013e-02,
			2.161214992172538524e-04,
			5.720644840858777846e-14,
			3.325474514168701167e-17,
			3.147380432198176412e-14,
			3.092399597147240624e-17,
			9.607595030230581153e-16,
			1.927589382528252344e-17,
			9.324480826368044019e-16,
			1.156507325466566521e-17}
		order := []int{1, 2, 3, 5, 7, 9, 4, 6, 8, 10}

		f := f455{}
		fc := &counter{f, 0}

		errqg := Qagiu(fc, 0.0, 0.0, 1.0e-3, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qagiu(f455) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qagiu(f455) smooth abserr")
		test.TestInt(t, fc.neval, exp_neval, "qagiu(f455) smooth neval")
		test.TestInt(t, w.size, exp_last, "qagiu(f455) smooth last")
		test.TestInt(t, status, exp_ier, "qagiu(f455) smooth status")

		for i := 0; i < 10; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qagiu(f455) smooth alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qagiu(f455) smooth blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-15, "qagiu(f455) smooth rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-4, "qagiu(f455) smooth elist")
			test.TestInt(t, w.order[i], order[i]-1, "qagiu(f455) smooth order")
		}

	}

	/* Test infinite range integral f15 using a relative error bound */

	{
		var result, abserr float64
		var status int
		w, _ := NewWorkspace(1000)
		exp_result := 6.553600000000024738e+04
		exp_abserr := 7.121667111456009280e-04
		exp_neval := 285
		exp_ier := 0
		exp_last := 10

		a := []float64{0.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01,
			6.250000000000000000e-02,
			3.125000000000000000e-02,
			1.562500000000000000e-02,
			7.812500000000000000e-03,
			3.906250000000000000e-03,
			1.953125000000000000e-03}
		b := []float64{1.953125000000000000e-03,
			1.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01,
			6.250000000000000000e-02,
			3.125000000000000000e-02,
			1.562500000000000000e-02,
			7.812500000000000000e-03,
			3.906250000000000000e-03}
		r := []float64{1.099297665754340292e+00,
			3.256176475185617591e-01,
			8.064694554185326325e+00,
			8.873128656118993263e+01,
			6.977679035845269482e+02,
			4.096981198511257389e+03,
			1.574317583220441520e+04,
			2.899418134793237914e+04,
			1.498314766425578091e+04,
			9.225251570832365360e+02}
		e := []float64{7.101865971621337814e-04,
			1.912660677170175771e-08,
			9.167763417119923333e-08,
			3.769501719163865578e-07,
			6.973493131275552509e-07,
			1.205653952340679711e-07,
			1.380003928453846583e-07,
			1.934652413547325474e-07,
			3.408933028357320364e-07,
			2.132473175465897029e-09}
		order := []int{1, 5, 4, 9, 8, 7, 6, 3, 2, 10}

		alpha := 5.0

		f := f15{alpha}
		fc := &counter{f, 0}

		errqg := Qagiu(fc, 0.0, 0.0, 1.0e-7, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qagiu(f15) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qagiu(f15) smooth abserr")
		test.TestInt(t, fc.neval, exp_neval, "qagiu(f15) smooth neval")
		test.TestInt(t, w.size, exp_last, "qagiu(f15) smooth last")
		test.TestInt(t, status, exp_ier, "qagiu(f15) smooth status")

		for i := 0; i < 10; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qagiu(f15) smooth alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qagiu(f15) smooth blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-15, "qagiu(f15) smooth rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-4, "qagiu(f15) smooth elist")
			test.TestInt(t, w.order[i], order[i]-1, "qagiu(f15) smooth order")
		}

	}

	/* Test infinite range integral f16 using an absolute error bound */

	{
		var result, abserr float64
		var status int
		w, _ := NewWorkspace(1000)
		exp_result := 1.000000000006713292e-04
		exp_abserr := 3.084062020905636316e-09
		exp_neval := 165
		exp_ier := 0
		exp_last := 6

		a := []float64{0.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01,
			6.250000000000000000e-02,
			3.125000000000000000e-02}
		b := []float64{3.125000000000000000e-02,
			1.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01,
			6.250000000000000000e-02}
		r := []float64{7.633587786326674618e-05,
			9.900990099009899620e-07,
			1.922522349322310737e-06,
			3.629434715543053753e-06,
			6.501422186103209199e-06,
			1.062064387653501389e-05}
		e := []float64{3.084061858351569051e-09,
			3.112064814755089674e-17,
			4.543453652226561245e-17,
			4.908618166361344548e-17,
			3.014338672269481784e-17,
			6.795996738013555461e-18}
		order := []int{1, 4, 3, 2, 5, 6}
		alpha := 1.0
		f := f16{alpha}
		fc := &counter{f, 0}
		errqg := Qagiu(fc, 99.9, 1.0e-7, 0.0, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qagiu(f16) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qagiu(f16) smooth abserr")
		test.TestInt(t, fc.neval, exp_neval, "qagiu(f16) smooth neval")
		test.TestInt(t, w.size, exp_last, "qagiu(f16) smooth last")
		test.TestInt(t, status, exp_ier, "qagiu(f16) smooth status")

		for i := 0; i < 6; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qagiu(f16) smooth alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qagiu(f16) smooth blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-15, "qagiu(f16) smooth rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-4, "qagiu(f16) smooth elist")
			test.TestInt(t, w.order[i], order[i]-1, "qagiu(f16) smooth order")
		}

	}

	/* Test infinite range integral myfn1 using an absolute error bound */
	{
		var result, abserr float64
		var status int
		w, _ := NewWorkspace(1000)

		exp_result := 2.275875794468747770e+00
		exp_abserr := 7.436490118267390744e-09
		exp_neval := 270
		exp_ier := 0
		exp_last := 5

		a := []float64{1.250000000000000000e-01,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			0.000000000000000000e+00,
			3.750000000000000000e-01}
		b := []float64{2.500000000000000000e-01,
			1.000000000000000000e+00,
			3.750000000000000000e-01,
			1.250000000000000000e-01,
			5.000000000000000000e-01}
		r := []float64{4.639317228058405717e-04,
			1.691664195356748834e+00,
			1.146307471900291086e-01,
			4.379392477350953574e-20,
			4.691169201991640669e-01}
		e := []float64{3.169263960393051137e-09,
			4.265988974874425043e-09,
			1.231954072964969637e-12,
			8.360902986775307673e-20,
			5.208244060463541433e-15}
		order := []int{2, 1, 3, 5, 4}

		f := myfn1{}
		fc := &counter{f, 0}

		errqg := Qagi(fc, 1.0e-7, 0.0, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qagiu(myfn1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qagiu(myfn1) smooth abserr")
		test.TestInt(t, fc.neval, exp_neval, "qagiu(myfn1) smooth neval")
		test.TestInt(t, w.size, exp_last, "qagiu(myfn1) smooth last")
		test.TestInt(t, status, exp_ier, "qagiu(myfn1) smooth status")

		for i := 0; i < 5; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qagiu(myfn1) smooth alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qagiu(myfn1) smooth blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-14, "qagiu(myfn1) smooth rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-4, "qagiu(myfn1) smooth elist")
			test.TestInt(t, w.order[i], order[i]-1, "qagiu(myfn1) smooth order")
		}

	}

	/* Test infinite range integral myfn2 using an absolute error bound */

	{
		var result, abserr float64
		var status int
		w, _ := NewWorkspace(1000)
		exp_result := 2.718281828459044647e+00
		exp_abserr := 1.588185109253204805e-10
		exp_neval := 135
		exp_ier := 0
		exp_last := 5

		a := []float64{0.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01,
			6.250000000000000000e-02}
		b := []float64{6.250000000000000000e-02,
			1.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01}
		r := []float64{8.315287189746029816e-07,
			1.718281828459045091e+00,
			8.646647167633871867e-01,
			1.328565310599463256e-01,
			2.477920647947255521e-03}
		e := []float64{1.533437090413525935e-10,
			4.117868247943567505e-12,
			7.802455785301941044e-13,
			5.395586026138397182e-13,
			3.713312434866150125e-14}
		order := []int{1, 2, 3, 4, 5}
		alpha := 1.0
		f := myfn2{alpha}
		fc := &counter{f, 0}

		errqg := Qagil(fc, 1.0, 1.0e-7, 0.0, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qagiu(myfn2) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qagiu(myfn2) smooth abserr")
		test.TestInt(t, fc.neval, exp_neval, "qagiu(myfn2) smooth neval")
		test.TestInt(t, w.size, exp_last, "qagiu(myfn2) smooth last")
		test.TestInt(t, status, exp_ier, "qagiu(myfn2) smooth status")

		for i := 0; i < 5; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qagiu(myfn2) smooth alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qagiu(myfn2) smooth blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-14, "qagiu(myfn2) smooth rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-4, "qagiu(myfn2) smooth elist")
			test.TestInt(t, w.order[i], order[i]-1, "qagiu(myfn2) smooth order")
		}

	}

	/*	QAGP  */
	/* Test integral f454 with integrable singular points */

	{
		w, _ := NewWorkspace(1000)
		/* All results are for GSL_IEEE_MODE=double-precision */

		exp_result := 5.274080611672716401e+01
		exp_abserr := 1.755703848687062418e-04
		exp_neval := 777
		exp_ier := 0
		exp_last := 20
		var result, abserr float64
		var status int

		a := []float64{9.687500000000000000e-01,
			1.401269388548935790e+00,
			1.414213562373095145e+00,
			1.000000000000000000e+00,
			0.000000000000000000e+00,
			2.207106781186547462e+00,
			1.810660171779821415e+00,
			1.207106781186547462e+00,
			5.000000000000000000e-01,
			1.103553390593273731e+00,
			1.612436867076458391e+00,
			1.310660171779821415e+00,
			7.500000000000000000e-01,
			1.051776695296636976e+00,
			1.513325214724776657e+00,
			1.362436867076458391e+00,
			8.750000000000000000e-01,
			1.463769388548935790e+00,
			1.388325214724776657e+00,
			9.375000000000000000e-01}
		b := []float64{1.000000000000000000e+00,
			1.414213562373095145e+00,
			1.463769388548935790e+00,
			1.051776695296636976e+00,
			5.000000000000000000e-01,
			3.000000000000000000e+00,
			2.207106781186547462e+00,
			1.310660171779821415e+00,
			7.500000000000000000e-01,
			1.207106781186547462e+00,
			1.810660171779821415e+00,
			1.362436867076458391e+00,
			8.750000000000000000e-01,
			1.103553390593273731e+00,
			1.612436867076458391e+00,
			1.388325214724776657e+00,
			9.375000000000000000e-01,
			1.513325214724776657e+00,
			1.401269388548935790e+00,
			9.687500000000000000e-01}
		r := []float64{-1.125078814079027711e-01,
			-1.565132123531515207e-01,
			-4.225328513207429193e-01,
			-1.830392049835374568e-01,
			6.575875041899758092e-03,
			4.873920540843067783e+01,
			6.032891565603589079e+00,
			-2.991531901645863023e-01,
			-7.326282608704996063e-03,
			-2.431894410706912923e-01,
			5.911661670635662835e-01,
			-2.236786562536174916e-01,
			-5.647871991778510847e-02,
			-1.305470403178642658e-01,
			-1.721363984401322045e-01,
			-1.589345454585119055e-01,
			-7.406626263352669715e-02,
			-2.208730668000830344e-01,
			-1.048692749517999567e-01,
			-6.302287584527696551e-02}
		e := []float64{2.506431410088378817e-02,
			2.730454695485963826e-02,
			1.017446081816190118e-01,
			3.252808038935910834e-02,
			7.300687878575027348e-17,
			5.411138804637469780e-13,
			6.697855121200013106e-14,
			3.321267596107916554e-15,
			1.417509685426979386e-16,
			2.699945168224041491e-15,
			6.573952690524728748e-15,
			2.483331942899818875e-15,
			6.270397525408045936e-16,
			1.449363299575615261e-15,
			1.911097929242846383e-15,
			1.764527917763735212e-15,
			8.223007012367522077e-16,
			2.452183642810224359e-15,
			1.164282836272345215e-15,
			6.996944784151910810e-16}
		order := []int{3, 4, 2, 1, 6, 7, 11, 8, 10, 12, 18,
			15, 16, 14, 19, 17, 20, 13, 9, 5}

		f := f454{}
		fc := &counter{f, 0}

		pts := make([]float64, 4)

		pts[0] = 0.0
		pts[1] = 1.0
		pts[2] = math.Sqrt(2.0)
		pts[3] = 3.0

		errqg := Qagp(fc, pts, 4, 0.0, 1.0e-3, w.limit, w, &result, &abserr, Qk21)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qagp(f454) singular result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qagp(f454) singular abserr")
		test.TestInt(t, fc.neval, exp_neval, "qagp(f454) singular neval")
		test.TestInt(t, w.size, exp_last, "qagp(f454) singular last")
		test.TestInt(t, status, exp_ier, "qagp(f454) singular status")

		for i := 0; i < 20; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qagp(f454) singular alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qagp(f454) singular blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-14, "qagp(f454) singular rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-4, "qagp(f454) singular elist")
			test.TestInt(t, w.order[i], order[i]-1, "qagp(f454) singular order")
		}

	}
	/*	QAWC  */
	/* Test cauchy integration using a relative error bound */

	{
		var result, abserr float64
		var status int
		w, _ := NewWorkspace(1000)
		exp_result := -8.994400695837000137e-02
		exp_abserr := 1.185290176227023727e-06
		exp_neval := 215
		exp_ier := 0
		exp_last := 6

		a := []float64{-1.000000000000000000e+00,
			2.500000000000000000e+00,
			1.250000000000000000e+00,
			6.250000000000000000e-01,
			-5.000000000000000000e-01,
			-7.500000000000000000e-01}
		b := []float64{-7.500000000000000000e-01,
			5.000000000000000000e+00,
			2.500000000000000000e+00,
			1.250000000000000000e+00,
			6.250000000000000000e-01,
			-5.000000000000000000e-01}
		r := []float64{-1.234231128040012976e-01,
			3.579970394639702888e-03,
			2.249831615049339983e-02,
			7.214232992127905808e-02,
			2.079093855884046535e-02,
			-8.553244917962132821e-02}
		e := []float64{1.172832717970022565e-06,
			9.018232896137375412e-13,
			1.815172652101790755e-12,
			1.006998195150956048e-13,
			1.245463873006391609e-08,
			1.833082948207153514e-15}
		order := []int{1, 5, 3, 2, 4, 6}

		f := f459{}
		fc := &counter{f, 0}

		errqg := Qawc(fc, -1.0, 5.0, 0.0, 0.0, 1.0e-3, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qawc(f459) result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qawc(f459) abserr")
		test.TestInt(t, fc.neval, exp_neval, "qawc(f459) neval")
		test.TestInt(t, w.size, exp_last, "qawc(f459) last")
		test.TestInt(t, status, exp_ier, "qawc(f459) status")

		for i := 0; i < 6; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qawc(f459) alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qawc(f459) blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-14, "qawc(f459) rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-4, "qawc(f459) elist")
			test.TestInt(t, w.order[i], order[i]-1, "qawc(f459) order")
		}

		fc.neval = 0
		errqg = Qawc(fc, 5.0, -1.0, 0.0, 0.0, 1.0e-3, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, -exp_result, 1e-14, "qawc(f459) rev result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qawc(f459) rev abserr")
		test.TestInt(t, fc.neval, exp_neval, "qawc(f459) rev neval")
		test.TestInt(t, w.size, exp_last, "qawc(f459) rev last")
		test.TestInt(t, status, exp_ier, "qawc(f459) rev status")
	}

	/*	QAWS  */
	/* Test QAWS singular integration using a relative error bound */

	{
		var result, abserr float64
		var status int
		w, _ := NewWorkspace(1000)
		qt, _ := NewQAWSTable(0.0, 0.0, 1, 0)

		exp_result := -1.892751853489401670e-01
		exp_abserr := 1.129133712015747658e-08
		exp_neval := 280
		exp_ier := 0
		exp_last := 8

		a := []float64{0.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01,
			6.250000000000000000e-02,
			3.125000000000000000e-02,
			1.562500000000000000e-02,
			7.812500000000000000e-03}
		b := []float64{7.812500000000000000e-03,
			1.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01,
			6.250000000000000000e-02,
			3.125000000000000000e-02,
			1.562500000000000000e-02}
		r := []float64{-4.126317299834445824e-05,
			-1.076283950172247789e-01,
			-6.240573216173390947e-02,
			-1.456169844189576269e-02,
			-3.408925115926728436e-03,
			-8.914083918175634211e-04,
			-2.574191402137795482e-04,
			-8.034390712936630608e-05}
		e := []float64{1.129099387465713953e-08,
			3.423394967694403596e-13,
			6.928428071454762659e-16,
			1.616673288784094320e-16,
			3.784667152924835070e-17,
			9.896621209399419425e-18,
			2.857926564445496100e-18,
			8.919965558336773736e-19}
		order := []int{1, 2, 3, 4, 5, 6, 7, 8}

		f := &f458{}
		fc := &counter{f, 0}

		errqg := Qaws(fc, 0.0, 1.0, qt, 0.0, 1.0e-7, w.limit, w, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qaws(f458) ln(x-a) result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qaws(f458) ln(x-a) abserr")
		test.TestInt(t, fc.neval, exp_neval, "qaws(f458) ln(x-a) neval")
		test.TestInt(t, w.size, exp_last, "qaws(f458) ln(x-a) last")
		test.TestInt(t, status, exp_ier, "qaws(f458) ln(x-a) status")

		for i := 0; i < 6; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qaws(f458) ln(x-a) alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qaws(f458) ln(x-a) blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-14, "qaws(f458) ln(x-a) rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-4, "qaws(f458) ln(x-a) elist")
			test.TestInt(t, w.order[i], order[i]-1, "qaws(f458) ln(x-a) order")

		}

		/* Test without logs */
		qt.Set(-0.5, -0.3, 0, 0)

		errqg = Qaws(fc, 0.0, 1.0, qt, 0.0, 1.0e-7, w.limit, w, &result, &abserr)

		if errqg != nil {
			status = errqg.Status()
		}

		exp_result = 9.896686656601706433e-01
		exp_abserr = 5.888032513201251628e-08

		test.TestRel(t, result, exp_result, 1e-14, "qaws(f458) AB result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qaws(f458) AB abserr")

		/* Test with ln(x - a) */

		qt.Set(-0.5, -0.3, 1, 0)

		errqg = Qaws(fc, 0.0, 1.0, qt, 0.0, 1.0e-7, w.limit, w, &result, &abserr)

		if errqg != nil {
			status = errqg.Status()
		}

		exp_result = -3.636679470586539620e-01
		exp_abserr = 2.851348775257054093e-08

		test.TestRel(t, result, exp_result, 1e-14, "qaws(f458) AB ln(x-a) result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qaws(f458) AB ln(x-a) abserr")

		/* Test with ln(b - x) */
		qt.Set(-0.5, -0.3, 0, 1)

		errqg = Qaws(fc, 0.0, 1.0, qt, 0.0, 1.0e-7, w.limit, w, &result, &abserr)

		if errqg != nil {
			status = errqg.Status()
		}

		exp_result = -1.911489253363409802e+00
		exp_abserr = 9.854016753016499034e-09

		test.TestRel(t, result, exp_result, 1e-14, "qaws(f458) AB ln(b-x) result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qaws(f458) AB ln(b-x) abserr")

		/* Test with ln(x - a) ln(b - x) */

		qt.Set(-0.5, -0.3, 1, 1)

		errqg = Qaws(fc, 0.0, 1.0, qt, 0.0, 1.0e-7, w.limit, w, &result, &abserr)

		if errqg != nil {
			status = errqg.Status()
		}

		exp_result = 3.159922862811048172e-01
		exp_abserr = 2.336183482198144595e-08

		test.TestRel(t, result, exp_result, 1e-14, "qaws(f458) AB ln(x-a)ln(b-x) result")
		test.TestRel(t, abserr, exp_abserr, 1e-6, "qaws(f458) AB ln(x-a)ln(b-x) abserr")

	}
	/*	QAWO  */
	/* Test oscillatory integration using a relative error bound */

	{
		var result, abserr float64
		var status int
		w, _ := NewWorkspace(1000)
		wo, _ := NewQAWOTable(10.0*gsl.Pi, 1.0, GSL_INTEG_SINE, 1000)

		exp_result := -1.281368483991674190e-01
		exp_abserr := 6.875028324415666248e-12
		exp_neval := 305
		exp_ier := 0
		exp_last := 9

		a := []float64{0.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01,
			6.250000000000000000e-02,
			3.125000000000000000e-02,
			1.562500000000000000e-02,
			7.812500000000000000e-03,
			3.906250000000000000e-03}
		b := []float64{3.906250000000000000e-03,
			1.000000000000000000e+00,
			5.000000000000000000e-01,
			2.500000000000000000e-01,
			1.250000000000000000e-01,
			6.250000000000000000e-02,
			3.125000000000000000e-02,
			1.562500000000000000e-02,
			7.812500000000000000e-03}
		r := []float64{-1.447193692377651136e-03,
			2.190541162282139478e-02,
			-2.587726479625663753e-02,
			5.483209176363500886e-02,
			-3.081695575172510582e-02,
			-9.178321994387816929e-02,
			-3.886716016498160953e-02,
			-1.242306301902117854e-02,
			-3.659495117871544145e-03}
		e := []float64{8.326506625798146465e-07,
			1.302638552580516100e-13,
			7.259224351945759794e-15,
			1.249770395036711102e-14,
			7.832180081562836579e-16,
			1.018998440559284116e-15,
			4.315121611695628020e-16,
			1.379237060008662177e-16,
			4.062855738364339357e-17}
		order := []int{1, 2, 4, 3, 6, 5, 7, 8, 9}

		f := &f456{}
		fc := &counter{f, 0}

		errqg := Qawo(fc, 0.0, 0.0, 1e-7, w.limit, w, wo, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qawo(f456) result")
		test.TestRel(t, abserr, exp_abserr, 1e-3, "qawo(f456) abserr")
		test.TestInt(t, fc.neval, exp_neval, "qawo(f456) neval")
		test.TestInt(t, w.size, exp_last, "qawo(f456) last")
		test.TestInt(t, status, exp_ier, "qawo(f456) status")

		for i := 0; i < 9; i++ {
			test.TestRel(t, w.alist[i], a[i], 1e-15, "qawo(f456) alist")
			test.TestRel(t, w.blist[i], b[i], 1e-15, "qawo(f456) blist")
			test.TestRel(t, w.rlist[i], r[i], 1e-14, "qawo(f456) rlist")
			test.TestRel(t, w.elist[i], e[i], 1e-2, "qawo(f456) elist")
			test.TestInt(t, w.order[i], order[i]-1, "qawo(f456) order")
		}

		/* In reverse, flip limit and sign of length */
		wo.SetLength(-1.0)

		fc.neval = 0
		errqg = Qawo(fc, 1.0, 0.0, 1e-7, w.limit, w, wo, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, -exp_result, 1e-14, "qawo(f456) rev result")
		test.TestRel(t, abserr, exp_abserr, 1e-3, "qawo(f456) rev abserr")
		test.TestInt(t, fc.neval, exp_neval, "qawo(f456) rev neval")
		test.TestInt(t, w.size, exp_last, "qawo(f456) rev last")
		test.TestInt(t, status, exp_ier, "qawo(f456) rev status")

	}

	/*	QAWF  */
	/* Test fourier integration using an absolute error bound */

	{
		var result, abserr float64
		var status int
		w, _ := NewWorkspace(1000)
		wc, _ := NewWorkspace(1000)
		wo, _ := NewQAWOTable(gsl.Pi/2.0, 1.0, GSL_INTEG_COSINE, 1000)

		exp_result := 9.999999999279802765e-01
		exp_abserr := 1.556289974669056164e-08
		exp_neval := 590
		exp_ier := 0
		exp_last := 12

		r := []float64{1.013283128125232802e+00,
			-1.810857954748607349e-02,
			7.466754034900931897e-03,
			-4.360312526786496237e-03,
			2.950184068216192904e-03,
			-2.168238443073697373e-03,
			1.680910783140869081e-03,
			-1.352797860944863345e-03,
			1.119354921991485901e-03,
			-9.462367583691360827e-04,
			8.136341270731781887e-04,
			-7.093931338504278145e-04}
		e := []float64{1.224798040766472695e-12,
			1.396565155187268456e-13,
			1.053844511655910310e-16,
			6.505213034913026604e-19,
			7.155734338404329264e-18,
			1.105886215935214523e-17,
			9.757819552369539906e-18,
			5.854691731421723944e-18,
			4.553649124439220312e-18,
			7.643625316022806260e-18,
			2.439454888092388058e-17,
			2.130457268934021451e-17}

		f := &f457{}
		fc := &counter{f, 0}

		errqg := Qawf(fc, 0.0, 1e-7, w.limit, w, wc, wo, &result, &abserr)
		if errqg != nil {
			status = errqg.Status()
		}

		test.TestRel(t, result, exp_result, 1e-14, "qawf(f457) result")
		test.TestRel(t, abserr, exp_abserr, 1e-3, "qawf(f457) abserr")
		test.TestInt(t, fc.neval, exp_neval, "qawf(f457) neval")
		test.TestInt(t, w.size, exp_last, "qawf(f457) last")
		test.TestInt(t, status, exp_ier, "qawf(f457) status")

		for i := 0; i < 9; i++ {
			test.TestRel(t, w.rlist[i], r[i], 1e-12, "qawf(f457) rlist")
			/* We can only get within two orders of magnitude on the error
			   here, which is very sensitive to the floating point precision */
			test.TestRel(t, w.elist[i], e[i], 50.0, "qawf(f457) elist")
		}

	}
}
