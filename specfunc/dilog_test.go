package specfunc

import (
	"fmt"
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestDilog(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{-3.0, -1.9393754207667089531, TEST_TOL0},
		{-0.5, -0.4484142069236462024, TEST_TOL0},
		{-0.001, -0.0009997501110486510834, TEST_TOL0},
		{0.1, 0.1026177910993911, TEST_TOL0},
		{0.7, 0.8893776242860387386, TEST_TOL0},
		{1.0, 1.6449340668482260, TEST_TOL0},
		{1.5, 2.3743952702724802007, TEST_TOL0},
		{2.0, 2.4674011002723397, TEST_TOL0},
		{5.0, 1.7837191612666306277, TEST_TOL0},
		{11.0, 0.3218540439999117111, TEST_TOL1},
		{12.59, 0.0010060918167266208634, TEST_TOL3},
		{12.595, 0.00003314826006436236810, TEST_TOL5},
		{13.0, -0.07806971248458575855, TEST_TOL2},
		{20.0, -1.2479770861745251168, TEST_TOL2},
		{150.0, -9.270042702348657270, TEST_TOL0},
		{1100.0, -21.232504073931749553, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Dilog_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Dilog_e(%v)", c.x))
		})
	}
}

func TestComplexDilog(t *testing.T) {
	r1, r2 := new(Result), new(Result)
	cases := []struct {
		r, theta, expected1, tol1, expected2, tol2 float64
	}{
		{0.99999, gsl.Pi / 2.0, -0.20561329262779687646, TEST_TOL0,
			0.91595774018131512060, TEST_TOL0},

		{0.991, gsl.Pi / 2.0, -0.20250384721077806127, TEST_TOL0,
			0.90888544355846447810, TEST_TOL0},

		{0.98, gsl.Pi / 2.0, -0.19871638377785918403, TEST_TOL2,
			0.90020045882981847610, TEST_TOL2},

		{0.98, -gsl.Pi / 2.0, -0.19871638377785918403, TEST_TOL2,
			-0.90020045882981847610, TEST_TOL2},

		{0.95, gsl.Pi / 2.0, -0.18848636456893572091, TEST_TOL1,
			0.87633754133420277830, TEST_TOL1},

		{0.8, gsl.Pi / 2.0, -0.13980800855429037810, TEST_TOL0,
			0.75310609092419884460, TEST_TOL0},

		{0.8, -gsl.Pi / 2.0, -0.13980800855429037810, TEST_TOL0,
			-0.75310609092419884460, TEST_TOL0},

		{0.5, gsl.Pi / 2.0, -0.05897507442156586346, TEST_TOL1,
			0.48722235829452235710, TEST_TOL1},

		{0.5, -gsl.Pi / 2.0, -0.05897507442156586346, TEST_TOL1,
			-0.48722235829452235710, TEST_TOL1},

		{0.01, gsl.Pi / 2.0, -0.000024999375027776215378, TEST_TOL3,
			0.009999888892888684820, TEST_TOL3},

		{0.01, -gsl.Pi / 2.0, -0.000024999375027776215378, TEST_TOL3,
			-0.009999888892888684820, TEST_TOL3},

		{0.99, gsl.Pi / 4.0, 0.56273366219795547757, TEST_TOL3,
			0.97009284079274560384, TEST_TOL3},

		{0.99, -gsl.Pi / 4.0, 0.56273366219795547757, TEST_TOL3,
			-0.97009284079274560384, TEST_TOL3},

		{0.99, 3.0 * gsl.Pi / 4.0, -0.66210902664245926235, TEST_TOL1,
			0.51995305609998319025, TEST_TOL1},

		{0.99, 5.0 * gsl.Pi / 4.0, -0.66210902664245926235, TEST_TOL1,
			-0.51995305609998319025, TEST_TOL1},

		{0.99, 3.0 * gsl.Pi / 2.0, -0.20215874509123277909, TEST_TOL1,
			-0.90809733095648731408, TEST_TOL1},

		{0.25, 3.0 * gsl.Pi / 2.0, -0.01538741178141053563, TEST_TOL1,
			-0.24830175098230686908, TEST_TOL1},

		{0.25, 15.0 / 8.0 * gsl.Pi, 0.24266162342377302235, TEST_TOL1,
			-0.10860883369274445067, TEST_TOL1},

		{0.99, gsl.Pi / 8.0, 1.0571539648820244720, TEST_TOL0,
			0.7469145254610851318, TEST_TOL0},

		{0.99, gsl.Pi / 64.0, 1.5381800285902999666, TEST_TOL0,
			0.1825271634987756651, TEST_TOL0},

		{0.99, -gsl.Pi / 8.0, 1.05715396488202447202, TEST_TOL1,
			-0.74691452546108513176, TEST_TOL1},

		{1.00001, gsl.Pi / 2.0, -0.20562022409960237363, TEST_TOL1,
			0.91597344814458309320, TEST_TOL1},

		{10.0, gsl.Pi / 2.0, -3.0596887943287347304, TEST_TOL0,
			3.7167814930680685900, TEST_TOL0},

		{100.0, gsl.Pi / 2.0, -11.015004738293824854, TEST_TOL0,
			7.2437843013083534970, TEST_TOL0},

		/** tests brought up by Jim McElwaine bug report */

		{1.1, -gsl.Pi / 2.0, -0.24099184177382733037, TEST_TOL1,
			-0.99309132538137822631, TEST_TOL1},

		{1.1, 3.0 * gsl.Pi / 2.0, -0.24099184177382733037, TEST_TOL1,
			-0.99309132538137822631, TEST_TOL1},

		{1.1, -3.0 * gsl.Pi / 2.0, -0.24099184177382733037, TEST_TOL1,
			0.99309132538137822631, TEST_TOL1},

		{1.1, -gsl.Pi - 0.25*gsl.Pi, -0.72908565537087935118, TEST_TOL1,
			0.56225783937234862649, TEST_TOL1},

		{1.1, gsl.Pi + 0.25*gsl.Pi, -0.72908565537087935118, TEST_TOL1,
			-0.56225783937234862649, TEST_TOL1},

		{1.1, -gsl.Pi / 128.0, 1.8881719454909716580, TEST_TOL1,
			-0.3556738764969238976, TEST_TOL1},

		{1.1, gsl.Pi / 128.0, 1.8881719454909716580, TEST_TOL1,
			0.3556738764969238976, TEST_TOL1},

		{1.5, gsl.Pi / 8.0, 1.3498525763442498343, TEST_TOL1,
			1.4976532712229749493, TEST_TOL1},

		{1.5, -gsl.Pi / 8.0, 1.3498525763442498343, TEST_TOL1,
			-1.4976532712229749493, TEST_TOL1},

		{1.5, 2.0*gsl.Pi + gsl.Pi/8.0, 1.3498525763442498343, TEST_TOL1,
			1.4976532712229749493, TEST_TOL1},

		{1.5, 2.0*gsl.Pi - gsl.Pi/8.0, 1.3498525763442498343, TEST_TOL1,
			-1.4976532712229749493, TEST_TOL1},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Complex_dilog_e(c.r, c.theta, r1, r2)
			run_test_sf_2(t, stat, r1, c.expected1, c.tol1, r2, c.expected2, c.tol2, err.SUCCESS, fmt.Sprintf("Complex_dilog_e(%v,%v)", c.r, c.theta))
		})
	}
}

func TestComplexDilogXY(t *testing.T) {
	r1, r2 := new(Result), new(Result)
	cases := []struct {
		x, y, expected1, tol1, expected2, tol2 float64
	}{
		{0.0, 0.5, -0.05897507442156586346, TEST_TOL1,
			0.48722235829452235710, TEST_TOL1},

		{0.0, -0.5, -0.05897507442156586346, TEST_TOL1,
			-0.48722235829452235710, TEST_TOL1},

		{0.91464073718617389108, 0.37885659804143889673, 1.0571539648820244720, TEST_TOL0,
			0.7469145254610851318, TEST_TOL0},

		{0.91464073718617389108, -0.37885659804143889673, 1.05715396488202447202, TEST_TOL1,
			-0.74691452546108513176, TEST_TOL1},

		{-1.5, 0.0, -1.1473806603755707541, TEST_TOL1,
			0.0, TEST_TOL1},

		{0.5, 0.0, 0.58224052646501250590, TEST_TOL1,
			0.0, TEST_TOL1},

		{1.5, 0.0, 2.3743952702724802007, TEST_TOL1,
			-1.2738062049196005309, TEST_TOL1},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Complex_dilog_xy_e(c.x, c.y, r1, r2)
			run_test_sf_2(t, stat, r1, c.expected1, c.tol1, r2, c.expected2, c.tol2, err.SUCCESS, fmt.Sprintf("Complex_dilog_xy_e(%v,%v)", c.x, c.y))
		})
	}
}

func TestComplexSpenceXY(t *testing.T) {
	r1, r2 := new(Result), new(Result)
	cases := []struct {
		x, y, expected1, tol1, expected2, tol2 float64
	}{

		{1.5, 0.0, -0.44841420692364620244, TEST_TOL1,
			0.0, TEST_TOL1},

		{0.5, 0.0, 0.58224052646501250590, TEST_TOL1,
			0.0, TEST_TOL1},

		{0.0, 0.0, 1.6449340668482264365, TEST_TOL1,
			0.0, TEST_TOL1},

		{-0.5, 0.0, 2.3743952702724802007, TEST_TOL1,
			-1.2738062049196005309, TEST_TOL1},

		{-0.5, 1.0 / 1024.0, 2.3723507455234125018, TEST_TOL1,
			-1.2742581376517839070, TEST_TOL1},

		{-0.5, -1.0 / 1024.0, 2.3723507455234125018, TEST_TOL1,
			1.2742581376517839070, TEST_TOL1},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Complex_spence_xy_e(c.x, c.y, r1, r2)
			run_test_sf_2(t, stat, r1, c.expected1, c.tol1, r2, c.expected2, c.tol2, err.SUCCESS, fmt.Sprintf("Complex_spence_xy_e(%v,%v)", c.x, c.y))
		})
	}
}
