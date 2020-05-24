package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestConicalPhalf(t *testing.T) {
	r := new(Result)
	cases := []struct {
		lambda, x, expected, tol float64
	}{
		{0.0, -0.5, 0.8573827581049917129, TEST_TOL0},
		{0.0, 0.5, 0.8573827581049917129, TEST_TOL0},
		{0.0, 2.0, 0.6062611623284649811, TEST_TOL0},
		{0.0, 100.0, 0.07979045091636735635, TEST_TOL0},

		{10.0, -0.5, 5.345484922591867188e+08, TEST_TOL1},
		{10.0, 0.5, 15137.910380385258370, TEST_TOL1},
		{10.0, 2.0, 0.4992680691891618544, TEST_TOL1},
		{10.0, 100.0, -0.07272008163718195685, TEST_TOL2},

		{200.0, -1.0e-3, 1.3347639529084185010e+136, TEST_TOL2},
		{200.0, 1.0e-8, 1.0928098010940058507e+136, TEST_TOL2},
		{200.0, 0.5, 3.895546021611205442e+90, TEST_TOL2},
		{200.0, 10.0, -0.04308567180833581268, TEST_TOL3},
		{200.0, 100.0, -0.04694669186576399194, TEST_TOL3},
		{200.0, 1000.0, 0.023698140704121273277, TEST_TOL3},
		{200.0, 1.0e+8, -0.00006790983312124277891, TEST_TOL3},

		{1.0e+8, 1.1, 1.1599311133054742944, TEST_SQRT_TOL0},
		{1.0e+8, 100.0, 0.07971967557381557875, TEST_SQRT_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := ConicalP_half_e(c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("ConicalP_half_e(%v,%v)", c.lambda, c.x))

		})
	}
}

func TestConicalPmhalf(t *testing.T) {
	r := new(Result)
	cases := []struct {
		lambda, x, expected, tol float64
	}{
		{0.0, -0.5, 1.7956982494514644808, TEST_TOL0},
		{0.0, 0.5, 0.8978491247257322404, TEST_TOL0},
		{0.0, 2.0, 0.7984204253272901551, TEST_TOL0},
		{0.0, 100.0, 0.4227531369388072584, TEST_TOL0},

		{10.0, -0.5, 5.345484922591867181e+07, TEST_TOL1},
		{10.0, 0.5, 1513.7910356104985334, TEST_TOL1},
		{10.0, 2.0, 0.03439243987215615642, TEST_TOL1},
		{10.0, 100.0, 0.003283756665952609624, TEST_TOL2},

		{200.0, -0.5, 1.7699538115312304280e+179, TEST_TOL2},
		{200.0, 1.0e-8, 5.464049005470029253e+133, TEST_TOL2},
		{200.0, 0.5, 1.9477730108056027211e+88, TEST_TOL2},
		{200.0, 10.0, 0.0012462575917716355362, TEST_TOL2},
		{200.0, 100.0, -0.0003225881344802625149, TEST_TOL2},
		{200.0, 1000.0, -0.00004330652890886567623, TEST_TOL3},
		{200.0, 1.0e+8, 2.0943091278037078483e-07, TEST_TOL3},

		{1.0e+8, 1.1, 2.092320445620989618e-09, 16.0 * TEST_SQRT_TOL0},
		{1.0e+8, 100.0, -3.359967833599016923e-11, 256.0 * TEST_SQRT_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := ConicalP_mhalf_e(c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("ConicalP_mhalf_e(%v,%v)", c.lambda, c.x))

		})
	}
}

func TestConicalP0(t *testing.T) {
	r := new(Result)
	cases := []struct {
		lambda, x, expected, tol float64
	}{
		{0.0, -0.5, 1.3728805006183501647, TEST_TOL0},
		{0.0, 0.5, 1.0731820071493643751, TEST_TOL0},
		{0.0, 2.0, 0.9012862993604472987, TEST_TOL0},
		{0.0, 100.0, 0.30091748588199264556, TEST_TOL0},

		{10.0, -0.5, 1.6795592815421804669e+08, TEST_TOL1},
		{10.0, 0.5, 4826.034132009618240, TEST_TOL1},
		{10.0, 2.0, 0.18798468917758716146, TEST_TOL2},
		{10.0, 100.0, -0.008622130749987962529, TEST_TOL2},

		{200.0, -0.5, 2.502194818646823e+180, TEST_TOL4},

		{1000.0, 100.0, 0.0017908817653497715844, TEST_TOL3},
		{1000.0, 1000.0, -0.0006566893804926284301, TEST_TOL3},
		{1000.0, 1.0e+8, 2.3167213561756390068e-06, TEST_TOL4},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := ConicalP_0_e(c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("ConicalP_0_e(%v,%v)", c.lambda, c.x))

		})
	}
}

func TestConicalP1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		lambda, x, expected, tol float64
	}{
		{0.0, -0.5, 0.4939371126656998499, TEST_TOL1},
		{0.0, 0.5, 0.14933621085538265636, TEST_TOL1},
		{0.0, 2.0, -0.13666874968871549533, TEST_TOL1},
		{0.0, 100.0, -0.10544528203156629098, TEST_TOL2},

		{10.0, -0.5, 1.7253802958788312520e+09, TEST_TOL2},
		{10.0, 0.5, 46781.02294059967988, TEST_TOL1},
		{10.0, 2.0, 0.26613342643657444400, TEST_TOL2},
		{10.0, 100.0, -0.23281959695501029796, TEST_TOL2},

		/* FIXME: Mathematica gets some brain-damaged numbers for
		 * these x < 0 points. I have checked what I am doing in detail,
		 * and it must be right because you can do it by summing
		 * manifestly positive definite quantities.
		 */
		{200.0, -0.999, 2.71635193199341135e+270, TEST_TOL2},
		{200.0, -0.9, 4.2952493176812905e+234, TEST_TOL2},
		{200.0, -0.5, 5.01159205956053439e+182, TEST_TOL3},
		{200.0, 0.999, 195733.0396081538, TEST_TOL2},
		{200.0, 10.0, -2.9272610662414349553, TEST_TOL2},

		{1000.0, 100.0, -1.7783258105862399857, TEST_TOL6},
		{1000.0, 1000.0, 0.4535161075156427179, TEST_TOL4},
		{1000.0, 1.0e+8, 0.0009983414549874888478, TEST_SQRT_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := ConicalP_1_e(c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("ConicalP_1_e(%v,%v)", c.lambda, c.x))

		})
	}
}

func TestConicalPsphReg(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l                        int
		lambda, x, expected, tol float64
	}{
		{2, 1.0, -0.5, 1.6406279287008789526, TEST_TOL0},
		{10, 1.0, -0.5, 0.000029315266725049129448, TEST_TOL1},
		{20, 1.0, -0.5, 7.335769429462034431e-15, TEST_TOL1},
		{30, 1.0, -0.5, 1.3235612394267378871e-26, TEST_TOL2},
		{10, 1.0, 0.5, 2.7016087199857873954e-10, TEST_TOL1},
		{20, 1.0, 0.5, 1.1782569701435933399e-24, TEST_TOL1},
		{30, 1.0, 0.5, 3.636240588303797919e-41, TEST_TOL1},
		{10, 1.0, 2.0, 2.4934929626284934483e-10, TEST_TOL1},
		{20, 1.0, 2.0, 1.1284762488012616191e-24, TEST_TOL2},
		{30, 100.0, 100.0, -1.6757772087159526048e-64, TEST_TOL6},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := ConicalP_sph_reg_e(c.l, c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("ConicalP_sph_reg_e(%v,%v,%v)", c.l, c.lambda, c.x))

		})
	}
}

func TestConicalPCylReg(t *testing.T) {
	r := new(Result)
	cases := []struct {
		l                        int
		lambda, x, expected, tol float64
	}{
		{2, 1.0, -0.5, 2.2048510472375258708, TEST_TOL0},
		{10, 1.0, -0.5, 0.00007335034531618655690, TEST_TOL1},
		{20, 1.0, -0.5, 2.5419860619212164696e-14, TEST_TOL1},
		{30, 1.0, -0.5, 5.579714972260536827e-26, TEST_TOL2},
		{10, 1.0, 0.5, 1.1674078819646475282e-09, TEST_TOL0},
		{20, 1.0, 0.5, 7.066408031229072207e-24, TEST_TOL1},
		{30, 1.0, 0.5, 2.6541973286862588488e-40, TEST_TOL1},
		{10, 1.0, 2.0, 1.0736109751890863051e-09, TEST_TOL2},
		{20, 1.0, 2.0, 6.760965304863386741e-24, TEST_TOL2},
		{30, 100.0, 100.0, -4.268753482520651007e-63, TEST_TOL4},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := ConicalP_cyl_reg_e(c.l, c.lambda, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("ConicalP_cyl_reg_e(%v,%v,%v)", c.l, c.lambda, c.x))

		})
	}
}
