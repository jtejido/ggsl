package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestHyperg2f1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, b, c, x, expected float64
		tol                  float64
	}{
		{1, 1, 1, 0.5, 2.0, TEST_TOL0},
		{8, 8, 1, 0.5, 12451584.0, TEST_TOL0},
		{8, -8, 1, 0.5, 0.13671875, TEST_TOL0},
		{8, -8.1, 1, 0.5, 0.14147385378899930422, TEST_TOL4},
		{8, -8, 1, -0.5, 4945.136718750000000, TEST_TOL0},
		{8, -8, -5.5, 0.5, -906.6363636363636364, TEST_TOL0},
		{8, -8, -5.5, -0.5, 24565.363636363636364, TEST_TOL0},
		{8, 8, 1, -0.5, -0.006476312098196747669, TEST_TOL2},
		{8, 8, 5, 0.5, 4205.714285714285714, TEST_TOL0},
		{8, 8, 5, -0.5, 0.0028489656290296436616, TEST_TOL2},
		{9, 9, 1, 0.99, 1.2363536673577259280e+38, TEST_TOL2},
		{9, 9, -1.5, 0.99, 3.796186436458346579e+46, TEST_TOL2},
		{9, 9, -1.5, -0.99, 0.14733409946001025146, TEST_TOL1},
		{9, 9, -8.5, 0.99, -1.1301780432998743440e+65, TEST_TOL2},
		{9, 9, -8.5, -0.99, -8.856462606575344483, TEST_TOL1},
		{9, 9, -21.5, 0.99, 2.0712920991876073253e+95, TEST_TOL3},
		{9, 9, -21.5, -0.99, -74.30517015382249216, TEST_TOL2},
		{9, 9, -100.5, 0.99, -3.186778061428268980e+262, TEST_TOL3},
		{9, 9, -100.5, -0.99, 2.4454358338375677520, TEST_TOL1},
		{25, 25, 1, -0.5, -2.9995530823639545027e-06, TEST_SQRT_TOL0},

		{1.5, 0.5, 2.0, 1.0 - 1.0/64.0, 3.17175539044729373926, TEST_TOL3},
		{1.5, 0.5, 2.0, 1.0 - 1.0/128.0, 3.59937243502024563424, TEST_TOL2},
		{1.5, 0.5, 2.0, 1.0 - 1.0/256.0, 4.03259299524392504369, TEST_TOL1},
		{1.5, 0.5, 2.0, 1.0 - 1.0/1024.0, 4.90784159359675398250, TEST_TOL1},
		{1.5, 0.5, 2.0, 1.0 - 1.0/65536.0, 7.552266033399683914, TEST_TOL1},
		{1.5, 0.5, 2.0, 1.0 - 1.0/16777216.0, 11.08235454026043830363, TEST_TOL1},

		{1.5, 0.5, 2.0, -1.0 + 1.0/1024.0, 0.762910940909954974527, TEST_TOL0},
		{1.5, 0.5, 2.0, -1.0 + 1.0/65536.0, 0.762762124908845424449, TEST_TOL0},
		{1.5, 0.5, 2.0, -1.0 + 1.0/1048576.0, 0.762759911089064738044, TEST_TOL0},

		/* added special handling with x == 1.0 , Richard J. Mathar, 2008-01-09 */

		{1.5, 0.5, 3.0, 1.0, 1.6976527263135502482014268, TEST_TOL2},
		{1.5, -4.2, 3.0, 1.0, .15583601560025710649555254, TEST_TOL2},
		{-7.4, 0.7, -1.5, 1.0, -.34478866959246584996859, TEST_TOL2},
		{0.1, -2.7, -1.5, 1.0, 1.059766766063610122925, TEST_TOL2},

		/* Taylor Binnington a = 0 */

		{0, -2, -4, 0.5, 1.0, TEST_TOL2},

		/* Andrew Benson <abenson@caltech.edu> bug #24812
		   in Pari:
		   poch(a,x) = { gamma(a+x)/gamma(a) }
		   t(a,b,c,x,k) = { (poch(a,k)*poch(b,k)/poch(c,k)) * (x^k)/(k!) }
		   suminf(k=0,t(-10.34, 2.05, 3.05, 0.1725,k))  */

		{-10.34, 2.05, 3.05, 0.1725, 0.310473552213130010351006093079548, TEST_TOL2},
		{-9.99999999999, 2.05, 3.05, 0.1725, 0.32141934630197487540298837643890, TEST_TOL2},

		/* Didier Pinchon also bug #24812 */
		{11, -1, 11.0 / 2.0, 0.125, 0.75, TEST_TOL2},

		/* Bill Maier - bug #45926 */
		{-0.2, 8.8, 10.0, 0.8, 0.77998971427681563, TEST_TOL1},
		{-0.2, 9.8, 11.0, 0.8, 0.77574573497387267, TEST_TOL0},

		// #if 0 /* XXX - bug #39056 */
		/* Test case from Hatef Monajemi <monajemi@stanford.edu> */

		//  {3.5, -0.5, 5.0, 0.9,  0.5923981284370653465208973272, TEST_TOL2},

		/* Test case from Robert L Wolpert <Wolpert@stat.duke.edu> */

		//  {-1.0, -10.0, 1.0, 0.5,  6.0, TEST_TOL0},

		/* Test case from ldnlwm@163.com */

		//  {3.23191, -4.0229, 8.02291, 0.5,  0.4300243900348170646, TEST_TOL2},
		// #endif
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Hyperg_2F1_e(c.a, c.b, c.c, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Hyperg_2F1_e(%v,%v,%v,%v)", c.a, c.b, c.c, c.x))
		})
	}
}

func TestHyperg2F1Conj(t *testing.T) {
	r := new(Result)
	cases := []struct {
		aR, aI, c, x, expected float64
		tol                    float64
	}{
		{1, 1, 1, 0.5, 3.352857095662929028, TEST_TOL0},
		{8, 8, 1, 0.5, 1.7078067538891293983e+09, TEST_TOL0},
		{8, 8, 5, 0.5, 285767.15696901140627, TEST_TOL1},
		{8, 8, 1, -0.5, 0.007248196261471276276, TEST_TOL3},
		{8, 8, 5, -0.5, 0.00023301916814505902809, TEST_TOL3},
		{25, 25, 1, -0.5, 5.1696944096e-06, TEST_SQRT_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Hyperg_2F1_conj_e(c.aR, c.aI, c.c, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Hyperg_2F1_conj_e(%v,%v,%v,%v)", c.aR, c.aI, c.c, c.x))
		})
	}
}

func TestHyperg2F1Renorm(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, b, c, x, expected float64
		tol                  float64
	}{
		{1, 1, 1, 0.5, 2.0, TEST_TOL0},
		{8, 8, 1, 0.5, 12451584.0, TEST_TOL0},
		{8, -8, 1, 0.5, 0.13671875, TEST_TOL0},
		{8, -8, 1, -0.5, 4945.13671875, TEST_TOL0},
		{8, -8, -5.5, 0.5, -83081.19167659493609, TEST_TOL2},
		{8, -8, -5.5, -0.5, 2.2510895952730178518e+06, TEST_TOL2},
		{8, 8, 5, 0.5, 175.2380952380952381, TEST_TOL1},
		{9, 9, -1.5, 0.99, 1.6063266334913066551e+46, TEST_TOL2},
		{9, 9, -1.5, -0.99, 0.06234327316254516616, TEST_TOL2},
		{5, 5, -1, 0.5, 4949760.0, TEST_TOL1},
		{5, 5, -10, 0.5, 139408493229637632000.0, TEST_TOL2},
		{5, 5, -100, 0.5, 3.0200107544594411315e+206, TEST_TOL3},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Hyperg_2F1_renorm_e(c.a, c.b, c.c, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Hyperg_2F1_renorm_e(%v,%v,%v,%v)", c.a, c.b, c.c, c.x))
		})
	}
}

func TestHyperg2F1ConjRenorm(t *testing.T) {
	r := new(Result)
	cases := []struct {
		aR, aI, c, x, expected float64
		tol                    float64
	}{
		{9, 9, -1.5, 0.99, 5.912269095984229412e+49, TEST_TOL2},
		{9, 9, -1.5, -0.99, 0.10834020229476124874, TEST_TOL2},
		{5, 5, -1, 0.5, 1.4885106335357933625e+08, TEST_TOL2},
		{5, 5, -10, 0.5, 7.968479361426355095e+21, TEST_TOL2},
		{5, 5, -100, 0.5, 3.1113180227052313057e+208, TEST_TOL3},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Hyperg_2F1_conj_renorm_e(c.aR, c.aI, c.c, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Hyperg_2F1_conj_renorm_e(%v,%v,%v,%v)", c.aR, c.aI, c.c, c.x))
		})
	}
}
