package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"math"
	"strconv"
	"testing"
)

func TestCoupling3j(t *testing.T) {
	r := new(Result)
	cases := []struct {
		two_ja, two_jb, two_jc, two_ma, two_mb, two_mc int
		expected, tol                                  float64
		status                                         int
	}{
		{0, 1, 1, 0, 1, -1, math.Sqrt(1.0 / 2.0), TEST_TOL0, err.SUCCESS},
		{1, 1, 2, 1, -1, 0, math.Sqrt(1.0 / 6.0), TEST_TOL0, err.SUCCESS},
		{2, 4, 6, 0, 2, -2, math.Sqrt(8.0 / 105.0), TEST_TOL0, err.SUCCESS},
		{4, 4, 8, 0, 0, 0, math.Sqrt(2.0 / 35.0), TEST_TOL0, err.SUCCESS},
		{4, 4, 8, 2, -2, 0, 2.0 / 3.0 * math.Sqrt(2.0/35.0), TEST_TOL2, err.SUCCESS},
		{4, 4, 8, 4, -4, 0, 1.0 / (3.0 * math.Sqrt(70.0)), TEST_TOL2, err.SUCCESS},

		/* Test 3j error checking */

		{-1, 1, 2, 1, -1, 0, math.NaN(), math.NaN(), err.EDOM},
		{1, -1, 2, 1, -1, 0, math.NaN(), math.NaN(), err.EDOM},
		{1, 1, -2, 1, -1, 0, math.NaN(), math.NaN(), err.EDOM},

		/* Test |m_i|<=j_i */

		{1, 1, 2, 2, -1, 0, 0, 0, err.SUCCESS},
		{1, 1, 2, 1, -2, 0, 0, 0, err.SUCCESS},
		{1, 1, 2, 1, -1, 3, 0, 0, err.SUCCESS},

		/* Test triangle condition j1 + j2 >= j, j >= j2 - j1, j>= j1 - j2 */

		{1, 1, 3, 1, -1, 0, 0, 0, err.SUCCESS},
		{1, 4, 2, 1, -1, 0, 0, 0, err.SUCCESS},
		{4, 1, 2, 1, -1, 0, 0, 0, err.SUCCESS},

		/* Test m1=m2=m3=0 with j1+j2+j3=odd*/

		{2 * 13, 2 * 13, 2 * 13, 0, 0, 0, 0, 0, err.SUCCESS},
		{2 * 2, 2 * 17, 2 * 18, 0, 0, 0, 0, 0, err.SUCCESS},
		{2 * 203, 2 * 203, 2 * 203, 0, 0, 0, 0, 0, err.SUCCESS},

		/* Test l1=249 l2=248, l3=2, m1=5, m2=-6, m3=1 */
		{2 * 249.0, 2 * 248.0, 2 * 2.0, 2 * 5.0, 2 * (-6.0), 2 * 1.0, 0.0228787564223517967033998, TEST_TOL3, err.SUCCESS},
		{2 * 248.0, 2 * 247.0, 2 * 2.0, 2 * 5.0, 2 * (-6.0), 2 * 1.0, -0.022926660587726369939271424097, TEST_TOL3, err.SUCCESS},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Coupling_3j_e(c.two_ja, c.two_jb, c.two_jc, c.two_ma, c.two_mb, c.two_mc, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Coupling_3j_e(%v,%v,%v,%v,%v,%v)", c.two_ja, c.two_jb, c.two_jc, c.two_ma, c.two_mb, c.two_mc))
		})
	}
}

func TestCoupling6j(t *testing.T) {
	r := new(Result)
	cases := []struct {
		two_ja, two_jb, two_jc, two_jd, two_je, two_jf int
		expected, tol                                  float64
		status                                         int
	}{
		{2, 2, 4, 2, 2, 2, 1.0 / 6.0, TEST_TOL0, err.SUCCESS},
		{4, 4, 2, 4, 4, 4, -1.0 / 10.0, TEST_TOL0, err.SUCCESS},
		{4, 4, 2, 4, 4, 2, 1.0 / 6.0, TEST_TOL0, err.SUCCESS},
		{4, 4, 2, 2, 2, 2, -0.5 / math.Sqrt(5.0), TEST_TOL0, err.SUCCESS},
		{4, 4, 4, 2, 2, 2, math.Sqrt(7.0/3.0) / 10.0, TEST_TOL0, err.SUCCESS},
		{6, 6, 6, 4, 4, 4, -math.Sqrt(3.0/5.0) / 14.0, TEST_TOL0, err.SUCCESS},
		{6, 6, 6, 4, 4, 2, -math.Sqrt(3.0/5.0) / 7.0, TEST_TOL0, err.SUCCESS},
		{1, 0, 1, 0, 1, 0, -math.Sqrt(1.0 / 2.0), TEST_TOL0, err.SUCCESS},
		{1, 0, 1, 1, 0, 1, -1.0 / 2.0, TEST_TOL0, err.SUCCESS},

		/* Test 6j error checking */

		{-2, 2, 4, 2, 2, 2, math.NaN(), math.NaN(), err.EDOM},
		{2, -2, 4, 2, 2, 2, math.NaN(), math.NaN(), err.EDOM},
		{2, 2, -4, 2, 2, 2, math.NaN(), math.NaN(), err.EDOM},
		{2, 2, 4, -2, 2, 2, math.NaN(), math.NaN(), err.EDOM},
		{2, 2, 4, 2, -2, 2, math.NaN(), math.NaN(), err.EDOM},
		{2, 2, 4, 2, 2, -2, math.NaN(), math.NaN(), err.EDOM},

		/* Test 6j triangle conditions */

		{2, 2, 4, 2, 2, 7, 0, 0, err.SUCCESS},
		{2, 2, 4, 2, 7, 2, 0, 0, err.SUCCESS},
		{2, 2, 4, 7, 2, 2, 0, 0, err.SUCCESS},
		{2, 2, 7, 2, 2, 2, 0, 0, err.SUCCESS},
		{2, 7, 4, 2, 2, 2, 0, 0, err.SUCCESS},
		{7, 2, 4, 2, 2, 2, 0, 0, err.SUCCESS},

		/* Test 6j half-integer/integer coupling conditions */

		{0, 2, 2, 44, 43, 43, 0, 0, err.SUCCESS},
		{1, 1, 1, 0, 1, 1, 0, 0, err.SUCCESS},
		{1, 1, 1, 1, 0, 1, 0, 0, err.SUCCESS},
		{1, 1, 1, 1, 1, 0, 0, 0, err.SUCCESS},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Coupling_6j_e(c.two_ja, c.two_jb, c.two_jc, c.two_jd, c.two_je, c.two_jf, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Coupling_6j_e(%v,%v,%v,%v,%v,%v)", c.two_ja, c.two_jb, c.two_jc, c.two_jd, c.two_je, c.two_jf))
		})
	}
}

func TestCoupling9j(t *testing.T) {
	r := new(Result)
	cases := []struct {
		two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji int
		expected, tol                                                          float64
		status                                                                 int
	}{
		{4, 2, 4, 3, 3, 2, 1, 1, 2, -math.Sqrt(1.0/6.0) / 10.0, TEST_TOL2, err.SUCCESS},
		{8, 4, 10, 7, 3, 8, 1, 1, 2, math.Sqrt(7.0/3.0) / 60.0, TEST_TOL2, err.SUCCESS},

		/* Test 9j error checking */

		{-4, 2, 4, 3, 3, 2, 1, 1, 2, math.NaN(), math.NaN(), err.EDOM},
		{4, -2, 4, 3, 3, 2, 1, 1, 2, math.NaN(), math.NaN(), err.EDOM},
		{4, 2, -4, 3, 3, 2, 1, 1, 2, math.NaN(), math.NaN(), err.EDOM},
		{4, 2, 4, -3, 3, 2, 1, 1, 2, math.NaN(), math.NaN(), err.EDOM},
		{4, 2, 4, 3, -3, 2, 1, 1, 2, math.NaN(), math.NaN(), err.EDOM},
		{4, 2, 4, 3, 3, -2, 1, 1, 2, math.NaN(), math.NaN(), err.EDOM},
		{4, 2, 4, 3, 3, 2, -1, 1, 2, math.NaN(), math.NaN(), err.EDOM},
		{4, 2, 4, 3, 3, 2, 1, -1, 2, math.NaN(), math.NaN(), err.EDOM},
		{4, 2, 4, 3, 3, 2, 1, 1, -2, math.NaN(), math.NaN(), err.EDOM},

		{10, 2, 4, 3, 3, 2, 1, 1, 2, 0, 0, err.SUCCESS},
		{4, 10, 4, 3, 3, 2, 1, 1, 2, 0, 0, err.SUCCESS},
		{4, 2, 10, 3, 3, 2, 1, 1, 2, 0, 0, err.SUCCESS},
		{4, 2, 4, 10, 3, 2, 1, 1, 2, 0, 0, err.SUCCESS},
		{4, 2, 4, 3, 10, 2, 1, 1, 2, 0, 0, err.SUCCESS},
		{4, 2, 4, 3, 3, 10, 1, 1, 2, 0, 0, err.SUCCESS},
		{4, 2, 4, 3, 3, 2, 10, 1, 2, 0, 0, err.SUCCESS},
		{4, 2, 4, 3, 3, 2, 1, 10, 2, 0, 0, err.SUCCESS},
		{4, 2, 4, 3, 3, 2, 1, 1, 10, 0, 0, err.SUCCESS},

		/* Test 9j half-integer/integer coupling conditions */

		{1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, err.SUCCESS},
		{1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, err.SUCCESS},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Coupling_9j_e(c.two_ja, c.two_jb, c.two_jc, c.two_jd, c.two_je, c.two_jf, c.two_jg, c.two_jh, c.two_ji, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Coupling_9j_e(%v,%v,%v,%v,%v,%v,%v,%v,%v)", c.two_ja, c.two_jb, c.two_jc, c.two_jd, c.two_je, c.two_jf, c.two_jg, c.two_jh, c.two_ji))
		})
	}
}
