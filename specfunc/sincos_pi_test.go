package specfunc

import (
	"fmt"
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"math"
	"strconv"
	"testing"
)

var (
	bigdbl = 2.0 / gsl.Float64Eps
)

func TestSinPi(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.5, 1, TEST_TOL0},
		{(-0.5), -1, TEST_TOL0},

		{1.5, -1, TEST_TOL0},
		{(-1.5), 1, TEST_TOL0},

		{2.5, 1, TEST_TOL0},
		{(-2.5), -1, TEST_TOL0},

		{3.5, -1, TEST_TOL0},
		{(-3.5), 1, TEST_TOL0},

		{0.375, 0.923879532511286756128183189397, TEST_TOL0},

		{(-0.375), -0.923879532511286756128183189397, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Sin_pi_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", c.x))
		})
	}

	fx := 0.0
	exact := 0.0
	ix := 0.0
	kmax := 12
	var k int

	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	fx = 0.5
	exact = 1.0

	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	fx = 0.03125
	exact = 0.0980171403295606019941955638886

	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	fx = 0.0625
	exact = 0.195090322016128267848284868477

	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	fx = 0.75
	exact = 0.707106781186547524400844362105

	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	fx = 0.0078125
	exact = 0.0245412285229122880317345294593

	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	var x float64
	fx = 0.0625
	exact = 0.195090322016128267848284868477
	ix = float64(math.MaxInt64) + 1.0
	ix += math.Abs(math.Mod(ix, 2.0)) /* make sure of even number */

	for k = 0; k < kmax; k++ {
		x = ix + fx
		x -= ix /* careful with compiler optimization */
		if (x != fx) || (math.Abs(ix+fx) >= bigdbl) {
			break
		}

		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix += 101.0
		exact = -exact
	}

	fx = -0.0625
	exact = -0.195090322016128267848284868477
	ix = float64(math.MaxInt64) - 1.0
	ix -= math.Abs(math.Mod(ix, 2.0)) /* make sure of even number */

	for k = 0; k < kmax; k++ {
		x = ix + fx
		x -= ix /* careful with compiler optimization */
		if (x != fx) || (math.Abs(ix+fx) >= bigdbl) {
			break
		}
		stat := Sin_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Sin_pi_e(%v)", ix+fx))
		ix -= 101.0
		exact = -exact
	}
}

func TestCosPi(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0, 1, TEST_TOL0},
		{1, -1, TEST_TOL0},
		{-1, -1, TEST_TOL0},
		{2, 1, TEST_TOL0},
		{-2, 1, TEST_TOL0},
		{3, -1, TEST_TOL0},
		{-3, -1, TEST_TOL0},
		{0.375, 0.382683432365089771728459984030, TEST_TOL0},
		{-0.375, 0.382683432365089771728459984030, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Cos_pi_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", c.x))
		})
	}

	fx := 0.0
	exact := 1.0
	ix := 0.0
	kmax := 12
	var k int

	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	fx = 0.5
	exact = 0

	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	fx = 0.0625
	exact = 0.980785280403230449126182236134

	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	fx = 0.4375
	exact = 0.195090322016128267848284868477

	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	fx = 0.4921875
	exact = 0.0245412285229122880317345294593

	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(3.0, float64(k+1))
		if k == 0 {
			exact = -exact
		}
	}

	exact = math.Abs(exact)
	ix = 0.0
	for k = 0; k < kmax; k++ {
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix = math.Pow(10.0, float64(k+1))
	}

	var x float64
	fx = 0.0625
	exact = 0.980785280403230449126182236134
	ix = float64(math.MaxInt64) + 1.0
	ix += math.Abs(math.Mod(ix, 2.0)) /* make sure of even number */

	for k = 0; k < kmax; k++ {
		x = ix + fx
		x -= ix /* careful with compiler optimization */
		if (x != fx) || (math.Abs(ix+fx) >= bigdbl) {
			break
		}

		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix += 101.0
		exact = -exact
	}

	fx = -0.0625
	exact = 0.980785280403230449126182236134
	ix = float64(math.MaxInt64) - 1.0
	ix -= math.Abs(math.Mod(ix, 2.0)) /* make sure of even number */

	for k = 0; k < kmax; k++ {
		x = ix + fx
		x -= ix /* careful with compiler optimization */
		if (x != fx) || (math.Abs(ix+fx) >= bigdbl) {
			break
		}
		stat := Cos_pi_e(ix+fx, r)
		run_test_sf(t, stat, r, exact, TEST_TOL0, err.SUCCESS, fmt.Sprintf("Cos_pi_e(%v)", ix+fx))
		ix -= 101.0
		exact = -exact
	}
}
