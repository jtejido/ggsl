package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestSynchrotron1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.01, 0.444972504114210632, TEST_TOL0},
		{1.0, 0.651422815355364504, TEST_TOL1},
		{10.0, 0.000192238264300868882, TEST_TOL1},
		{100.0, 4.69759366592220221e-43, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Synchrotron_1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Synchrotron_1_e(%v)", c.x))
		})
	}
}

func TestSynchrotron2(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{0.01, 0.23098077342226277732, TEST_TOL2},
		{1.0, 0.4944750621042082670, TEST_TOL1},
		{10.0, 0.00018161187569530204281, TEST_TOL1},
		{256.0, 1.3272635474353774058e-110, TEST_TOL4}, /* exp()... not my fault */
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Synchrotron_2_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Synchrotron_2_e(%v)", c.x))
		})
	}
}
