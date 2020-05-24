package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestDawson(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{1.0e-15, 1.0e-15, TEST_TOL0},
		{0.5, 0.4244363835020222959, TEST_TOL0},
		{2.0, 0.30134038892379196603, TEST_TOL0},
		{1000.0, 0.0005000002500003750009, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Dawson_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Dawson_e(%v)", c.x))
		})
	}
}
