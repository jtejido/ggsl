package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestHyperg0f1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		c, x, expected float64
		tol            float64
	}{
		{1, 0.5, 1.5660829297563505373, TEST_TOL0},
		{5, 0.5, 1.1042674404828684574, TEST_TOL1},
		{100, 30, 1.3492598639485110176, TEST_TOL2},
		{-0.5, 3, -39.29137997543434276, TEST_TOL1},
		{-100.5, 50, 0.6087930289227538496, TEST_TOL3},
		{1, -5.0, -0.3268752818235339109, TEST_TOL0},
		{-0.5, -5.0, -4.581634759005381184, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Hyperg_0F1_e(c.c, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Hyperg_0F1_e(%v,%v)", c.c, c.x))
		})
	}
}
