package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestTransport2(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0e-10, 9.9999999999999999999e-11, TEST_TOL0},
		{1.0, 0.97303256135517012845, TEST_TOL0},
		{3.0, 2.41105004901695346199, TEST_TOL0},
		{10.0, 3.28432911449795173575, TEST_TOL0},
		{100.0, 3.28986813369645287294, TEST_TOL0},
		{1.0e+05, 3.28986813369645287294, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Transport_2_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Transport_2_e(%v)", c.x))
		})
	}
}

func TestTransport3(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0e-10, 4.999999999999999999997e-21, TEST_TOL0},
		{1.0, 0.479841006572417499939, TEST_TOL0},
		{3.0, 3.210604662942246772338, TEST_TOL0},
		{5.0, 5.614386613842273228585, TEST_TOL0},
		{10.0, 7.150322712008592975030, TEST_TOL0},
		{30.0, 7.212341416160946511930, TEST_TOL0},
		{100.0, 7.212341418957565712398, TEST_TOL0},
		{1.0e+05, 7.212341418957565712398, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Transport_3_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Transport_3_e(%v)", c.x))
		})
	}
}

func TestTransport4(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0e-10, 3.33333333333333333333e-31, TEST_TOL0},
		{1.0e-07, 3.33333333333333166666e-22, TEST_TOL0},
		{1.0e-04, 3.33333333166666666726e-13, TEST_TOL0},
		{0.1, 0.000333166726172109903824, TEST_TOL0},
		{1.0, 0.31724404523442648241, TEST_TOL0},
		{3.0, 5.96482239737147652446, TEST_TOL0},
		{5.0, 15.3597843168821829816, TEST_TOL0},
		{10.0, 25.2736676770304417334, TEST_TOL0},
		{30.0, 25.9757575220840937469, TEST_TOL0},
		{100.0, 25.9757576090673165963, TEST_TOL1},
		{1.0e+05, 25.9757576090673165963, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Transport_4_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Transport_4_e(%v)", c.x))
		})
	}
}

func TestTransport5(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0e-10, 2.49999999999999999999e-41, TEST_TOL0},
		{1.0e-07, 2.49999999999999861111e-29, TEST_TOL0},
		{1.0e-04, 2.49999999861111111163e-17, TEST_TOL0},
		{0.1, 0.000024986116317791487410, TEST_TOL0},
		{1.0, 0.236615879239094789259153, TEST_TOL0},
		{3.0, 12.77055769104415951115760, TEST_TOL0},
		{5.0, 50.26309221817518778543615, TEST_TOL0},
		{10.0, 116.3807454024207107698556, TEST_TOL0},
		{30.0, 124.4313279083858954839911, TEST_TOL0},
		{100.0, 124.4313306172043911597639, TEST_TOL0},
		{1.0e+05, 124.43133061720439115976, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Transport_5_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Transport_5_e(%v)", c.x))
		})
	}
}
