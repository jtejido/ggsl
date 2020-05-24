package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestDebye1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{0.1, 0.975277750004723276, TEST_TOL0},
		{1.0, 0.777504634112248239, TEST_TOL0},
		{10.0, 0.164443465679946027, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Debye_1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Debye_1_e(%v)", c.x))
		})
	}
}

func TestDebye2(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{0.1, 0.967083287045302664, TEST_TOL0},
		{1.0, 0.70787847562782924, TEST_TOL0},
		{10.0, 0.0479714980201218708, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Debye_2_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Debye_2_e(%v)", c.x))
		})
	}
}

func TestDebye3(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{0.1, 0.962999940487211048, TEST_TOL0},
		{1.0, 0.674415564077814667, TEST_TOL0},
		{10.0, 0.0192957656903454886, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Debye_3_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Debye_3_e(%v)", c.x))
		})
	}
}

func TestDebye4(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{0.1, 0.960555486124335944, TEST_TOL0},
		{1.0, 0.654874068886737049, TEST_TOL0},
		{10.0, 0.00967367556027115896, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Debye_4_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Debye_4_e(%v)", c.x))
		})
	}
}

func TestDebye5(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{0.1, 0.95892849428310568745, TEST_TOL0},
		{1.0, 0.6421002580217790246, TEST_TOL0},
		{10.0, 0.005701535852992908538, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Debye_5_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Debye_5_e(%v)", c.x))
		})
	}
}

func TestDebye6(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, tol float64
	}{
		{0.1, 0.95776777382605465878, TEST_TOL0},
		{1.0, 0.63311142583495107588, TEST_TOL0},
		{10.0, 3.7938493294615955279e-3, TEST_TOL0},
	}
	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Debye_6_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Debye_6_e(%v)", c.x))
		})
	}
}
