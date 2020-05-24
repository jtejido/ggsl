package specfunc

import (
	"fmt"
	"github.com/jtejido/ggsl/err"
	"strconv"
	"testing"
)

func TestExpintE1(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-1.0, -1.8951178163559367555, TEST_TOL0},
		{1.0e-10, 22.448635265138923980, TEST_TOL0},
		{1.0e-05, 10.935719800043695615, TEST_TOL0},
		{0.1, 1.82292395841939066610, TEST_TOL0},
		{1.0, 0.21938393439552027368, TEST_TOL0},
		{10.0, 4.156968929685324277e-06, TEST_TOL1},
		{50.0, 3.783264029550459019e-24, TEST_TOL2},
		{300.0, 1.710384276804510115e-133, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expint_E1_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expint_E1_e(%v)", c.x))
		})
	}
}

func TestExpintE2(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-1.0, 0.8231640121031084799, TEST_TOL1},
		{0.0, 1.0, TEST_TOL0},
		{1.0 / 4294967296.0, 0.9999999947372139168, TEST_TOL0},
		{1.0 / 65536.0, 0.9998243233207178845, TEST_TOL0},
		{0.1, 0.7225450221940205066, TEST_TOL0},
		{1.0, 0.14849550677592204792, TEST_TOL0},
		{10.0, 3.830240465631608762e-06, TEST_TOL1},
		{50.0, 3.711783318868827367e-24, TEST_TOL2},
		{300.0, 1.7047391998483433998e-133, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expint_E2_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expint_E2_e(%v)", c.x))
		})
	}
}

func TestExpintEn(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{1, -1.0, -1.8951178163559367555, TEST_TOL0},
		{1, 1.0e-10, 22.448635265138923980, TEST_TOL0},
		{1, 1.0e-05, 10.935719800043695615, TEST_TOL0},
		{1, 0.1, 1.82292395841939066610, TEST_TOL0},
		{1, 1.0, 0.21938393439552027368, TEST_TOL0},
		{1, 10.0, 4.156968929685324277e-06, TEST_TOL1},
		{1, 50.0, 3.783264029550459019e-24, TEST_TOL2},
		{1, 300.0, 1.710384276804510115e-133, TEST_TOL2},
		{2, -1.0, 0.8231640121031084799, TEST_TOL1},
		{2, 0.0, 1.0, TEST_TOL0},
		{2, 1.0 / 4294967296.0, 0.9999999947372139168, TEST_TOL0},
		{2, 1.0 / 65536.0, 0.9998243233207178845, TEST_TOL0},
		{2, 0.1, 0.7225450221940205066, TEST_TOL0},
		{2, 1.0, 0.14849550677592204792, TEST_TOL0},
		{2, 10.0, 3.830240465631608762e-06, TEST_TOL1},
		{2, 50.0, 3.711783318868827367e-24, TEST_TOL2},
		{2, 300.0, 1.7047391998483433998e-133, TEST_TOL2},
		{3, 0.0, 0.5, TEST_TOL0},
		{3, 1.0 / 4294967296.0, 0.499999999767169356972, TEST_TOL1},
		{3, 1.0 / 65536.0, 0.4999847426094515610, TEST_TOL0},
		{3, 0.1, 0.4162914579082787612543, TEST_TOL0},
		{3, 1.0, 0.10969196719776013683858, TEST_TOL1},
		{3, 10.0, 0.000003548762553084381959981, TEST_TOL1},
		{3, 50.0, 3.6429094264752049812e-24, TEST_TOL2},
		{3, 300.0, 1.699131143349179084e-133, TEST_TOL2},
		{10, 0.0, 0.111111111111111111, TEST_TOL0},
		{10, 1.0 / 4294967296.0, 0.111111111082007280658, TEST_TOL2},
		{10, 1.0 / 65536.0, 0.11110920377910896018606, TEST_TOL1},
		{10, 0.1, 0.099298432000896813567905, TEST_TOL1},
		{10, 1.0, 0.036393994031416401634164534, TEST_TOL1},
		{10, 10.0, 0.00000232530265702821081778968, TEST_TOL1},
		{10, 50.0, 3.223296586749110919572e-24, TEST_TOL2},
		{10, 300.0, 1.6608815083360041367294736e-133, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expint_En_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expint_En_e(%v, %v)", c.n, c.x))
		})
	}
}

func TestExpintEi(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-1.0, -0.21938393439552027368, TEST_TOL0},
		{1.0 / 4294967296.0, -21.603494112783886397, TEST_TOL0},
		{1.0, 1.8951178163559367555, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expint_Ei_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expint_Ei_e(%v)", c.x))
		})
	}
}

func TestExpintE1Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10000.0, -0.00010001000200060024012, TEST_TOL0},
		{-1000.0, -0.0010010020060241207251, TEST_TOL0},
		{-10.0, -0.11314702047341077803, TEST_TOL0},
		{-1.0, -0.69717488323506606877, TEST_TOL0},
		{1.0e-10, 22.448635267383787506, TEST_TOL0},
		{1.0e-05, 10.935829157788483865, TEST_TOL0},
		{0.1, 2.0146425447084516791, TEST_TOL0},
		{1.0, 0.59634736232319407434, TEST_TOL0},
		{10.0, 0.091563333939788081876, TEST_TOL0},
		{50.0, 0.019615109930114870365, TEST_TOL0},
		{300.0, 0.0033222955652707070644, TEST_TOL0},
		{1000.0, 0.00099900199402388071500, TEST_TOL0},
		{10000.0, 0.000099990001999400239880, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expint_E1_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expint_E1_scaled_e(%v)", c.x))
		})
	}
}

func TestExpintE2Scaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-10000.0, -0.00010002000600240120072, TEST_TOL3},
		{-1000.0, -0.0010020060241207250807, TEST_TOL3},
		{-10.0, -0.13147020473410778034, TEST_TOL1},
		{-1.0, 0.30282511676493393123, TEST_TOL1},
		{0.0, 1.0, TEST_TOL1},
		{1.0 / 4294967296.0, 0.99999999497004455927, TEST_TOL0},
		{1.0 / 65536.0, 0.99983957954556245453, TEST_TOL0},
		{0.1, 0.79853574552915483209, TEST_TOL0},
		{1.0, 0.40365263767680592566, TEST_TOL0},
		{10.0, 0.084366660602119181239, TEST_TOL1},
		{50.0, 0.019244503494256481735, TEST_TOL2},
		{300.0, 0.0033113304187878806691, TEST_TOL0},
		{1000.0, 0.00099800597611928500004, TEST_TOL0},
		{10000.0, 0.000099980005997601199281, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expint_E2_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expint_E2_scaled_e(%v)", c.x))
		})
	}
}

func TestExpintEnScaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{1, -10000.0, -0.00010001000200060024012, TEST_TOL0},
		{1, -1000.0, -0.0010010020060241207251, TEST_TOL0},
		{1, -10.0, -0.11314702047341077803, TEST_TOL0},
		{1, -1.0, -0.69717488323506606877, TEST_TOL0},
		{1, 1.0e-10, 22.448635267383787506, TEST_TOL0},
		{1, 1.0e-05, 10.935829157788483865, TEST_TOL0},
		{1, 0.1, 2.0146425447084516791, TEST_TOL0},
		{1, 1.0, 0.59634736232319407434, TEST_TOL0},
		{1, 10.0, 0.091563333939788081876, TEST_TOL0},
		{1, 50.0, 0.019615109930114870365, TEST_TOL0},
		{1, 300.0, 0.0033222955652707070644, TEST_TOL0},
		{1, 1000.0, 0.00099900199402388071500, TEST_TOL0},
		{1, 10000.0, 0.000099990001999400239880, TEST_TOL0},

		{2, -10000.0, -0.00010002000600240120072, TEST_TOL3},
		{2, -1000.0, -0.0010020060241207250807, TEST_TOL3},
		{2, -10.0, -0.13147020473410778034, TEST_TOL1},
		{2, -1.0, 0.30282511676493393123, TEST_TOL1},
		{2, 0.0, 1.0, TEST_TOL1},
		{2, 1.0 / 4294967296.0, 0.99999999497004455927, TEST_TOL0},
		{2, 1.0 / 65536.0, 0.99983957954556245453, TEST_TOL0},
		{2, 0.1, 0.79853574552915483209, TEST_TOL0},
		{2, 1.0, 0.40365263767680592566, TEST_TOL0},
		{2, 10.0, 0.084366660602119181239, TEST_TOL1},
		{2, 50.0, 0.019244503494256481735, TEST_TOL2},
		{2, 300.0, 0.0033113304187878806691, TEST_TOL0},
		{2, 1000.0, 0.00099800597611928500004, TEST_TOL0},
		{2, 10000.0, 0.000099980005997601199281, TEST_TOL0},

		{3, 0.0, 0.5, TEST_TOL0},
		{3, 1.0 / 4294967296.0, 0.4999999998835846787586, TEST_TOL1},
		{3, 1.0 / 65536.0, 0.4999923718293796877864492, TEST_TOL0},
		{3, 0.1, 0.4600732127235422583955, TEST_TOL0},
		{3, 1.0, 0.298173681161597037170539, TEST_TOL1},
		{3, 10.0, 0.07816669698940409380349, TEST_TOL1},
		{3, 50.0, 0.0188874126435879566345, TEST_TOL2},
		{3, 300.0, 0.00330043718181789963028657675, TEST_TOL2},

		{10, 0.0, 0.111111111111111111, TEST_TOL0},
		{10, 1.0 / 4294967296.0, 0.11111111110787735217158, TEST_TOL2},
		{10, 1.0 / 65536.0, 0.1111108991839472074435, TEST_TOL1},
		{10, 0.1, 0.1097417392579033988025, TEST_TOL1},
		{10, 1.0, 0.09892913264064615521915, TEST_TOL1},
		{10, 10.0, 0.0512181994376050593314159875, TEST_TOL1},
		{10, 50.0, 0.0167118436335939556034579, TEST_TOL2},
		{10, 300.0, 0.0032261400811599644878615, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expint_En_scaled_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expint_En_scaled_e(%v, %v)", c.n, c.x))
		})
	}
}

func TestExpintEiScaled(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-1000.0, -0.00099900199402388071500, TEST_TOL0},
		{-1.0, -0.59634736232319407434, TEST_TOL0},
		{1.0 / 4294967296.0, -21.603494107753930958, TEST_TOL0},
		{1.0, 0.69717488323506606877, TEST_TOL0},
		{1000.0, 0.0010010020060241207251, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expint_Ei_scaled_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expint_Ei_scaled_e(%v)", c.x))
		})
	}
}

func TestExpint3(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0e-10, 1.0e-10, TEST_TOL0},
		{1.0e-05, 9.9999999999999975e-06, TEST_TOL0},
		{0.1, 0.09997500714119079665122, TEST_TOL0},
		{0.5, 0.48491714311363971332427, TEST_TOL0},
		{1.0, 0.80751118213967145285833, TEST_TOL0},
		{2.0, 0.89295351429387631138208, TEST_TOL0},
		{5.0, 0.89297951156924921121856, TEST_TOL0},
		{10.0, 0.89297951156924921121856, TEST_TOL0},
		{100.0, 0.89297951156924921121856, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Expint_3_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Expint_3_e(%v)", c.x))
		})
	}
}
