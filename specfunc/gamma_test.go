package specfunc

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/err"
	"strconv"
	"testing"
)

func TestLnGamma(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{-0.1, 2.368961332728788655, TEST_TOL0},
		{-1.0 / 256.0, 5.547444766967471595, TEST_TOL0},
		{1.0e-08, 18.420680738180208905, TEST_TOL0},
		{0.1, 2.252712651734205, TEST_TOL0},
		{1.0 + 1.0/256.0, -0.0022422226599611501448, TEST_TOL0},
		{2.0 + 1.0/256.0, 0.0016564177556961728692, TEST_TOL0},
		{100.0, 359.1342053695753, TEST_TOL0},
		{-1.0 - 1.0/65536.0, 11.090348438090047844, TEST_TOL0},
		{-1.0 - 1.0/268435456.0, 19.408121054103474300, TEST_TOL0},
		{-100.5, -364.9009683094273518, TEST_TOL0},
		{-100 - 1.0/65536.0, -352.6490910117097874, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Lngamma_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Lngamma_e(%v)", c.x))
		})
	}
}

func TestLnGammaSgn(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected, expSign float64
		tol                  float64
	}{
		{0.7, 0.26086724653166651439, 1., TEST_TOL1},
		{0.1, 2.2527126517342059599, 1., TEST_TOL0},
		{-0.1, 2.368961332728788655, -1., TEST_TOL0},
		{-1.0 - 1.0/65536.0, 11.090348438090047844, 1., TEST_TOL0},
		{-2.0 - 1.0/256.0, 4.848447725860607213, -1., TEST_TOL0},
		{-2.0 - 1.0/65536.0, 10.397193628164674967, -1., TEST_TOL0},
		{-3.0 - 1.0/8.0, 0.15431112768404182427, 1., TEST_TOL2},
		{-100.5, -364.9009683094273518, -1., TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			var sgn float64
			stat := Lngamma_sgn_e(c.x, r, &sgn)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Lngamma_sgn_e(%v)", c.x))

		})
	}
}

func TestGamma(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0 + 1.0/4096.0, 0.9998591371459403421, TEST_TOL0},
		{1.0 + 1.0/32.0, 0.9829010992836269148, TEST_TOL0},
		{2.0 + 1.0/256.0, 1.0016577903733583299, TEST_TOL0},
		{9.0, 40320.0, TEST_TOL0},
		{10.0, 362880.0, TEST_TOL0},
		{100.0, 9.332621544394415268e+155, TEST_TOL2},
		{170.0, 4.269068009004705275e+304, TEST_TOL2},
		{171.0, 7.257415615307998967e+306, TEST_TOL2},
		{-10.5, -2.640121820547716316e-07, TEST_TOL0},
		{-11.25, 6.027393816261931672e-08, TEST_TOL0},
		{-1.0 + 1.0/65536.0, -65536.42280587818970, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gamma_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gamma_e(%v)", c.x))
		})
	}
}

func TestGammaStar(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0e-08, 3989.423555759890865, TEST_TOL1},
		{1.0e-05, 126.17168469882690233, TEST_TOL0},
		{0.001, 12.708492464364073506, TEST_TOL0},
		{1.5, 1.0563442442685598666, TEST_TOL0},
		{3.0, 1.0280645179187893045, TEST_TOL0},
		{9.0, 1.0092984264218189715, TEST_TOL0},
		{11.0, 1.0076024283104962850, TEST_TOL0},
		{100.0, 1.0008336778720121418, TEST_TOL0},
		{1.0e+05, 1.0000008333336805529, TEST_TOL0},
		{1.0e+20, 1.0, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gammastar_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gammastar_e(%v)", c.x))
		})
	}
}

func TestGammaInv(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, expected float64
		tol         float64
	}{
		{1.0, 1.0, TEST_TOL0},
		{2.0, 1.0, TEST_TOL0},
		{3.0, 0.5, TEST_TOL0},
		{4.0, 1.0 / 6.0, TEST_TOL0},
		{10.0, 1.0 / 362880.0, TEST_TOL0},
		{100.0, 1.0715102881254669232e-156, TEST_TOL2},
		{0.0, 0.0, TEST_TOL0},
		{-1.0, 0.0, TEST_TOL0},
		{-2.0, 0.0, TEST_TOL0},
		{-3.0, 0.0, TEST_TOL0},
		{-4.0, 0.0, TEST_TOL0},
		{-10.5, -1.0 / 2.640121820547716316e-07, TEST_TOL2},
		{-11.25, 1.0 / 6.027393816261931672e-08, TEST_TOL1},
		{-1.0 + 1.0/65536.0, -1.0 / 65536.42280587818970, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gammainv_e(c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gammainv_e(%v)", c.x))
		})
	}
}

func TestTaylorCoeff(t *testing.T) {
	r := new(Result)
	cases := []struct {
		n           int
		x, expected float64
		tol         float64
	}{
		{10, 1.0 / 1048576.0, 1.7148961854776073928e-67, TEST_TOL0},
		{10, 1.0 / 1024.0, 2.1738891788497900281e-37, TEST_TOL0},
		{10, 1.0, 2.7557319223985890653e-07, TEST_TOL0},
		{10, 5.0, 2.6911444554673721340, TEST_TOL0},
		{10, 500.0, 2.6911444554673721340e+20, TEST_TOL0},
		{100, 100.0, 1.0715102881254669232e+42, TEST_TOL1},
		{1000, 200.0, 2.6628790558154746898e-267, TEST_TOL1},
		{1000, 500.0, 2.3193170139740855074e+131, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Taylorcoeff_e(c.n, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Taylorcoeff_e(%v, %v)", c.n, c.x))
		})
	}
}

func TestGammaInc(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, x, expected float64
		tol            float64
	}{
		{-1.0 / 1048576.0, 1.0 / 1048576.0, 13.285819596290624271, TEST_TOL0},
		{-0.001, 1.0 / 1048576.0, 13.381275128625328858, TEST_TOL0},
		{-1.0, 1.0 / 1048576.0, 1.0485617142715768655e+06, TEST_TOL0},
		{-0.00001, 0.001, 6.3317681434563592142, TEST_TOL0},
		{-0.0001, 0.001, 6.3338276439767189385, TEST_TOL0},
		{-0.001, 0.001, 6.3544709102510843793, TEST_TOL0},
		{-0.5, 0.001, 59.763880515942196981, TEST_TOL0},
		{-1.0, 0.001, 992.66896046923884234, TEST_TOL0},
		{-3.5, 0.001, 9.0224404490639003706e+09, TEST_TOL1},
		{-10.5, 0.001, 3.0083661558184815656e+30, TEST_TOL2},
		{-0.001, 0.1, 1.8249109609418620068, TEST_TOL0},
		{-0.5, 0.1, 3.4017693366916154163, TEST_TOL0},
		{-10.0, 0.1, 8.9490757483586989181e+08, TEST_TOL1},
		{-10.5, 0.1, 2.6967403834226421766e+09, TEST_TOL1},
		{-0.001, 1.0, 0.21928612679072766340, TEST_TOL1},
		{-0.5, 1.0, 0.17814771178156069019, TEST_TOL1},
		{-1.0, 1.0, 0.14849550677592204792, TEST_TOL1},
		{-2.5, 1.0, 0.096556648631275160264, TEST_TOL1},
		{-1.0, 10.0, 3.8302404656316087616e-07, TEST_TOL1},
		{-0.001, 10.0, 4.1470562324807320961e-06, TEST_TOL1},
		{-0.5, 10.0, 1.2609042613241570681e-06, TEST_TOL0},
		{-1.0, 10.0, 3.8302404656316087616e-07, TEST_TOL1},
		{-10.5, 10.0, 6.8404927328441566785e-17, TEST_TOL1},
		{-100.0, 10.0, 4.1238327669858313997e-107, TEST_TOL2},
		{-200.0, 10.0, 2.1614091830529343423e-207, TEST_TOL2},
		{0.0, 0.001, 6.3315393641361493320, TEST_TOL0},
		{0.001, 0.001, 6.3087159394864007261, TEST_TOL0},
		{1.0, 0.001, 0.99900049983337499167, TEST_TOL0},
		{10.0, 0.001, 362880.0, TEST_TOL0},
		{0.0, 1.0, 0.21938393439552027368, TEST_TOL0},
		{0.001, 1.0, 0.21948181320730279613, TEST_TOL1},
		{1.0, 1.0, 0.36787944117144232160, TEST_TOL0},
		{10.0, 1.0, 362879.95956592242045, TEST_TOL0},
		{100.0, 1.0, 9.3326215443944152682e+155, TEST_TOL0},
		{0.0, 100.0, 3.6835977616820321802e-46, TEST_TOL2},
		{0.001, 100.0, 3.7006367674063550631e-46, TEST_TOL2},
		{1.0, 100.0, 3.7200759760208359630e-44, TEST_TOL2},
		{10.0, 100.0, 4.0836606309106112723e-26, TEST_TOL2},
		{100.0, 100.0, 4.5421981208626694294e+155, TEST_TOL1},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gamma_inc_e(c.a, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gamma_inc_e(%v, %v)", c.a, c.x))
		})
	}
}

func TestGammaIncP(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, x, expected float64
		tol            float64
	}{
		{1e-100, 0.001, 1.0, TEST_TOL0},
		{0.001, 0.001, 0.9936876467088602902, TEST_TOL0},
		{0.001, 1.0, 0.9997803916424144436, TEST_TOL0},
		{0.001, 10.0, 0.9999999958306921828, TEST_TOL0},
		{1.0, 0.001, 0.0009995001666250083319, TEST_TOL0},
		{1.0, 1.01, 0.6357810204284766802, TEST_TOL0},
		{1.0, 10.0, 0.9999546000702375151, TEST_TOL0},
		{10.0, 10.01, 0.5433207586693410570, TEST_TOL0},
		{10.0, 20.0, 0.9950045876916924128, TEST_TOL0},
		{1000.0, 1000.1, 0.5054666401440661753, TEST_TOL2},
		{1000.0, 2000.0, 1.0, TEST_TOL0},
		{34.0, 32.0, 0.3849626436463866776322932129, TEST_TOL2},
		{37.0, 3.499999999999999289e+01, 0.3898035054195570860969333039, TEST_TOL2},
		{10, 1e-16, 2.755731922398588814734648067e-167, TEST_TOL2},
		{1263131.0, 1261282.3637, 0.04994777516935182963821362168, TEST_TOL4},
		{1263131.0, 1263131.0, 0.500118321758657770672882362502514254, TEST_TOL4},
		{100, 99.0, 0.4733043303994607, TEST_TOL2},
		{200, 199.0, 0.4811585880878718, TEST_TOL2},
		{5670, 4574, 3.063972328743934e-55, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gamma_inc_P_e(c.a, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gamma_inc_P_e(%v, %v)", c.a, c.x))
		})
	}
}

func TestGammaIncQ(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, x, expected float64
		tol            float64
	}{
		{0.0, 0.001, 0.0, TEST_TOL0},
		{0.001, 0.001, 0.006312353291139709793, TEST_TOL0},
		{0.001, 1.0, 0.00021960835758555639171, TEST_TOL1},
		{0.001, 2.0, 0.00004897691783098147880, TEST_TOL2},
		{0.001, 5.0, 1.1509813397308608541e-06, TEST_TOL1},
		{1.0, 0.001, 0.9990004998333749917, TEST_TOL0},
		{1.0, 1.01, 0.3642189795715233198, TEST_TOL0},
		{1.0, 10.0, 0.00004539992976248485154, TEST_TOL0},
		{10.0, 10.01, 0.4566792413306589430, TEST_TOL0},
		{10.0, 100.0, 1.1253473960842733885e-31, TEST_TOL2},
		{1000.0, 1000.1, 0.4945333598559338247, TEST_TOL2},
		{1000.0, 2000.0, 6.847349459614753180e-136, TEST_TOL2},
		{100, 99.0, 0.5266956696005394, TEST_TOL2},
		{200, 199.0, 0.5188414119121281, TEST_TOL2},
		{5670, 4574, 1.0000000000000000, TEST_TOL2},
		{1.0e+06 - 1.0, 1.0e+06 - 2.0, 0.50026596175224547004, TEST_TOL3},
		{1.0e+06 + 2.0, 1.0e+06 + 1.0, 0.50026596135330304336, TEST_TOL2},
		{1.0e+06, 1.0e+06 - 2.0, 0.50066490399940144811, TEST_TOL2},
		{1.0e+07, 1.0e+07 - 2.0, 0.50021026104978614908, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Gamma_inc_Q_e(c.a, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Gamma_inc_Q_e(%v, %v)", c.a, c.x))
		})
	}
}

func TestPoch(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, x, expected float64
		tol            float64
	}{
		{5, 0.0, 1.0, TEST_TOL0},
		{7, 3, 504.0, TEST_TOL0},
		{5, 2, 30.0, TEST_TOL1},
		{5, 1.0 / 256.0, 1.0059023106151364982, TEST_TOL0},
		{-9.0, -4.0, 1.0 / 17160.0, TEST_TOL0},
		{-9.0, -3.0, -1.0 / 1320.0, TEST_TOL0},
		{-9.0, -3.5, 0, TEST_TOL0},
		{-9.0, 4.0, 3024.0, TEST_TOL0},
		{-9.0, 3.0, -504.0, TEST_TOL0},
		{-9.0, 3.5, 0.0, TEST_TOL0},
		{-9.0, 0.0, 1.0, TEST_TOL0},
		{-8.0, -4.0, 1.0 / 11880.0, TEST_TOL0},
		{-8.0, -3.0, -1.0 / 990.0, TEST_TOL0},
		{-8.0, +4.0, 1680.0, TEST_TOL0},
		{-8.0, +3.0, -336.0, TEST_TOL0},
		{-3.0, +4.0, 0.0, TEST_TOL0},
		{-3.0, +3.0, -6.0, TEST_TOL2},
		{-4.0, +4.0, 24.0, TEST_TOL2},
		{-3.0, +100.0, 0.0, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Poch_e(c.a, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Poch_e(%v, %v)", c.a, c.x))
		})
	}
}

func TestPochrel(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, x, expected float64
		tol            float64
	}{
		{5, 0, 1.506117668431800472, TEST_TOL1},
		{7, 3, 503.0 / 3.0, TEST_TOL0},
		{5, 2, 29.0 / 2.0, TEST_TOL1},
		{5, 0.01, 1.5186393661368275330, TEST_TOL2},
		{-5.5, 0.01, 1.8584945633829063516, TEST_TOL1},
		{-5.5, -1.0 / 8.0, 1.0883319303552135488, TEST_TOL1},
		{-5.5, -1.0 / 256.0, 1.7678268037726177453, TEST_TOL1},
		{-5.5, -11.0, 0.09090909090939652475, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Pochrel_e(c.a, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Pochrel_e(%v, %v)", c.a, c.x))
		})
	}
}

func TestLnPoch(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, x, expected float64
		tol            float64
	}{
		{5, 0.0, 0.0, TEST_TOL0},
		{5, 1.0 / 65536.0, 0.000022981557571259389129, TEST_TOL0},
		{5, 1.0 / 256.0, 0.005884960217985189004, TEST_TOL2},
		{7, 3, 6.222576268071368616, TEST_TOL0},
		{5, 2, 3.401197381662155375, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Lnpoch_e(c.a, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Lnpoch_e(%v, %v)", c.a, c.x))
		})
	}
}

func TestLnBeta(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, y, expected float64
		tol            float64
	}{
		{1.0e-8, 1.0e-8, 19.113827924512310617, TEST_TOL0},
		{1.0e-8, 0.01, 18.420681743788563403, TEST_TOL0},
		{1.0e-8, 1.0, 18.420680743952365472, TEST_TOL0},
		{1.0e-8, 10.0, 18.420680715662683009, TEST_TOL0},
		{1.0e-8, 1000.0, 18.420680669107656949, TEST_TOL0},
		{0.1, 0.1, 2.9813614810376273949, TEST_TOL1},
		{0.1, 1.0, 2.3025850929940456840, TEST_TOL1},
		{0.1, 100.0, 1.7926462324527931217, TEST_TOL0},
		{0.1, 1000, 1.5619821298353164928, TEST_TOL0},
		{1.0, 1.00025, -0.0002499687552073570, TEST_TOL4},
		{1.0, 1.01, -0.009950330853168082848, TEST_TOL3},
		{1.0, 1000.0, -6.907755278982137052, TEST_TOL0},
		{100.0, 100.0, -139.66525908670663927, TEST_TOL2},
		{100.0, 1000.0, -336.4348576477366051, TEST_TOL0},
		{100.0, 1.0e+8, -1482.9339185256447309, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Lnbeta_e(c.x, c.y, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Lnbeta_e(%v, %v)", c.x, c.y))
		})
	}
}

func TestBeta(t *testing.T) {
	r := new(Result)
	cases := []struct {
		x, y, expected float64
		tol            float64
	}{
		{1.0, 1.0, 1.0, TEST_TOL0},
		{1.0, 1.001, 0.9990009990009990010, TEST_TOL0},
		{1.0, 5.0, 0.2, TEST_TOL1},
		{1.0, 100.0, 0.01, TEST_TOL1},
		{10.0, 100.0, 2.3455339739604649879e-15, TEST_TOL2},

		/* Test negative arguments */
		{2.5, -0.1, -11.43621278354402041480, TEST_TOL2},
		{2.5, -1.1, 14.555179906328753255202, TEST_TOL2},
		{-0.25, -0.1, -13.238937960945229110, TEST_TOL2},
		{-1.25, -0.1, -14.298052997820847439, TEST_TOL2},
		{-100.1, -99.1, -1.005181917797644630375787297e60, TEST_TOL3},
		{-100.1, 99.3, 0.0004474258199579694011200969001, TEST_TOL2},
		{100.1, -99.3, 1.328660939628876472028853747, TEST_TOL2},
		{-100.1, 1.2, 0.00365530364287960795444856281, TEST_TOL3},
		{100.1, -1.2, 1203.895236907821059270698160, TEST_TOL2},
		{-100.1, -1.2, -3236.073671884748847700283841, TEST_TOL2},
		{-100.001, 0.0099, -853.946649365611147996495177, TEST_TOL4},

		/* Other test cases */
		{1e-32, 1.5, 1e32, TEST_TOL2},
		{1e-6, 0.5, 1000001.386293677092419390336, TEST_TOL2},

		{-1.5, 0.5, 0.0, TEST_TOL0},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Beta_e(c.x, c.y, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Beta_e(%v, %v)", c.x, c.y))
		})
	}
}

func TestBetaInc(t *testing.T) {
	r := new(Result)
	cases := []struct {
		a, b, x, expected float64
		tol               float64
	}{
		{1.0, 1.0, 0.0, 0.0, TEST_TOL2},
		{1.0, 1.0, 1.0, 1.0, TEST_TOL2},
		{0.1, 0.1, 1.0, 1.0, TEST_TOL2},
		{1.0, 1.0, 0.5, 0.5, TEST_TOL2},
		{0.1, 1.0, 0.5, 0.9330329915368074160, TEST_TOL2},
		{10.0, 1.0, 0.5, 0.0009765625000000000000, TEST_TOL2},
		{50.0, 1.0, 0.5, 8.881784197001252323e-16, TEST_TOL2},
		{1.0, 0.1, 0.5, 0.06696700846319258402, TEST_TOL2},
		{1.0, 10.0, 0.5, 0.99902343750000000000, TEST_TOL2},
		{1.0, 50.0, 0.5, 0.99999999999999911180, TEST_TOL2},
		{1.0, 1.0, 0.1, 0.10, TEST_TOL2},
		{1.0, 2.0, 0.1, 0.19, TEST_TOL2},
		{1.0, 2.0, 0.9, 0.99, TEST_TOL2},
		{50.0, 60.0, 0.5, 0.8309072939016694143, TEST_TOL2},
		{90.0, 90.0, 0.5, 0.5, TEST_TOL2},
		{500.0, 500.0, 0.6, 0.9999999999157549630, TEST_TOL2},
		{5000.0, 5000.0, 0.4, 4.518543727260666383e-91, TEST_TOL5},
		{5000.0, 5000.0, 0.6, 1.0, TEST_TOL2},
		{5000.0, 2000.0, 0.6, 8.445388773903332659e-89, TEST_TOL5},

		{-0.1, -0.1, 1.0, 1.0, TEST_TOL2},
		{-0.1, -0.2, 1.0, 1.0, TEST_TOL2},
		{-0.2, -0.1, 1.0, 1.0, TEST_TOL2},

		{-0.1, -0.2, 0.5, 0.675252001958389971991335, TEST_TOL2},
		{-0.2, -0.1, 0.5, 0.324747998041610028008665, TEST_TOL2},

		{-0.1, -0.1, 0.0, 0.0, TEST_TOL2},
		{-0.1, -0.2, 0.0, 0.0, TEST_TOL2},
		{-0.2, -0.1, 0.0, 0.0, TEST_TOL2},

		{-0.1, -0.2, 0.3, 0.7469186777964287252, TEST_TOL2},
		{-0.2, -0.1, 0.3, 0.3995299653262016818, TEST_TOL2},

		/* Bug report from Thomas Tanner <tanner@gmx.de> */

		{0.5, 101.5, 0.999457, 1.0, TEST_TOL2},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Beta_inc_e(c.a, c.b, c.x, r)
			run_test_sf(t, stat, r, c.expected, c.tol, err.SUCCESS, fmt.Sprintf("Beta_inc_e(%v, %v, %v)", c.a, c.b, c.x))
		})
	}
}
