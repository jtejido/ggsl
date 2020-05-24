package specfunc

import (
	"fmt"
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/err"
	"math"
	"strconv"
	"testing"
)

func TestEllintKcomp(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		k        float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		{0.99, mode, 3.3566005233611923760, TEST_TOL0, err.SUCCESS},
		{0.50, mode, 1.6857503548125960429, TEST_TOL0, err.SUCCESS},
		{0.010, mode, 1.5708355989121522360, TEST_TOL0, err.SUCCESS},

		/* Bug report from Thies Heidecke */
		{0.99999999906867742538, mode, 11.4369284843320018031, TEST_SNGL, err.SUCCESS},
		{0.99, mode, 3.3566005233611923760, TEST_TOL0, err.SUCCESS},
		{0.50, mode, 1.6857503548125960429, TEST_TOL0, err.SUCCESS},
		{0.010, mode, 1.5708355989121522360, TEST_TOL0, err.SUCCESS},

		/* Bug report from Thies Heidecke */
		{0.99999999906867742538, mode, 11.4369284843320018031, TEST_SNGL, err.SUCCESS},
		/* Bug report from Will M. Farr bug #31362 */
		/* FIXME: we are accepting MAXITER as the return code, but really
		   this should be changed to EINVAL in the routine itself */

		{math.NaN(), mode, math.NaN(), math.NaN(), err.EMAXITER},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_Kcomp_e(c.k, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_Kcomp_e(%v, %v)", c.k, c.mode))
		})
	}
}

func TestEllintEcomp(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		k        float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		{0.99, mode, 1.0284758090288040010, TEST_TOL0, err.SUCCESS},
		{0.50, mode, 1.4674622093394271555, TEST_TOL0, err.SUCCESS},
		{0.01, mode, 1.5707570561503852873, TEST_TOL0, err.SUCCESS},
		/* Bug report from Will M. Farr bug #31362 */
		/* FIXME: we are accepting MAXITER as the return code, but really
		   this should be changed to EINVAL in the routine itself */

		{math.NaN(), mode, math.NaN(), math.NaN(), err.EMAXITER},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_Ecomp_e(c.k, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_Ecomp_e(%v, %v)", c.k, c.mode))
		})
	}
}

func TestEllintPcomp(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		k, n     float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		{0.99, 0.1, mode, 3.13792612351836506315593, TEST_TOL0, err.SUCCESS},
		{0.50, 0.1, mode, 1.60455249360848890075108, TEST_TOL0, err.SUCCESS},
		{0.01, 0.1, mode, 1.49773208536003801277453, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_Pcomp_e(c.k, c.n, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_Pcomp_e(%v, %v, %v)", c.k, c.n, c.mode))
		})
	}
}

func TestEllintDcomp(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		k        float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		{0.99, mode, 2.375395076351788975665323192, TEST_TOL0, err.SUCCESS},
		{0.50, mode, 0.8731525818926755496456335628, TEST_TOL0, err.SUCCESS},
		{0.01, mode, 0.7854276176694868932799393751, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_Dcomp_e(c.k, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_Dcomp_e(%v, %v)", c.k, c.mode))
		})
	}
}

func TestEllintF(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		phi, k   float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		{gsl.Pi / 3.0, 0.99, mode, 1.3065333392738766762, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 3.0, 0.50, mode, 1.0895506700518854093, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 3.0, 0.01, mode, 1.0472129063770918952, TEST_TOL0, err.SUCCESS},
		/* F, argument phi > pi/2  */

		{gsl.Pi / 2.0, 0.99, mode, 3.35660052336119237603347, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 2.0, 0.50, mode, 1.68575035481259604287120, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 2.0, 0.01, mode, 1.57083559891215223602641, TEST_TOL0, err.SUCCESS},

		{2 * gsl.Pi / 3.0, 0.99, mode, 5.40666770744850807588478, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi / 3.0, 0.50, mode, 2.28195003957330667648585, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi / 3.0, 0.01, mode, 2.09445829144721257687207, TEST_TOL0, err.SUCCESS},

		{gsl.Pi, 0.99, mode, 6.71320104672238475206694, TEST_TOL0, err.SUCCESS},
		{gsl.Pi, 0.50, mode, 3.37150070962519208574241, TEST_TOL0, err.SUCCESS},
		{gsl.Pi, 0.01, mode, 3.14167119782430447205281, TEST_TOL0, err.SUCCESS},

		{4 * gsl.Pi / 3, 0.99, mode, 8.01973438599626142824910, TEST_TOL0, err.SUCCESS},
		{4 * gsl.Pi / 3, 0.50, mode, 4.46105137967707749499897, TEST_TOL0, err.SUCCESS},
		{4 * gsl.Pi / 3, 0.01, mode, 4.18888410420139636723356, TEST_TOL0, err.SUCCESS},

		{3 * gsl.Pi / 2.0, 0.99, mode, 10.0698015700835771281004, TEST_TOL0, err.SUCCESS},
		{3 * gsl.Pi / 2.0, 0.50, mode, 5.05725106443778812861361, TEST_TOL0, err.SUCCESS},
		{3 * gsl.Pi / 2.0, 0.01, mode, 4.71250679673645670807922, TEST_TOL0, err.SUCCESS},

		{5 * gsl.Pi / 3, 0.99, mode, 12.1198687541708928279517, TEST_TOL0, err.SUCCESS},
		{5 * gsl.Pi / 3, 0.50, mode, 5.65345074919849876222825, TEST_TOL0, err.SUCCESS},
		{5 * gsl.Pi / 3, 0.01, mode, 5.23612948927151704892488, TEST_TOL0, err.SUCCESS},

		{2 * gsl.Pi, 0.99, mode, 13.4264020934447695041339, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi, 0.50, mode, 6.74300141925038417148481, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi, 0.01, mode, 6.28334239564860894410562, TEST_TOL0, err.SUCCESS},

		{7 * gsl.Pi / 3.0, 0.99, mode, 14.7329354327186461803160, TEST_TOL0, err.SUCCESS},
		{7 * gsl.Pi / 3.0, 0.50, mode, 7.83255208930226958074138, TEST_TOL0, err.SUCCESS},
		{7 * gsl.Pi / 3.0, 0.01, mode, 7.33055530202570083928637, TEST_TOL0, err.SUCCESS},

		/* F, negative argument phi < 0  */

		{-gsl.Pi / 2.0, 0.99, mode, -3.35660052336119237603347, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi / 2.0, 0.50, mode, -1.68575035481259604287120, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi / 2.0, 0.01, mode, -1.57083559891215223602641, TEST_TOL0, err.SUCCESS},

		{-2 * gsl.Pi / 3.0, 0.99, mode, -5.40666770744850807588478, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi / 3.0, 0.50, mode, -2.28195003957330667648585, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi / 3.0, 0.01, mode, -2.09445829144721257687207, TEST_TOL0, err.SUCCESS},

		{-gsl.Pi, 0.99, mode, -6.71320104672238475206694, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi, 0.50, mode, -3.37150070962519208574241, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi, 0.01, mode, -3.14167119782430447205281, TEST_TOL0, err.SUCCESS},

		{-4 * gsl.Pi / 3, 0.99, mode, -8.01973438599626142824910, TEST_TOL0, err.SUCCESS},
		{-4 * gsl.Pi / 3, 0.50, mode, -4.46105137967707749499897, TEST_TOL0, err.SUCCESS},
		{-4 * gsl.Pi / 3, 0.01, mode, -4.18888410420139636723356, TEST_TOL0, err.SUCCESS},

		{-3 * gsl.Pi / 2.0, 0.99, mode, -10.0698015700835771281004, TEST_TOL0, err.SUCCESS},
		{-3 * gsl.Pi / 2.0, 0.50, mode, -5.05725106443778812861361, TEST_TOL0, err.SUCCESS},
		{-3 * gsl.Pi / 2.0, 0.01, mode, -4.71250679673645670807922, TEST_TOL0, err.SUCCESS},

		{-5 * gsl.Pi / 3, 0.99, mode, -12.1198687541708928279517, TEST_TOL0, err.SUCCESS},
		{-5 * gsl.Pi / 3, 0.50, mode, -5.65345074919849876222825, TEST_TOL0, err.SUCCESS},
		{-5 * gsl.Pi / 3, 0.01, mode, -5.23612948927151704892488, TEST_TOL0, err.SUCCESS},

		{-2 * gsl.Pi, 0.99, mode, -13.4264020934447695041339, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi, 0.50, mode, -6.74300141925038417148481, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi, 0.01, mode, -6.28334239564860894410562, TEST_TOL0, err.SUCCESS},

		{-7 * gsl.Pi / 3.0, 0.99, mode, -14.7329354327186461803160, TEST_TOL0, err.SUCCESS},
		{-7 * gsl.Pi / 3.0, 0.50, mode, -7.83255208930226958074138, TEST_TOL0, err.SUCCESS},
		{-7 * gsl.Pi / 3.0, 0.01, mode, -7.33055530202570083928637, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_F_e(c.phi, c.k, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_F_e(%v, %v, %v)", c.phi, c.k, c.mode))
		})
	}
}

func TestEllintE(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		phi, k   float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		{gsl.Pi / 3.0, 0.99, mode, 0.8704819220377943536, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 3.0, 0.50, mode, 1.0075555551444720293, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 3.0, 0.01, mode, 1.0471821963889481104, TEST_TOL0, err.SUCCESS},
		/* E, argument phi > pi/2  */

		{gsl.Pi / 2.0, 0.99, mode, 1.02847580902880400098389, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 2.0, 0.50, mode, 1.46746220933942715545980, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 2.0, 0.01, mode, 1.57075705615038528733708, TEST_TOL0, err.SUCCESS},

		{2 * gsl.Pi / 3.0, 0.99, mode, 1.18646969601981364833972, TEST_TOL1, err.SUCCESS},
		{2 * gsl.Pi / 3.0, 0.50, mode, 1.92736886353438228163734, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi / 3.0, 0.01, mode, 2.09433191591182246425715, TEST_TOL0, err.SUCCESS},

		{gsl.Pi, 0.99, mode, 2.05695161805760800196777, TEST_TOL0, err.SUCCESS},
		{gsl.Pi, 0.50, mode, 2.93492441867885431091959, TEST_TOL0, err.SUCCESS},
		{gsl.Pi, 0.01, mode, 3.14151411230077057467416, TEST_TOL0, err.SUCCESS},

		{4 * gsl.Pi / 3, 0.99, mode, 2.92743354009540235559582, TEST_TOL0, err.SUCCESS},
		{4 * gsl.Pi / 3, 0.50, mode, 3.94247997382332634020184, TEST_TOL0, err.SUCCESS},
		{4 * gsl.Pi / 3, 0.01, mode, 4.18869630868971868509117, TEST_TOL0, err.SUCCESS},

		{3 * gsl.Pi / 2.0, 0.99, mode, 3.08542742708641200295166, TEST_TOL0, err.SUCCESS},
		{3 * gsl.Pi / 2.0, 0.50, mode, 4.40238662801828146637939, TEST_TOL0, err.SUCCESS},
		{3 * gsl.Pi / 2.0, 0.01, mode, 4.71227116845115586201123, TEST_TOL0, err.SUCCESS},

		{5 * gsl.Pi / 3, 0.99, mode, 3.24342131407742165030750, TEST_TOL1, err.SUCCESS},
		{5 * gsl.Pi / 3, 0.50, mode, 4.86229328221323659255693, TEST_TOL0, err.SUCCESS},
		{5 * gsl.Pi / 3, 0.01, mode, 5.23584602821259303893130, TEST_TOL0, err.SUCCESS},

		{2 * gsl.Pi, 0.99, mode, 4.11390323611521600393555, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi, 0.50, mode, 5.86984883735770862183918, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi, 0.01, mode, 6.28302822460154114934831, TEST_TOL0, err.SUCCESS},

		{7 * gsl.Pi / 3.0, 0.99, mode, 4.98438515815301035756360, TEST_TOL0, err.SUCCESS},
		{7 * gsl.Pi / 3.0, 0.50, mode, 6.87740439250218065112143, TEST_TOL0, err.SUCCESS},
		{7 * gsl.Pi / 3.0, 0.01, mode, 7.33021042099048925976532, TEST_TOL0, err.SUCCESS},

		/* Test some negative arguments, phi < 0 */

		{-gsl.Pi / 2.0, 0.99, mode, -1.02847580902880400098389, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi / 2.0, 0.50, mode, -1.46746220933942715545980, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi / 2.0, 0.01, mode, -1.57075705615038528733708, TEST_TOL0, err.SUCCESS},

		{-2 * gsl.Pi / 3.0, 0.99, mode, -1.18646969601981364833972, TEST_TOL1, err.SUCCESS},
		{-2 * gsl.Pi / 3.0, 0.50, mode, -1.92736886353438228163734, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi / 3.0, 0.01, mode, -2.09433191591182246425715, TEST_TOL0, err.SUCCESS},

		{-gsl.Pi, 0.99, mode, -2.05695161805760800196777, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi, 0.50, mode, -2.93492441867885431091959, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi, 0.01, mode, -3.14151411230077057467416, TEST_TOL0, err.SUCCESS},

		{-4 * gsl.Pi / 3, 0.99, mode, -2.92743354009540235559582, TEST_TOL0, err.SUCCESS},
		{-4 * gsl.Pi / 3, 0.50, mode, -3.94247997382332634020184, TEST_TOL0, err.SUCCESS},
		{-4 * gsl.Pi / 3, 0.01, mode, -4.18869630868971868509117, TEST_TOL0, err.SUCCESS},

		{-3 * gsl.Pi / 2.0, 0.99, mode, -3.08542742708641200295166, TEST_TOL0, err.SUCCESS},
		{-3 * gsl.Pi / 2.0, 0.50, mode, -4.40238662801828146637939, TEST_TOL0, err.SUCCESS},
		{-3 * gsl.Pi / 2.0, 0.01, mode, -4.71227116845115586201123, TEST_TOL0, err.SUCCESS},

		{-5 * gsl.Pi / 3, 0.99, mode, -3.24342131407742165030750, TEST_TOL1, err.SUCCESS},
		{-5 * gsl.Pi / 3, 0.50, mode, -4.86229328221323659255693, TEST_TOL0, err.SUCCESS},
		{-5 * gsl.Pi / 3, 0.01, mode, -5.23584602821259303893130, TEST_TOL0, err.SUCCESS},

		{-2 * gsl.Pi, 0.99, mode, -4.11390323611521600393555, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi, 0.50, mode, -5.86984883735770862183918, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi, 0.01, mode, -6.28302822460154114934831, TEST_TOL0, err.SUCCESS},

		{-7 * gsl.Pi / 3.0, 0.99, mode, -4.98438515815301035756360, TEST_TOL0, err.SUCCESS},
		{-7 * gsl.Pi / 3.0, 0.50, mode, -6.87740439250218065112143, TEST_TOL0, err.SUCCESS},
		{-7 * gsl.Pi / 3.0, 0.01, mode, -7.33021042099048925976532, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_E_e(c.phi, c.k, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_E_e(%v, %v, %v)", c.phi, c.k, c.mode))
		})
	}
}

func TestEllintP(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		phi, k, n float64
		mode      gsl.MODE_T
		expected  float64
		tol       float64
		status    int
	}{
		{gsl.Pi / 3.0, 0.99, 0.5, mode, 1.1288726598764099882, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 3.0, 0.50, 0.5, mode, 0.9570574331323584890, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 3.0, 0.01, 0.5, mode, 0.9228868127118118465, TEST_TOL0, err.SUCCESS},
		/* P, argument phi > pi/2  */

		{gsl.Pi / 2.0, 0.99, -0.1, mode, 3.61678162163246646783050, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 2.0, 0.50, -0.1, mode, 1.78030349465454812629168, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 2.0, 0.01, -0.1, mode, 1.65580719756898353270922, TEST_TOL0, err.SUCCESS},

		{2 * gsl.Pi / 3.0, 0.99, -0.1, mode, 5.88008918207571119911983, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi / 3.0, 0.50, -0.1, mode, 2.43655207300356008717867, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi / 3.0, 0.01, -0.1, mode, 2.23211110528200554950903, TEST_TOL0, err.SUCCESS},

		{gsl.Pi, 0.99, -0.1, mode, 7.23356324326493293566099, TEST_TOL0, err.SUCCESS},
		{gsl.Pi, 0.50, -0.1, mode, 3.56060698930909625258336, TEST_TOL0, err.SUCCESS},
		{gsl.Pi, 0.01, -0.1, mode, 3.31161439513796706541844, TEST_TOL0, err.SUCCESS},

		{4 * gsl.Pi / 3, 0.99, -0.1, mode, 8.58703730445415467220216, TEST_TOL0, err.SUCCESS},
		{4 * gsl.Pi / 3, 0.50, -0.1, mode, 4.68466190561463241798805, TEST_TOL0, err.SUCCESS},
		{4 * gsl.Pi / 3, 0.01, -0.1, mode, 4.39111768499392858132786, TEST_TOL0, err.SUCCESS},

		{3 * gsl.Pi / 2.0, 0.99, -0.1, mode, 10.8503448648973994034915, TEST_TOL0, err.SUCCESS},
		{3 * gsl.Pi / 2.0, 0.50, -0.1, mode, 5.34091048396364437887504, TEST_TOL0, err.SUCCESS},
		{3 * gsl.Pi / 2.0, 0.01, -0.1, mode, 4.96742159270695059812767, TEST_TOL0, err.SUCCESS},

		{5 * gsl.Pi / 3, 0.99, -0.1, mode, 13.1136524253406441347808, TEST_TOL0, err.SUCCESS},
		{5 * gsl.Pi / 3, 0.50, -0.1, mode, 5.99715906231265633976204, TEST_TOL0, err.SUCCESS},
		{5 * gsl.Pi / 3, 0.01, -0.1, mode, 5.54372550041997261492747, TEST_TOL0, err.SUCCESS},

		{2 * gsl.Pi, 0.99, -0.1, mode, 14.4671264865298658713220, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi, 0.50, -0.1, mode, 7.12121397861819250516672, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi, 0.01, -0.1, mode, 6.62322879027593413083689, TEST_TOL0, err.SUCCESS},

		{7 * gsl.Pi / 3.0, 0.99, -0.1, mode, 15.8206005477190876078631, TEST_TOL0, err.SUCCESS},
		{7 * gsl.Pi / 3.0, 0.50, -0.1, mode, 8.24526889492372867057141, TEST_TOL0, err.SUCCESS},
		{7 * gsl.Pi / 3.0, 0.01, -0.1, mode, 7.70273208013189564674630, TEST_TOL0, err.SUCCESS},

		/* P, negative argument phi < 0  */

		{-gsl.Pi / 2.0, 0.99, -0.1, mode, -3.61678162163246646783050, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi / 2.0, 0.50, -0.1, mode, -1.78030349465454812629168, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi / 2.0, 0.01, -0.1, mode, -1.65580719756898353270922, TEST_TOL0, err.SUCCESS},

		{-2 * gsl.Pi / 3.0, 0.99, -0.1, mode, -5.88008918207571119911983, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi / 3.0, 0.50, -0.1, mode, -2.43655207300356008717867, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi / 3.0, 0.01, -0.1, mode, -2.23211110528200554950903, TEST_TOL0, err.SUCCESS},

		{-gsl.Pi, 0.99, -0.1, mode, -7.23356324326493293566099, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi, 0.50, -0.1, mode, -3.56060698930909625258336, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi, 0.01, -0.1, mode, -3.31161439513796706541844, TEST_TOL0, err.SUCCESS},

		{-4 * gsl.Pi / 3, 0.99, -0.1, mode, -8.58703730445415467220216, TEST_TOL0, err.SUCCESS},
		{-4 * gsl.Pi / 3, 0.50, -0.1, mode, -4.68466190561463241798805, TEST_TOL0, err.SUCCESS},
		{-4 * gsl.Pi / 3, 0.01, -0.1, mode, -4.39111768499392858132786, TEST_TOL0, err.SUCCESS},

		{-3 * gsl.Pi / 2.0, 0.99, -0.1, mode, -10.8503448648973994034915, TEST_TOL0, err.SUCCESS},
		{-3 * gsl.Pi / 2.0, 0.50, -0.1, mode, -5.34091048396364437887504, TEST_TOL0, err.SUCCESS},
		{-3 * gsl.Pi / 2.0, 0.01, -0.1, mode, -4.96742159270695059812767, TEST_TOL0, err.SUCCESS},

		{-5 * gsl.Pi / 3, 0.99, -0.1, mode, -13.1136524253406441347808, TEST_TOL0, err.SUCCESS},
		{-5 * gsl.Pi / 3, 0.50, -0.1, mode, -5.99715906231265633976204, TEST_TOL0, err.SUCCESS},
		{-5 * gsl.Pi / 3, 0.01, -0.1, mode, -5.54372550041997261492747, TEST_TOL0, err.SUCCESS},

		{-2 * gsl.Pi, 0.99, -0.1, mode, -14.4671264865298658713220, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi, 0.50, -0.1, mode, -7.12121397861819250516672, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi, 0.01, -0.1, mode, -6.62322879027593413083689, TEST_TOL0, err.SUCCESS},

		{-7 * gsl.Pi / 3.0, 0.99, -0.1, mode, -15.8206005477190876078631, TEST_TOL0, err.SUCCESS},
		{-7 * gsl.Pi / 3.0, 0.50, -0.1, mode, -8.24526889492372867057141, TEST_TOL0, err.SUCCESS},
		{-7 * gsl.Pi / 3.0, 0.01, -0.1, mode, -7.70273208013189564674630, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_P_e(c.phi, c.k, c.n, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_P_e(%v, %v, %v, %v)", c.phi, c.k, c.n, c.mode))
		})
	}
}

func TestEllintRF(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		x, y, z  float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		{5.0e-11, 1.0e-10, 1.0, mode, 12.36441982979439, TEST_TOL0, err.SUCCESS},
		{1.0, 2.0, 3.0, mode, 0.7269459354689082, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_RF_e(c.x, c.y, c.z, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_RF_e(%v, %v, %v, %v)", c.x, c.y, c.z, c.mode))
		})
	}
}

func TestEllintRD(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		x, y, z  float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		{5.0e-11, 1.0e-10, 1.0, mode, 34.0932594919337362, TEST_TOL0, err.SUCCESS},
		{1.0, 2.0, 3.0, mode, 0.2904602810289906, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_RD_e(c.x, c.y, c.z, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_RD_e(%v, %v, %v, %v)", c.x, c.y, c.z, c.mode))
		})
	}
}

func TestEllintRC(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		x, y     float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		{1.0, 2.0, mode, 0.7853981633974482, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_RC_e(c.x, c.y, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_RC_e(%v, %v, %v)", c.x, c.y, c.mode))
		})
	}
}

func TestEllintRJ(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		x, y, z, p float64
		mode       gsl.MODE_T
		expected   float64
		tol        float64
		status     int
	}{
		{2.0, 3.0, 4.0, 5.0, mode, 0.1429757966715675, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_RJ_e(c.x, c.y, c.z, c.p, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_RJ_e(%v, %v, %v, %v, %v)", c.x, c.y, c.z, c.p, c.mode))
		})
	}
}

func TestEllintD(t *testing.T) {
	r := new(Result)
	mode := gsl.MODE_DEFAULT
	cases := []struct {
		phi, k   float64
		mode     gsl.MODE_T
		expected float64
		tol      float64
		status   int
	}{
		/* D, argument phi > pi/2  */

		{gsl.Pi / 2.0, 0.99, mode, 2.375395076351788975665323192, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 2.0, 0.50, mode, 0.8731525818926755496456335628, TEST_TOL0, err.SUCCESS},
		{gsl.Pi / 2.0, 0.01, mode, 0.7854276176694868932799393751, TEST_TOL0, err.SUCCESS},

		{2 * gsl.Pi / 3.0, 0.99, mode, 4.305885125424644860264320635, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi / 3.0, 0.50, mode, 1.418324704155697579394036402, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi / 3.0, 0.01, mode, 1.263755353901126149206022061, TEST_TOL0, err.SUCCESS},

		{gsl.Pi, 0.99, mode, 4.750790152703577951330646444, TEST_TOL0, err.SUCCESS},
		{gsl.Pi, 0.50, mode, 1.746305163785351099291267125, TEST_TOL0, err.SUCCESS},
		{gsl.Pi, 0.01, mode, 1.570855235338973786559878750, TEST_TOL0, err.SUCCESS},

		{4 * gsl.Pi / 3, 0.99, mode, 5.195695179982511042396972113, TEST_TOL0, err.SUCCESS},
		{4 * gsl.Pi / 3, 0.50, mode, 2.074285623415004619188497818, TEST_TOL0, err.SUCCESS},
		{4 * gsl.Pi / 3, 0.01, mode, 1.877955116776821423913735408, TEST_TOL0, err.SUCCESS},

		{3 * gsl.Pi / 2.0, 0.99, mode, 7.126185229055366926995969476, TEST_TOL0, err.SUCCESS},
		{3 * gsl.Pi / 2.0, 0.50, mode, 2.619457745678026648936900687, TEST_TOL0, err.SUCCESS},
		{3 * gsl.Pi / 2.0, 0.01, mode, 2.356282853008460679839818125, TEST_TOL0, err.SUCCESS},

		{5 * gsl.Pi / 3, 0.99, mode, 9.056675278128222811594967044, TEST_TOL0, err.SUCCESS},
		{5 * gsl.Pi / 3, 0.50, mode, 3.164629867941048678685303509, TEST_TOL0, err.SUCCESS},
		{5 * gsl.Pi / 3, 0.01, mode, 2.834610589240099935765900794, TEST_TOL0, err.SUCCESS},

		{2 * gsl.Pi, 0.99, mode, 9.501580305407155902661292832, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi, 0.50, mode, 3.492610327570702198582534249, TEST_TOL0, err.SUCCESS},
		{2 * gsl.Pi, 0.01, mode, 3.141710470677947573119757500, TEST_TOL0, err.SUCCESS},

		{7 * gsl.Pi / 3.0, 0.99, mode, 9.946485332686088993727618315, TEST_TOL0, err.SUCCESS},
		{7 * gsl.Pi / 3.0, 0.50, mode, 3.820590787200355718479764901, TEST_TOL0, err.SUCCESS},
		{7 * gsl.Pi / 3.0, 0.01, mode, 3.448810352115795210473614120, TEST_TOL0, err.SUCCESS},

		/* P, negative argument phi < 0  */

		{-gsl.Pi / 2.0, 0.99, mode, -2.375395076351788975665323192, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi / 2.0, 0.50, mode, -0.8731525818926755496456335628, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi / 2.0, 0.01, mode, -0.7854276176694868932799393751, TEST_TOL0, err.SUCCESS},

		{-2 * gsl.Pi / 3.0, 0.99, mode, -4.305885125424644860264320635, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi / 3.0, 0.50, mode, -1.418324704155697579394036402, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi / 3.0, 0.01, mode, -1.263755353901126149206022061, TEST_TOL0, err.SUCCESS},

		{-gsl.Pi, 0.99, mode, -4.750790152703577951330646444, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi, 0.50, mode, -1.746305163785351099291267125, TEST_TOL0, err.SUCCESS},
		{-gsl.Pi, 0.01, mode, -1.570855235338973786559878750, TEST_TOL0, err.SUCCESS},

		{-4 * gsl.Pi / 3, 0.99, mode, -5.195695179982511042396972113, TEST_TOL0, err.SUCCESS},
		{-4 * gsl.Pi / 3, 0.50, mode, -2.074285623415004619188497818, TEST_TOL0, err.SUCCESS},
		{-4 * gsl.Pi / 3, 0.01, mode, -1.877955116776821423913735408, TEST_TOL0, err.SUCCESS},

		{-3 * gsl.Pi / 2.0, 0.99, mode, -7.126185229055366926995969476, TEST_TOL0, err.SUCCESS},
		{-3 * gsl.Pi / 2.0, 0.50, mode, -2.619457745678026648936900687, TEST_TOL0, err.SUCCESS},
		{-3 * gsl.Pi / 2.0, 0.01, mode, -2.356282853008460679839818125, TEST_TOL0, err.SUCCESS},

		{-5 * gsl.Pi / 3, 0.99, mode, -9.056675278128222811594967044, TEST_TOL0, err.SUCCESS},
		{-5 * gsl.Pi / 3, 0.50, mode, -3.164629867941048678685303509, TEST_TOL0, err.SUCCESS},
		{-5 * gsl.Pi / 3, 0.01, mode, -2.834610589240099935765900794, TEST_TOL0, err.SUCCESS},

		{-2 * gsl.Pi, 0.99, mode, -9.501580305407155902661292832, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi, 0.50, mode, -3.492610327570702198582534249, TEST_TOL0, err.SUCCESS},
		{-2 * gsl.Pi, 0.01, mode, -3.141710470677947573119757500, TEST_TOL0, err.SUCCESS},

		{-7 * gsl.Pi / 3.0, 0.99, mode, -9.946485332686088993727618315, TEST_TOL0, err.SUCCESS},
		{-7 * gsl.Pi / 3.0, 0.50, mode, -3.820590787200355718479764901, TEST_TOL0, err.SUCCESS},
		{-7 * gsl.Pi / 3.0, 0.01, mode, -3.448810352115795210473614120, TEST_TOL0, err.SUCCESS},
	}

	for i, c := range cases {
		t.Run(strconv.Itoa(i), func(t *testing.T) {
			stat := Ellint_D_e(c.phi, c.k, c.mode, r)
			run_test_sf(t, stat, r, c.expected, c.tol, c.status, fmt.Sprintf("Ellint_D_e(%v, %v, %v)", c.phi, c.k, c.mode))
		})
	}
}
