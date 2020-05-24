package integration

import (
	"fmt"
	gsl "github.com/lucky-se7en/ggsl"
	"github.com/lucky-se7en/ggsl/specfunc"
	"github.com/lucky-se7en/ggsl/test"
	"math"
	"testing"
)

func test_fixed_quadrature(t *testing.T, T FixedType, n int, a, b, alpha, beta, tol, exact float64, f gsl.Function, desc string) int {
	var status int
	w, _ := NewFixedWorkspace(T, n, a, b, alpha, beta)
	var result float64

	str := fmt.Sprintf("%s a=%g b=%g alpha=%g beta=%g", desc, a, b, alpha, beta)

	Fixed(f, &result, w)
	test.TestRel(t, result, exact, tol, str)
	return status
}

func TestFixed(t *testing.T) {
	initTesting()
	/* test fixed quadrature */
	{

		var exact, a, b float64
		deg := 5 /* monomial degree */
		var dterm float64
		if (deg % 2) == 0 {
			dterm = 1.
		} else {
			dterm = -1.
		}

		params := &monomial_params{deg, 1.0}
		f := f_monomial{params}
		n := 15
		for b = 1.1; b <= 4.0; b += 0.1 {
			/* test with a < b */
			a = b - 1.0
			/* Legendre quadrature */
			exact = (math.Pow(b, float64(params.degree)+1.0) - math.Pow(a, float64(params.degree)+1.0)) / (float64(params.degree) + 1.0)
			test_fixed_quadrature(t, Legendre{}, n, a, b, 0.0, 0.0, 1.0e-12, exact, &f, "legendre monomial")

			/* Chebyshev type 1 quadrature */
			exact = gsl.Sign(b-a) * gsl.Pi * math.Pow(0.5*(a+b), float64(params.degree)) * specfunc.Hyperg_2F1(0.5*(1-float64(params.degree)), -0.5*float64(params.degree), 1.0, (b-a)*(b-a)/((b+a)*(b+a)))
			test_fixed_quadrature(t, ChebyshevType1{}, n, a, b, 0.0, 0.0, 1.0e-12, exact, &f, "chebyshev monomial")

			/* Laguerre quadrature */
			exact = math.Pow(b, -1.0-float64(deg)) * math.Exp(a*b) * specfunc.Gamma_inc(1.0+float64(deg), a*b)
			test_fixed_quadrature(t, Laguerre{}, n, a, b, 0.0, 0.0, 1.0e-12, exact, &f, "laguerre monomial")

			/* Hermite quadrature */
			exact = 0.5 * math.Pow(b, -0.5*float64(deg)) * (-(-1.0+dterm)*a*float64(deg)*specfunc.Gamma(0.5*float64(deg))*specfunc.Hyperg_1F1(0.5-0.5*float64(deg), 1.5, -a*a*b) +
				(1.0+dterm)*specfunc.Gamma(0.5*(1.0+float64(deg)))*specfunc.Hyperg_1F1(-0.5*float64(deg), 0.5, -a*a*b)/math.Sqrt(b))
			test_fixed_quadrature(t, Hermite{}, n, a, b, 0.0, 0.0, 1.0e-12, exact, &f, "hermite monomial")

			/* Chebyshev type 2 quadrature */
			exact = gsl.Sign(b-a) * gsl.PiOver2 * math.Pow(0.5*(a+b), float64(params.degree)) * specfunc.Hyperg_2F1(0.5*(1-float64(params.degree)), -0.5*float64(params.degree), 2.0, (b-a)*(b-a)/((b+a)*(b+a))) * 0.25 * (b - a) * (b - a)
			test_fixed_quadrature(t, ChebyshevType2{}, n, a, b, 0.0, 0.0, 1.0e-12, exact, &f, "chebyshev2 monomial")

			/* now test with a > b */
			a = b + 1.0

			/* Legendre quadrature */
			exact = (math.Pow(b, float64(params.degree)+1.0) - math.Pow(a, float64(params.degree)+1.0)) / (float64(params.degree) + 1.0)
			test_fixed_quadrature(t, Legendre{}, n, a, b, 0.0, 0.0, 1.0e-12, exact, &f, "legendre monomial")

			/* Laguerre quadrature */
			exact = math.Pow(b, -1.0-float64(deg)) * math.Exp(a*b) * specfunc.Gamma_inc(1.0+float64(deg), a*b)
			test_fixed_quadrature(t, Laguerre{}, n, a, b, 0.0, 0.0, 1.0e-12, exact, &f, "laguerre monomial")

			/* Hermite quadrature */
			exact = 0.5 * math.Pow(b, -0.5*float64(deg)) * (-(-1.0+dterm)*a*float64(deg)*specfunc.Gamma(0.5*float64(deg))*specfunc.Hyperg_1F1(0.5-0.5*float64(deg), 1.5, -a*a*b) +
				(1.0+dterm)*specfunc.Gamma(0.5*(1.0+float64(deg)))*specfunc.Hyperg_1F1(-0.5*float64(deg), 0.5, -a*a*b)/math.Sqrt(b))
			test_fixed_quadrature(t, Hermite{}, n, a, b, 0.0, 0.0, 1.0e-12, exact, &f, "hermite monomial")

		}

		/* now test on myfn1 */
		f2 := myfn1{}
		n = 200

		test_fixed_quadrature(t, Legendre{}, n, 1.2, 1.6, 0.0, 0.0, 1.0e-12, 0.01505500344456001, &f2, "legendre myfn1")
		test_fixed_quadrature(t, ChebyshevType1{}, n, 1.2, 2.6, 0.0, 0.0, 1.0e-12, 0.0582346516219999, &f2, "chebyshev myfn1")
		test_fixed_quadrature(t, Gegenbauer{}, n, 1.2, 1.6, 2.0, 0.0, 1.0e-12, 1.2279468957162412661311711271e-5, &f2, "gegenbauer myfn1")
		test_fixed_quadrature(t, Gegenbauer{}, n, 1.2, 1.6, -0.5, 0.0, 1.0e-12, 1.228256086101808986e-1, &f2, "gegenbauer myfn1")
		test_fixed_quadrature(t, Jacobi{}, n, 1.2, 1.6, 2.0, 1.5, 1.0e-12, 3.173064776410033e-5, &f2, "jacobi myfn1")
		test_fixed_quadrature(t, Jacobi{}, n, 1.2, 1.6, -0.5, -0.5, 1.0e-12, 1.228256086101808986e-1, &f2, "jacobi myfn1")
		test_fixed_quadrature(t, Laguerre{}, n, 1.2, 0.6, 0.5, 0.0, 1.0e-12, 0.006604180366378123, &f2, "laguerre myfn1")
		test_fixed_quadrature(t, Hermite{}, n, 1.2, 0.6, 1.0, 0.0, 1.0e-12, 0.6542819629825344, &f2, "hermite myfn1")
		test_fixed_quadrature(t, Exponential{}, n, 1.2, 1.6, 2.0, 0.0, 1.0e-12, 2.1315535492168832898083633e-4, &f2, "exponential myfn1")
		test_fixed_quadrature(t, Rational{}, 15, 1.2, 1.6, 2.0, -33.4, 1.0e-9, 4.8457468060064844e-20, &f2, "rational myfn1")
		test_fixed_quadrature(t, ChebyshevType2{}, n, 1.2, 2.6, 0.0, 0.0, 1.0e-12, 0.0081704088896491, &f2, "chebyshev2 myfn1")
	}

	/* test Gegenbauer quadrature */
	{
		exactarray := []float64{4.15933612154155020161400717857e-7, 744697.808572324010134504819452, 55.2024994284578980512106835228, 7.95574829722734114107142857143, 0.00179653588816666666666666666667}
		aarray := []float64{0.123, 7.747, 1.47, -1.47, 0.0}
		barray := []float64{0.456, 12.0, 2.0, 2.0, 0.47}
		alphaarray := []float64{2.0, 0.5, -0.5, 1.0, 0.0}

		params := &monomial_params{5, 1.0}
		f := f_monomial{params}

		n := 50
		for k := 0; k < 5; k++ {
			test_fixed_quadrature(t, Gegenbauer{}, n, aarray[k], barray[k], alphaarray[k], 0.0, 1.0e-12, exactarray[k], &f, "gegenbauer monomial")
		}
	}

	/* test Jacobi quadrature */
	{

		exactarray := []float64{9.052430592016123480501898e-7, 3.131716150347619771233591755e6, 0.04435866422797298224404592896, 5.287059602300844442782407, 2.5337038518475893688512749675e-6}
		aarray := []float64{0.123, 7.747, 1.47, -1.47, 0.0}
		barray := []float64{0.456, 12.0, 2.0, 2.0, 0.47}
		var alpha, beta float64

		params := &monomial_params{5, 1.0}
		f := f_monomial{params}
		alpha = 2.0
		beta = 1.5
		n := 50
		for k := 0; k < 5; k++ {
			test_fixed_quadrature(t, Jacobi{}, n, aarray[k], barray[k], alpha, beta, 1.0e-12, exactarray[k], &f, "jacobi monomial")
		}
	}

	/* test Exponential quadrature */
	{
		exactarray := []float64{1.598864206823942764921875e-4, 624615.81848571833291063083819, 0.222578063871903188095238095238, 28.8968950008739567709168294271, 4.62725113500425479890950520833e-7}
		aarray := []float64{0.123, 7.747, 1.47, -1.47, 0.0}
		barray := []float64{0.456, 12.0, 2.0, 2.0, 0.47}
		alphaarray := []float64{1.0, 1.5, 2.0, 3.0, 5.0}

		params := &monomial_params{5, 1.0}
		f := f_monomial{params}

		n := 50
		for k := 0; k < 5; k++ {
			test_fixed_quadrature(t, Exponential{}, n, aarray[k], barray[k], alphaarray[k], 0.0, 1.0e-12, exactarray[k], &f, "exponential monomial")
		}
	}

	/* test Rational quadrature */
	{
		exactarray := []float64{1.312245361412108703130374957e-10, 0.0170362044485924082779613124672, 8.93065131938394658578136414201e-11, 7.17990217357447544326794457270e-13, -11.0760676986664098133970869634, 0.00290392485414197833688178206557}
		aarray := []float64{0.0, 0.123, 7.747, 1.47, -1.47, 0.0}
		barray := []float64{2.0, 0.456, 12.0, 2.0, 2.0, 0.47}
		alphaarray := []float64{0.0, 1.0, 1.5, 2.0, 3.0, 5.0}
		betaarray := []float64{-21.0, -12.0, -13.0, -22.0, -21.0, -16.0}

		params := &monomial_params{5, 1.0}
		f := f_monomial{params}

		n := 5
		for k := 0; k < 5; k++ {
			test_fixed_quadrature(t, Rational{}, n, aarray[k], barray[k], alphaarray[k], betaarray[k], 1.0e-12, exactarray[k], &f, "rational monomial")
		}
	}
}
