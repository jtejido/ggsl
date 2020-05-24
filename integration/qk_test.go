package integration

import (
	"github.com/jtejido/ggsl/test"
	"testing"
)

func TestQk(t *testing.T) {
	initTesting()
	/* Test the basic Gauss-Kronrod rules with a smooth positive function. */
	{
		var result, abserr, resabs, resasc float64
		exp_result := 7.716049357767090777e-02
		exp_abserr := 2.990224871000550874e-06
		exp_resabs := 7.716049357767090777e-02
		exp_resasc := 4.434273814139995384e-02

		alpha := 2.6
		f := f1{alpha}

		Qk15(&f, 0.0, 1.0, &result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk15(f1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk15(f1) smooth abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk15(f1) smooth resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk15(f1) smooth resasc")

		Qk15(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)

		test.TestRel(t, result, -exp_result, 1e-15, "qk15(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk15(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk15(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk15(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 7.716049379303084599e-02
		exp_abserr := 9.424302194248481445e-08
		exp_resabs := 7.716049379303084599e-02
		exp_resasc := 4.434311425038358484e-02

		alpha := 2.6
		f := f1{alpha}

		Qk21(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk21(f1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk21(f1) smooth abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk21(f1) smooth resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk21(f1) smooth resasc")

		Qk21(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk21(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk21(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk21(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk21(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 7.716049382494900855e-02
		exp_abserr := 1.713503193600029893e-09
		exp_resabs := 7.716049382494900855e-02
		exp_resasc := 4.427995051868838933e-02

		alpha := 2.6
		f := f1{alpha}

		Qk31(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk31(f1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk31(f1) smooth abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk31(f1) smooth resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk31(f1) smooth resasc")

		Qk31(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk31(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk31(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk31(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk31(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 7.716049382681375302e-02
		exp_abserr := 9.576386660975511224e-11
		exp_resabs := 7.716049382681375302e-02
		exp_resasc := 4.421521169637691873e-02

		alpha := 2.6
		f := f1{alpha}

		Qk41(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk41(f1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk41(f1) smooth abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk41(f1) smooth resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk41(f1) smooth resasc")

		Qk41(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk41(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk41(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk41(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk41(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 7.716049382708510540e-02
		exp_abserr := 1.002079980317363772e-11
		exp_resabs := 7.716049382708510540e-02
		exp_resasc := 4.416474291216854892e-02

		alpha := 2.6
		f := f1{alpha}

		Qk51(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk51(f1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qk51(f1) smooth abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk51(f1) smooth resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk51(f1) smooth resasc")

		Qk51(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk51(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qk51(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk51(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk51(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 7.716049382713800753e-02
		exp_abserr := 1.566060362296155616e-12
		exp_resabs := 7.716049382713800753e-02
		exp_resasc := 4.419287685934316506e-02

		alpha := 2.6
		f := f1{alpha}

		Qk61(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk61(f1) smooth result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qk61(f1) smooth abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk61(f1) smooth resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk61(f1) smooth resasc")

		Qk61(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk61(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qk61(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk61(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk61(f1) reverse resasc")
	}

	/* Now test the basic rules with a positive function that has a
	   singularity. This should give large values of abserr which would
	   find discrepancies in the abserr calculation. */

	{
		var result, abserr, resabs, resasc float64
		exp_result := 1.555688196612745777e+01
		exp_abserr := 2.350164577239293706e+01
		exp_resabs := 1.555688196612745777e+01
		exp_resasc := 2.350164577239293706e+01

		alpha := -0.9
		f := f1{alpha}

		Qk15(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk15(f1) singular result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk15(f1) singular abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk15(f1) singular resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk15(f1) singular resasc")

		Qk15(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk15(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk15(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk15(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk15(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 1.799045317938126232e+01
		exp_abserr := 2.782360287710622515e+01
		exp_resabs := 1.799045317938126232e+01
		exp_resasc := 2.782360287710622515e+01

		alpha := -0.9
		f := f1{alpha}

		Qk21(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk21(f1) singular result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk21(f1) singular abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk21(f1) singular resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk21(f1) singular resasc")

		Qk21(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk21(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk21(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk21(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk21(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 2.081873305159121657e+01
		exp_abserr := 3.296500137482590276e+01
		exp_resabs := 2.081873305159121301e+01
		exp_resasc := 3.296500137482590276e+01

		alpha := -0.9
		f := f1{alpha}

		Qk31(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk31(f1) singular result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk31(f1) singular abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk31(f1) singular resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk31(f1) singular resasc")

		Qk31(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk31(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk31(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk31(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk31(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 2.288677623903126701e+01
		exp_abserr := 3.671538820274916048e+01
		exp_resabs := 2.288677623903126701e+01
		exp_resasc := 3.671538820274916048e+01

		alpha := -0.9
		f := f1{alpha}

		Qk41(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk41(f1) singular result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk41(f1) singular abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk41(f1) singular resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk41(f1) singular resasc")

		Qk41(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk41(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk41(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk41(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk41(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 2.449953612016972215e+01
		exp_abserr := 3.967771249391228849e+01
		exp_resabs := 2.449953612016972215e+01
		exp_resasc := 3.967771249391228849e+01

		alpha := -0.9
		f := f1{alpha}

		Qk51(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk51(f1) singular result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk51(f1) singular abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk51(f1) singular resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk51(f1) singular resasc")

		Qk51(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk51(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk51(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk51(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk51(f1) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := 2.583030240976628988e+01
		exp_abserr := 4.213750493076978643e+01
		exp_resabs := 2.583030240976628988e+01
		exp_resasc := 4.213750493076978643e+01

		alpha := -0.9
		f := f1{alpha}

		Qk61(&f, 0.0, 1.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk61(f1) singular result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk61(f1) singular abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk61(f1) singular resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk61(f1) singular resasc")

		Qk61(&f, 1.0, 0.0,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk61(f1) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk61(f1) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk61(f1) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk61(f1) reverse resasc")
	}

	/* Test the basic Gauss-Kronrod rules with a smooth oscillating
	   function, over an unsymmetric range. This should find any
	   discrepancies in the abscissae. */

	{
		var result, abserr, resabs, resasc float64
		exp_result := -7.238969575483799046e-01
		exp_abserr := 8.760080200939757174e-06
		exp_resabs := 1.165564172429140788e+00
		exp_resasc := 9.334560307787327371e-01

		alpha := 1.3
		f := f3{alpha}

		Qk15(&f, 0.3, 2.71,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk15(f3) oscill result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk15(f3) oscill abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk15(f3) oscill resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk15(f3) oscill resasc")

		Qk15(&f, 2.71, 0.3,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk15(f3) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk15(f3) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk15(f3) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk15(f3) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := -7.238969575482959717e-01
		exp_abserr := 7.999213141433641888e-11
		exp_resabs := 1.150829032708484023e+00
		exp_resasc := 9.297591249133687619e-01

		alpha := 1.3
		f := f3{alpha}

		Qk21(&f, 0.3, 2.71,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk21(f3) oscill result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qk21(f3) oscill abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk21(f3) oscill resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk21(f3) oscill resasc")

		Qk21(&f, 2.71, 0.3,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk21(f3) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-5, "qk21(f3) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk21(f3) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk21(f3) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := -7.238969575482959717e-01
		exp_abserr := 1.285805464427459261e-14
		exp_resabs := 1.158150602093290571e+00
		exp_resasc := 9.277828092501518853e-01

		alpha := 1.3
		f := f3{alpha}

		Qk31(&f, 0.3, 2.71,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk31(f3) oscill result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk31(f3) oscill abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk31(f3) oscill resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk31(f3) oscill resasc")

		Qk31(&f, 2.71, 0.3,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk31(f3) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk31(f3) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk31(f3) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk31(f3) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := -7.238969575482959717e-01
		exp_abserr := 1.286535726271015626e-14
		exp_resabs := 1.158808363486595328e+00
		exp_resasc := 9.264382258645686985e-01

		alpha := 1.3
		f := f3{alpha}

		Qk41(&f, 0.3, 2.71,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk41(f3) oscill result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk41(f3) oscill abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk41(f3) oscill resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk41(f3) oscill resasc")

		Qk41(&f, 2.71, 0.3,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk41(f3) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk41(f3) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk41(f3) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk41(f3) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := -7.238969575482961938e-01
		exp_abserr := 1.285290995039385778e-14
		exp_resabs := 1.157687209264406381e+00
		exp_resasc := 9.264666884071264263e-01

		alpha := 1.3
		f := f3{alpha}

		Qk51(&f, 0.3, 2.71,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk51(f3) oscill result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk51(f3) oscill abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk51(f3) oscill resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk51(f3) oscill resasc")

		Qk51(&f, 2.71, 0.3,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk51(f3) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk51(f3) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk51(f3) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk51(f3) reverse resasc")
	}

	{
		var result, abserr, resabs, resasc float64
		exp_result := -7.238969575482959717e-01
		exp_abserr := 1.286438572027470736e-14
		exp_resabs := 1.158720854723590099e+00
		exp_resasc := 9.270469641771273972e-01

		alpha := 1.3
		f := f3{alpha}

		Qk61(&f, 0.3, 2.71,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, exp_result, 1e-15, "qk61(f3) oscill result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk61(f3) oscill abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk61(f3) oscill resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk61(f3) oscill resasc")

		Qk61(&f, 2.71, 0.3,
			&result, &abserr, &resabs, &resasc)
		test.TestRel(t, result, -exp_result, 1e-15, "qk61(f3) reverse result")
		test.TestRel(t, abserr, exp_abserr, 1e-7, "qk61(f3) reverse abserr")
		test.TestRel(t, resabs, exp_resabs, 1e-15, "qk61(f3) reverse resabs")
		test.TestRel(t, resasc, exp_resasc, 1e-15, "qk61(f3) reverse resasc")
	}
}
