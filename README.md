# ggsl

## Golang port of GNU Scientific Library

Status: last sync with master - 12/19/2019

The license follows GPLv3 (As per the source). 

The method names have been stripped off of 'gsl_xxxxx_' (e.g. gsl_sf_poch_e() will become Poch_e(), gsl_integration_fixed() will become Fixed()) as it
should be known already (not only was it handful to type) that it is GSL.

Machine and Math constants have been renamed to match that of what Golang's naming convention in math package looks like.


## Usage

**ggsl.Function** is an interface that implements:

```golang
type Function interface {
	Evaluate(x float64) float64
}
```

This is the most basic type, as all packages would use this as much as possible, instead of func(float64) float64 (we'd love to have the freedom to use struct fields, embedding, etc.).

### integration package

#### [Quadrature](https://www.gnu.org/software/gsl/doc/html/integration.html#romberg-integration)

1. Qag(f ggsl.Function, a, b, epsabs, epsrel float64, limit int, workspace \*Workspace, result, abserr \*float64, q integration.QKFunction) err.GSLError
2. Qags(f ggsl.Function, a, b, epsabs, epsrel float64, limit int, workspace \*Workspace, result, abserr \*float64) err.GSLError
3. Qagi(f ggsl.Function, epsabs, epsrel float64, limit int, workspace \*Workspace, result, abserr \*float64) err.GSLError
4. Qagil(f ggsl.Function, b, epsabs, epsrel float64, limit int, workspace \*Workspace, result, abserr \*float64) err.GSLError
5. Qagiu(f ggsl.Function, a, epsabs, epsrel float64, limit int, workspace \*Workspace, result, abserr \*float64) err.GSLError
6. Qagp(f ggsl.Function, pts []float64, npts int, epsabs, epsrel float64, limit int, workspace \*Workspace, result, abserr \*float64, q integration.QKFunction) err.GSLError
7. Qawc(f ggsl.Function, a, b, c, epsabs, epsrel float64, limit int, workspace \*Workspace, result, abserr \*float64) err.GSLError
8. Qaws(f ggsl.Function, a, b float64, t \*QAWSTable, epsabs, epsrel float64, limit int, workspace \*Workspace, result, abserr \*float64) err.GSLError
9. Qawo(f ggsl.Function, a, epsabs, epsrel float64, limit int, workspace \*Workspace, wf \*QAWOTable, result, abserr \*float64) err.GSLError
10. Qawf(f ggsl.Function, a, epsabs float64, limit int, workspace, cycle_workspace \*Workspace, wf \*QAWOTable, result, abserr \*float64) err.GSLError


The low-level integration rules in ggsl are implementing:

```golang
	type QKFunction func(f gsl.Function, a, b float64, result, abserr, defabs, resabs *float64)
```

And has the following (abscissae of the N-point kronrod rule):

1. Qk15
2. Qk21
3. Qk31
4. Qk41
5. Qk51
6. Qk61


```golang
import (
	integ "github.com/lucky-se7en/ggsl/integration"
	"math"
)

type f1 struct {
	alpha float64
}

func (f f1) Evaluate(x float64) float64 {
	return math.Pow(x, f.alpha) * math.Log(1.0/x)
}

func main()  {
	alpha := 2.6
	f := &f1{alpha}

	workspace, _ := integ.NewWorkspace(1000)
	var result, abserr float64
	err := Qag(fc, 0.0, 1.0, 0.0, 1e-10, w.limit, w, &result, &abserr, integ.Qk15)
	fmt.Println(err)
}
```

#### [Romberg Integration](https://www.gnu.org/software/gsl/doc/html/integration.html#romberg-integration)
1. Romberg(f gsl.Function, a, b, epsabs, epsrel float64, result \*float64, neval \*int, w \*RombergWorkspace) err.GSLError


```golang
type f_sin struct {}

func (f f_sin) Evaluate(x float64) float64 {
	return math.Sin(x)
}

func main() {
	var result float64
	var neval int
	w, _ := integ.NewRombergWorkspace(20)
	f := f_sin{}
	err := integ.Romberg(&f, 0.0, gsl.PiOver2, 0.0, 1e-10, &result, &neval, w)
	fmt.Println(err)
}
```

#### [Gauss-Legendre Integration](https://www.gnu.org/software/gsl/doc/html/integration.html#gauss-legendre-integration)

1. GLFixed(f ggsl.Function, a, b float64, t \*GLFixedTable) float64
2. GLFixedPoint(a, b float64, i int, xi, wi \*float64, t \*GLFixedTable) err.GSLError

```golang
type f_sin struct {}

func (f f_sin) Evaluate(x float64) float64 {
	return math.Sin(x)
}

func main() {
    f := &f_sin{}
    a := 0.0
    b := gsl.Pi
    var result float64

	tbl, _ := integ.NewGLFixedTable(n)
	result = integ.GLFixed(f, a, b, tbl)
}
```

#### [Fixed point quadratures](https://www.gnu.org/software/gsl/doc/html/integration.html#fixed-point-quadratures)

1. Fixed(f gsl.Function, result \*float64, w \*FixedWorkspace) err.GSLError

The table below lists the weighting functions currently supported. They implement:
```golang
type FixedType interface {
	check(n int, params *fixedParams) err.GSLError
	init(n int, diag, subdiag []float64, params *fixedParams) err.GSLError
}
```
which is used internally by the integration method.
1. Legendre
2. ChebyshevType1
3. Gegenbauer
4. Jacobi
5. Laguerre
6. Hermite
7. Exponential
8. Rational
9. ChebyshevType2


```golang
type myfn1 struct {}

func (f myfn1) Evaluate(x float64) float64 {
	return math.Exp(-x - x*x)
}

func main() {
	var t integ.ChebyshevType2
	n := 200
	f := &myfn1{}
	w, _ := integ.NewFixedWorkspace(t, n, 1.2, 2.6, 0.0, 0.0)
	var result float64

	err := integ.Fixed(f2, &result, w)
	fmt.Println(err)
}
```

### specfunc package


All Special Functions are implemented except LegendreP (requires total reimplementation) and Mathieu functions (requires implementation of GSL's linear algebra).


### err package


Similarly to GSL, you can set your own error handler to prevent programs from terminating (see [err](https:github.com/lucky-se7en/ggsl/err) package):

```golang
func initTesting() {
	// If you'd like to know where it catches flu.
	err.SetErrorHandler(myErrorHandler)
	// or just disable it
	err.SetErrorHandlerOff()
	os.Setenv("GSL_TEST_VERBOSE", "1") // this is used for printing verbose test results
}

func myErrorHandler(reason, file string, line, gsl_errno int) {
	fmt.Printf("(error expected [%s:%d: %s (%d)])\n", file, line, reason, gsl_errno)
}
```

