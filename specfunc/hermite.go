/* specfunc/hermite.c
 *
 * Copyright (C) 2011, 2012, 2013, 2014, 2019 Konrad Griessinger (konradg(at)gmx.net)
 * Copyright (C) 2019 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
package specfunc

import (
	gsl "github.com/jtejido/ggsl"
	"github.com/jtejido/ggsl/err"
	"github.com/jtejido/ggsl/sys"
	"math"
)

/*----------------------------------------------------------------------*
 * "The purpose of computing is insight, not numbers." - R.W. Hamming   *
 * Hermite polynomials, Hermite functions                               *
 * and their respective arbitrary derivatives                           *
 *----------------------------------------------------------------------*/

/* TODO:
 * - array functions for derivatives of Hermite functions
 * - asymptotic approximation for derivatives of Hermite functions
 * - refine existing asymptotic approximations, especially around x=math.Sqrt(2*n+1) or x=math.Sqrt(2*n+1)*math.Sqrt(2), respectively
 */

func pow2(n int) float64 { return Pow_int(2, n) }
func rnd(x float64) float64 {
	if x >= 0 {
		return float64(int(x + 0.5))
	}

	return float64(int(x - 0.5))
}

/* evaluates the probabilists' Hermite polynomial of order n at position x */
func Hermite_prob_e(n int, x float64, result *Result) err.GSLError {
	if n < 0 {
		return DomainError(result)
	} else if n == 0 {
		result.val = 1.
		result.err = 0.
		return nil
	} else if n == 1 {
		result.val = x
		result.err = 0.
		return nil
	} else if x == 0.0 {
		if gsl.IsOdd(n) {
			result.val = 0.0
			result.err = 0.0
			return nil
		} else {
			/* for n even, He_n(0) = (-1)^{n/2} (n - 1)!! */
			var status err.GSLError
			if n-1 > DOUBLEFACT_NMAX { /* test if (n - 1)!! will overflow */
				status = err.Overflow()
				if gsl.IsOdd(n / 2) {
					result.val = math.Inf(-1)
				} else {
					result.val = math.Inf(1)
				}
				result.err = math.Inf(1)
			} else {
				Doublefact_e(uint(n-1), result)
				if gsl.IsOdd(n / 2) {
					result.val = -result.val
				}
			}

			return status
		}
	} else {
		/* upward recurrence: He_{n+1} = x He_n - n He_{n-1} */

		var status err.GSLError
		abs_x := math.Abs(x)
		thresh1 := gsl.MaxFloat64
		if abs_x > 1.0 {
			thresh1 = 0.9 * gsl.MaxFloat64 / abs_x
		}
		thresh2 := 0.9 * gsl.MaxFloat64

		p_n0 := 1.0 /* He_0(x) */
		p_n1 := x   /* He_1(x) */
		p_n := p_n1

		e_n0 := gsl.Float64Eps
		e_n1 := math.Abs(x) * gsl.Float64Eps
		e_n := e_n1

		var j int

		for j = 1; j < n; j++ {
			fj := float64(j)
			if math.Abs(p_n1) > thresh1 || math.Abs(p_n0) > thresh2/fj {
				status = err.Overflow()
				break
			}

			p_n = x*p_n1 - fj*p_n0
			p_n0 = p_n1
			p_n1 = p_n

			e_n = math.Abs(x)*e_n1 + fj*e_n0
			e_n0 = e_n1
			e_n1 = e_n
		}

		result.val = p_n
		result.err = e_n + math.Abs(result.val)*gsl.Float64Eps

		return status
	}
}

func Hermite_prob(n int, x float64) float64 {
	result := new(Result)
	status := Hermite_prob_e(n, x, result)
	return EvalResult(result, status)
}

/* Evaluates the m-th derivative of the probabilists' Hermite polynomial of order n at position x.
 * The direct formula He^{(m)}_n = n!/(n-m)!*He_{n-m}(x) (where He_j(x) is the j-th probabilists' Hermite polynomial and He^{(m)}_j(x) its m-th derivative) is employed. */
func Hermite_prob_deriv_e(m, n int, x float64, result *Result) err.GSLError {
	if n < 0 || m < 0 {
		return DomainError(result)
	} else if n < m {
		result.val = 0.
		result.err = 0.
		return nil
	} else {
		var status err.GSLError
		f := Choose(uint(n), uint(m)) * Fact(uint(m))
		var He Result

		status = Hermite_prob_e(n-m, x, &He)
		if status == nil {
			result.val = He.val * f
			result.err = He.err*f + gsl.Float64Eps*math.Abs(result.val)
		} else {
			result.val = He.val
			result.err = math.Inf(1)
		}

		return status
	}
}

func Hermite_prob_deriv(m, n int, x float64) float64 {
	result := new(Result)
	status := Hermite_prob_deriv_e(m, n, x, result)
	return EvalResult(result, status)
}

/* evaluates the physicists' Hermite polynomial of order n at position x */
func Hermite_e(n int, x float64, result *Result) err.GSLError {
	if n < 0 {
		return DomainError(result)
	} else if n == 0 {
		result.val = 1.
		result.err = 0.
		return nil
	} else if n == 1 {
		result.val = 2.0 * x
		result.err = 0.
		return nil
	} else if x == 0.0 {
		if gsl.IsOdd(n) {
			result.val = 0.0
			result.err = 0.0
			return nil
		} else {
			/* for n even, H_n(0) = (-2)^{n/2} (n - 1)!! */

			var status err.GSLError
			m := n >> 1

			if n-1 > DOUBLEFACT_NMAX { /* test if (n - 1)!! will overflow */
				status = err.Overflow()
				result.val = math.Inf(1)
				if gsl.IsOdd(m) {
					result.val = math.Inf(-1)
				}
				result.err = math.Inf(1)
			} else {
				f := Pow_int(2.0, m)

				Doublefact_e(uint(n-1), result)

				if result.val > 0.9*gsl.MaxFloat64/f { /* test if 2^{n/2} * (n-1)!! will overflow */
					status = err.Overflow()
					result.val = math.Inf(1)
					if gsl.IsOdd(m) {
						result.val = math.Inf(-1)
					}
					result.err = math.Inf(1)
				} else {
					result.val *= f
					result.err *= f
					if gsl.IsOdd(m) {
						result.val = -result.val
					}
				}
			}

			return status
		}
	} else {
		/* upward recurrence: H_{n+1} = 2x H_n - 2n H_{n-1} */

		var status err.GSLError
		two_x := 2.0 * x
		abs_two_x := math.Abs(two_x)
		thresh1 := gsl.MaxFloat64
		if abs_two_x > 1.0 {
			thresh1 = 0.9 * gsl.MaxFloat64 / abs_two_x
		}
		thresh2 := 0.9 * gsl.MaxFloat64 / 2.0

		p_n0 := 1.0   /* H_0(x) */
		p_n1 := two_x /* H_1(x) */
		p_n := p_n1

		e_n0 := gsl.Float64Eps
		e_n1 := 2. * math.Abs(x) * gsl.Float64Eps
		e_n := e_n1

		var j int

		for j = 1; j <= n-1; j++ {
			fj := float64(j)
			if math.Abs(p_n1) > thresh1 || math.Abs(p_n0) > thresh2/fj {
				status = err.Overflow()
				break
			}

			p_n = two_x*p_n1 - 2.0*fj*p_n0
			p_n0 = p_n1
			p_n1 = p_n

			e_n = 2. * (math.Abs(x)*e_n1 + fj*e_n0)
			e_n0 = e_n1
			e_n1 = e_n
		}

		result.val = p_n
		result.err = e_n + math.Abs(result.val)*gsl.Float64Eps

		return status
	}
}

func Hermite(n int, x float64) float64 {
	result := new(Result)
	status := Hermite_e(n, x, result)
	return EvalResult(result, status)
}

/* Evaluates the m-th derivative of the physicists' Hermite polynomial of order n at position x.
 * The direct formula H^{(m)}_n = 2**m*n!/(n-m)!*H_{n-m}(x) (where H_j(x) is the j-th physicists' Hermite polynomial and H^{(m)}_j(x) its m-th derivative) is employed. */
func Hermite_deriv_e(m, n int, x float64, result *Result) err.GSLError {
	if n < 0 || m < 0 {
		return DomainError(result)
	} else if n < m {
		result.val = 0.
		result.err = 0.
		return nil
	} else {
		var status err.GSLError
		f := Choose(uint(n), uint(m)) * Fact(uint(m)) * pow2(m)
		var H Result

		status = Hermite_e(n-m, x, &H)
		if status == nil {
			result.val = H.val * f
			result.err = H.err*f + gsl.Float64Eps*math.Abs(result.val)
		} else {
			result.val = H.val
			result.err = math.Inf(1)
		}

		return status
	}
}

func Hermite_deriv(m, n int, x float64) float64 {
	result := new(Result)
	status := Hermite_deriv_e(m, n, x, result)
	return EvalResult(result, status)
}

/* evaluates the Hermite function of order n at position x */
func Hermite_func_e(n int, x float64, result *Result) err.GSLError {
	if n < 0 {
		return DomainError(result)
	} else if x == 0.0 {
		if gsl.IsOdd(n) {
			result.val = 0.
			result.err = 0.
			return nil
		} else {
			f := 1.0
			if gsl.IsOdd(n / 2) {
				f = -1.0
			}
			var j int

			for j = 1; j < n; j += 2 {
				fj := float64(j)
				f *= math.Sqrt(fj / (fj + 1.0))
			}

			result.val = f / math.Sqrt(math.SqrtPi)
			result.err = gsl.Float64Eps * math.Abs(result.val)
			return nil
		}
	} else if n == 0 {
		result.val = math.Exp(-0.5*x*x) / math.Sqrt(math.SqrtPi)
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if n == 1 {
		result.val = math.Sqrt2 * x * math.Exp(-0.5*x*x) / math.Sqrt(math.SqrtPi)
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else {
		/*
		 * This algorithm is based on the modified recurrence algorithm
		 * found in the appendix of:
		 *
		 * B. Bunck, BIT Numerical Mathematics, 49, 281 (2009)
		 *
		 * Numerical tests showed that this algorithm is more stable for
		 * large x (x > 40) than the standard recurrence relation.
		 * Accuracy is comparable to the recurrence relation method
		 * for small x and all n.
		 *
		 * See:
		 *
		 * https://scicomp.stackexchange.com/questions/30896/generate-high-n-quantum-harmonic-oscillator-states-numerically
		 *
		 * for further discussion.
		 */

		hi2 := 1.0 / math.Sqrt(math.SqrtPi) /* \hat{h}_0 */
		hi1 := math.Sqrt2 * x * hi2         /* \hat{h}_1 */
		hi := 0.0
		sum_log_scale := 0.0
		var abshi float64
		var i int

		for i = 2; i <= n; i++ {
			fi := float64(i)
			hi = math.Sqrt(2.0/fi)*x*hi1 - math.Sqrt((fi-1.0)/fi)*hi2
			hi2 = hi1
			hi1 = hi

			abshi = math.Abs(hi)
			if abshi > 1.0 {
				log_scale := rnd(math.Log(abshi))
				scale := math.Exp(-log_scale)

				hi *= scale
				hi1 *= scale
				hi2 *= scale
				sum_log_scale += log_scale
			}
		}

		result.val = hi * math.Exp(-0.5*x*x+sum_log_scale)
		result.err = float64(n) * gsl.Float64Eps * math.Abs(result.val)

		return nil
	}
}

func Hermite_func(n int, x float64) float64 {
	result := new(Result)
	status := Hermite_func_e(n, x, result)
	return EvalResult(result, status)
}

/*
 * This algorithm is based on the contour integral algorithm of:
 *
 * B. Bunck, BIT Numerical Mathematics, 49, 281 (2009)
 *
 * It has O(math.Sqrt(n)) complexity
 */
func Hermite_func_fast_e(n int, x float64, result *Result) err.GSLError {
	if n < 1000 || x == 0.0 {
		/* for small n, the recurrence method is faster and more accurate */
		return Hermite_func_e(n, x, result)
	} else {
		var j int
		fn := float64(n)
		k := math.Sqrt(0.5 * fn)
		steps := int(math.Ceil(6.211 * math.Sqrt(fn)))
		dt := gsl.Pi / float64(steps)
		invn2 := 1.0 / (fn * fn)
		var ex, ex_e, cs, cs_e, sn, sn2, t float64
		var lngamma Result

		if n < 36 {
			Lnfact_e(uint(n), &lngamma)
			lngamma.val *= 0.5
			lngamma.err *= 0.5
			t = 0.5*fn*math.Log(fn) + 0.25*gsl.LnPi
			cs = 0.5 * fn
			lngamma.val += cs - t
			lngamma.err += (cs + t) * gsl.Float64Eps
		} else {
			/* approximate ln(gamma_{n,k}) using Stirling's formula */
			lngamma.val = 0.25 * math.Log(2*fn)
			lngamma.err = (lngamma.val + ((((invn2/3360+1.0/2520)*invn2+1.0/720)*invn2)+1.0/24)/fn) * gsl.Float64Eps
			lngamma.val -= ((((invn2/3360-1.0/2520)*invn2 + 1.0/720) * invn2) - 1.0/24) / fn
		}

		ex = math.Exp(lngamma.val - fn - 0.5*x*x - 2*x*k)
		cs = 1
		if gsl.IsOdd(n) {
			cs = -1
		}
		result.val = 0.5 * ex * cs
		result.err = 0.5 * ex * (lngamma.err + (fn+0.5*x*x+math.Abs(2*x*k)+1)*gsl.Float64Eps)
		ex = math.Exp(lngamma.val - fn - 0.5*x*x + 2*x*k)
		result.val += 0.5 * ex
		result.err += 0.5 * ex * (lngamma.err + (fn+0.5*x*x+math.Abs(2*x*k)+1)*gsl.Float64Eps)

		for j = 1; j < steps; j++ {
			fk := float64(k)
			t = float64(j) * dt
			cs = math.Cos(t)
			ex = math.Exp(lngamma.val - 0.5*x*x + (2*x*fk-fn*cs)*cs)
			ex_e = ex * (lngamma.err + gsl.Float64Eps*(1+0.5*x*x+(math.Abs(2*x*fk)+math.Abs(fn*cs))*math.Abs(cs)))
			sn = math.Sin(t)
			sn2 = math.Sin(2 * t)
			cs = math.Cos(2*x*fk*sn - 0.5*fn*sn2 - fn*t)
			cs_e = math.Min(1.0+math.Abs(cs), gsl.Float64Eps*(math.Abs(cs)+(math.Abs(2*x*fk*sn)+math.Abs(0.5*fn*sn2)+fn*t)*math.Abs(math.Sin(2*x*fk*sn-0.5*fn*sn2-fn*t))))
			result.val += ex * cs
			result.err += ex*cs_e + ex_e*math.Abs(cs) + gsl.Float64Eps*math.Abs(ex*cs)
		}

		result.val *= gsl.OneOverPi * dt
		result.err = gsl.OneOverPi*dt*result.err + gsl.Float64Eps*math.Abs(result.val)

		return nil
	}
}

func Hermite_func_fast(n int, x float64) float64 {
	result := new(Result)
	status := Hermite_func_fast_e(n, x, result)
	return EvalResult(result, status)
}

/* Evaluates all probabilists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
func Hermite_prob_array(nmax int, x float64, result_array []float64) err.GSLError {
	if nmax < 0 {
		return err.ERROR("domain error", err.EDOM)
	} else if nmax == 0 {
		result_array[0] = 1.0
		return nil
	} else if nmax == 1 {
		result_array[0] = 1.0
		result_array[1] = x
		return nil
	} else {
		/* upward recurrence: He_{n+1} = x He_n - n He_{n-1} */
		var status err.GSLError
		abs_x := math.Abs(x)
		thresh1 := gsl.MaxFloat64
		if abs_x > 1.0 {
			thresh1 = 0.9 * gsl.MaxFloat64 / abs_x
		}
		thresh2 := 0.9 * gsl.MaxFloat64
		p_n0 := 1.0 /* He_0(x) */
		p_n1 := x   /* He_1(x) */
		p_n := p_n1

		var j int

		result_array[0] = 1.0
		result_array[1] = x

		for j = 1; j < nmax; j++ {
			fj := float64(j)
			if math.Abs(p_n1) > thresh1 || math.Abs(p_n0) > thresh2/fj {
				status = err.Overflow()
				break
			}

			p_n = x*p_n1 - fj*p_n0
			p_n0 = p_n1
			p_n1 = p_n

			result_array[j+1] = p_n
		}

		return status
	}
}

/* Evaluates the m-th derivative of all probabilists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
func Hermite_prob_array_deriv(m, nmax int, x float64, result_array []float64) err.GSLError {
	if nmax < 0 || m < 0 {
		return err.ERROR("domain error", err.EDOM)
	} else if m == 0 {
		Hermite_prob_array(nmax, x, result_array)
		return nil
	} else if nmax < m {
		var j int

		for j = 0; j <= nmax; j++ {
			result_array[j] = 0.0
		}

		return nil
	} else if nmax == m {
		var j int

		for j = 0; j < m; j++ {
			result_array[j] = 0.0
		}

		result_array[nmax] = Fact(uint(m))
		return nil
	} else if nmax == m+1 {
		var j int

		for j = 0; j < m; j++ {
			result_array[j] = 0.0
		}

		result_array[nmax-1] = Fact(uint(m))
		result_array[nmax] = result_array[nmax-1] * float64(m+1) * x
		return nil
	} else {
		/* upward recurrence: He^{(m)}_{n+1} = (n+1)/(n-m+1)*(x He^{(m)}_n - n He^{(m)}_{n-1}) */

		p_n0 := Fact(uint(m))           /* He^{(m)}_{m}(x) */
		p_n1 := p_n0 * float64(m+1) * x /* He^{(m)}_{m+1}(x) */
		p_n := p_n1
		var j int

		for j = 0; j < m; j++ {
			result_array[j] = 0.0
		}

		result_array[m] = p_n0
		result_array[m+1] = p_n1

		for j = m + 1; j <= nmax-1; j++ {
			fj := float64(j)
			fm := float64(m)
			p_n = (x*p_n1 - fj*p_n0) * (fj + 1.0) / (fj - fm + 1)
			p_n0 = p_n1
			p_n1 = p_n
			result_array[j+1] = p_n
		}

		return nil
	}
}

/* Evaluates all derivatives (starting from 0) up to the mmax-th derivative of the probabilists' Hermite polynomial of order n at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
func Hermite_prob_deriv_array(mmax, n int, x float64, result_array []float64) err.GSLError {
	if n < 0 || mmax < 0 {
		return err.ERROR("domain error", err.EDOM)
	} else if n == 0 {
		var j int

		result_array[0] = 1.0
		for j = 1; j <= mmax; j++ {
			result_array[j] = 0.0
		}

		return nil
	} else if n == 1 && mmax > 0 {
		var j int

		result_array[0] = x
		result_array[1] = 1.0

		for j = 2; j <= mmax; j++ {
			result_array[j] = 0.0
		}

		return nil
	} else if mmax == 0 {
		result_array[0] = Hermite_prob(n, x)
		return nil
	} else if mmax == 1 {
		result_array[0] = Hermite_prob(n, x)
		result_array[1] = float64(n) * Hermite_prob(n-1, x)
		return nil
	} else {
		/* upward recurrence */
		k := int(gsl.Max(0, n-mmax))
		p_n0 := Hermite_prob(k, x)   /* He_k(x) */
		p_n1 := Hermite_prob(k+1, x) /* He_{k+1}(x) */
		p_n := p_n1
		var j int

		for j = n + 1; j <= mmax; j++ {
			result_array[j] = 0.0
		}

		result_array[int(gsl.Min(n, mmax))] = p_n0
		result_array[int(gsl.Min(n, mmax))-1] = p_n1

		for j = int(gsl.Min(mmax, n)) - 1; j > 0; j-- {
			k++
			p_n = x*p_n1 - float64(k)*p_n0
			p_n0 = p_n1
			p_n1 = p_n
			result_array[j-1] = p_n
		}

		p_n = 1.0
		for j = 1; j <= int(gsl.Min(n, mmax)); j++ {
			p_n = p_n * float64(n-j+1)
			result_array[j] = p_n * result_array[j]
		}

		return nil
	}
}

/* Evaluates the series sum_{j=0}^n a_j*He_j(x) with He_j being the j-th probabilists' Hermite polynomial.
 * For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to probabilists' Hermite polynomials is used. */
func Hermite_prob_series_e(n int, x float64, a []float64, result *Result) err.GSLError {
	if n < 0 {
		return DomainError(result)
	} else if n == 0 {
		result.val = a[0]
		result.err = 0.
		return nil
	} else if n == 1 {
		result.val = a[0] + a[1]*x
		result.err = 2. * gsl.Float64Eps * (math.Abs(a[0]) + math.Abs(a[1]*x))
		return nil
	} else {
		/* downward recurrence: b_n = a_n + x b_{n+1} - (n+1) b_{n+2} */
		b0 := 0.
		b1 := 0.
		btmp := 0.

		e0 := 0.
		e1 := 0.
		etmp := e1

		var j int

		for j = n; j >= 0; j-- {
			btmp = b0
			b0 = a[j] + x*b0 - float64(j+1)*b1
			b1 = btmp

			etmp = e0
			e0 = (gsl.Float64Eps*math.Abs(a[j]) + math.Abs(x)*e0 + float64(j+1)*e1)
			e1 = etmp
		}

		result.val = b0
		result.err = e0 + math.Abs(b0)*gsl.Float64Eps
		return nil
	}
}

func Hermite_prob_series(n int, x float64, a []float64) float64 {
	result := new(Result)
	status := Hermite_prob_series_e(n, x, a, result)
	return EvalResult(result, status)
}

/* Evaluates all physicists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
func Hermite_array(nmax int, x float64, result_array []float64) err.GSLError {
	if nmax < 0 {
		return err.ERROR("domain error", err.EDOM)
	} else if nmax == 0 {
		result_array[0] = 1.0
		return nil
	} else if nmax == 1 {
		result_array[0] = 1.0
		result_array[1] = 2.0 * x
		return nil
	} else {
		/* upward recurrence: H_{n+1} = 2x H_n - 2n H_{n-1} */
		var status err.GSLError
		two_x := 2.0 * x
		abs_two_x := math.Abs(two_x)
		thresh1 := gsl.MaxFloat64
		if abs_two_x > 1.0 {
			thresh1 = 0.9 * gsl.MaxFloat64 / abs_two_x
		}
		thresh2 := 0.9 * gsl.MaxFloat64 / 2.0

		p_n0 := 1.0   /* H_0(x) */
		p_n1 := two_x /* H_1(x) */
		p_n := p_n1

		var j int

		result_array[0] = 1.0
		result_array[1] = 2.0 * x

		for j = 1; j < nmax; j++ {
			fj := float64(j)
			if math.Abs(p_n1) > thresh1 || math.Abs(p_n0) > thresh2/fj {
				status = err.Overflow()
			}

			p_n = two_x*p_n1 - 2.0*fj*p_n0
			p_n0 = p_n1
			p_n1 = p_n

			result_array[j+1] = p_n
		}

		return status
	}
}

/* Evaluates the m-th derivative of all physicists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
func Hermite_array_deriv(m, nmax int, x float64, result_array []float64) err.GSLError {
	if nmax < 0 || m < 0 {
		return err.ERROR("domain error", err.EDOM)
	} else if m == 0 {
		Hermite_array(nmax, x, result_array)
		return nil
	} else if nmax < m {
		var j int

		for j = 0; j <= nmax; j++ {
			result_array[j] = 0.0
		}

		return nil
	} else if nmax == m {
		var j int

		for j = 0; j < m; j++ {
			result_array[j] = 0.0
		}

		result_array[nmax] = pow2(m) * Fact(uint(m))

		return nil
	} else if nmax == m+1 {
		var j int

		for j = 0; j < m; j++ {
			result_array[j] = 0.0
		}

		result_array[nmax-1] = pow2(m) * Fact(uint(m))
		result_array[nmax] = result_array[nmax-1] * 2 * float64(m+1) * x
		return nil
	} else {
		/* upward recurrence: H^{(m)}_{n+1} = 2(n+1)/(n-m+1)*(x H^{(m)}_n - n H^{(m)}_{n-1}) */

		p_n0 := pow2(m) * Fact(uint(m))     /* H^{(m)}_{m}(x) */
		p_n1 := p_n0 * 2 * float64(m+1) * x /* H^{(m)}_{m+1}(x) */
		var p_n float64
		var j int

		for j = 0; j < m; j++ {
			result_array[j] = 0.0
		}

		result_array[m] = p_n0
		result_array[m+1] = p_n1

		for j = m + 1; j < nmax; j++ {
			fj := float64(j)
			fm := float64(m)
			p_n = (x*p_n1 - fj*p_n0) * 2 * (fj + 1.0) / (fj - fm + 1.0)
			p_n0 = p_n1
			p_n1 = p_n
			result_array[j+1] = p_n
		}

		return nil
	}
}

/* Evaluates all derivatives (starting from 0) up to the mmax-th derivative of the physicists' Hermite polynomial of order n at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
func Hermite_deriv_array(mmax, n int, x float64, result_array []float64) err.GSLError {
	if n < 0 || mmax < 0 {
		return err.ERROR("domain error", err.EDOM)
	} else if n == 0 {
		var j int

		result_array[0] = 1.0
		for j = 1; j <= mmax; j++ {
			result_array[j] = 0.0
		}

		return nil
	} else if n == 1 && mmax > 0 {
		var j int

		result_array[0] = 2 * x
		result_array[1] = 1.0
		for j = 2; j <= mmax; j++ {
			result_array[j] = 0.0
		}

		return nil
	} else if mmax == 0 {
		result_array[0] = Hermite(n, x)
		return nil
	} else if mmax == 1 {
		result_array[0] = Hermite(n, x)
		result_array[1] = 2 * float64(n) * Hermite(n-1, x)
		return nil
	} else {
		/* upward recurrence */

		k := int(gsl.Max(0, n-mmax))
		p_n0 := Hermite(k, x)   /* H_k(x) */
		p_n1 := Hermite(k+1, x) /* H_{k+1}(x) */
		p_n := p_n1
		var j int

		for j = n + 1; j <= mmax; j++ {
			result_array[j] = 0.0
		}

		result_array[int(gsl.Min(n, mmax))] = p_n0
		result_array[int(gsl.Min(n, mmax)-1)] = p_n1

		for j = int(gsl.Min(mmax, n)) - 1; j > 0; j-- {
			k++
			p_n = 2*x*p_n1 - 2*float64(k)*p_n0
			p_n0 = p_n1
			p_n1 = p_n
			result_array[j-1] = p_n
		}

		p_n = 1.0
		for j = 1; j <= int(gsl.Min(n, mmax)); j++ {
			p_n *= 2.0 * float64(n-j+1)
			result_array[j] *= p_n
		}

		return nil
	}
}

/* Evaluates the series sum_{j=0}^n a_j*H_j(x) with H_j being the j-th physicists' Hermite polynomial.
 * For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to physicists' Hermite polynomials is used. */
func Hermite_series_e(n int, x float64, a []float64, result *Result) err.GSLError {
	if n < 0 {
		return DomainError(result)
	} else if n == 0 {
		result.val = a[0]
		result.err = 0.
		return nil
	} else if n == 1 {
		result.val = a[0] + a[1]*2.*x
		result.err = 2. * gsl.Float64Eps * (math.Abs(a[0]) + math.Abs(a[1]*2.*x))
		return nil
	} else {
		/* downward recurrence: b_n = a_n + 2x b_{n+1} - 2(n+1) b_{n+2} */

		b0 := 0.
		b1 := 0.
		btmp := 0.

		e0 := 0.
		e1 := 0.
		etmp := e1

		var j int

		for j = n; j >= 0; j-- {
			fj := float64(j)
			btmp = b0
			b0 = a[j] + 2.*x*b0 - 2.*(fj+1)*b1
			b1 = btmp

			etmp = e0
			e0 = (gsl.Float64Eps*math.Abs(a[j]) + math.Abs(2.*x)*e0 + 2.*(fj+1)*e1)
			e1 = etmp
		}

		result.val = b0
		result.err = e0 + math.Abs(b0)*gsl.Float64Eps
		return nil
	}
}

func Hermite_series(n int, x float64, a []float64) float64 {
	result := new(Result)
	status := Hermite_series_e(n, x, a, result)
	return EvalResult(result, status)
}

/* Evaluates all Hermite functions up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
func Hermite_func_array(nmax int, x float64, result_array []float64) err.GSLError {
	if nmax < 0 {
		return err.ERROR("domain error", err.EDOM)
	} else if nmax == 0 {
		result_array[0] = math.Exp(-0.5*x*x) / math.Sqrt(math.SqrtPi)
		return nil
	} else if nmax == 1 {
		result_array[0] = math.Exp(-0.5*x*x) / math.Sqrt(math.SqrtPi)
		result_array[1] = result_array[0] * math.Sqrt2 * x
		return nil
	} else {
		/* upward recurrence: Psi_{n+1} = math.Sqrt(2/(n+1))*x Psi_n - math.Sqrt(n/(n+1)) Psi_{n-1} */
		arg := -0.5 * x * x
		hi2 := 1.0 / math.Sqrt(math.SqrtPi)
		hi1 := math.Sqrt2 * x * hi2
		hi := 0.0
		sum_log_scale := 0.0
		var abshi float64
		var i int

		result_array[0] = math.Exp(arg) * hi2
		result_array[1] = result_array[0] * math.Sqrt2 * x

		for i = 2; i <= nmax; i++ {
			fi := float64(i)
			hi = math.Sqrt(2.0/fi)*x*hi1 - math.Sqrt((fi-1.0)/fi)*hi2
			hi2 = hi1
			hi1 = hi

			abshi = math.Abs(hi)
			if abshi > 1.0 {
				log_scale := rnd(math.Log(abshi))
				scale := math.Exp(-log_scale)

				hi *= scale
				hi1 *= scale
				hi2 *= scale
				sum_log_scale += log_scale
			}

			result_array[i] = hi * math.Exp(arg+sum_log_scale)
		}

		return nil
	}
}

/* Evaluates the series sum_{j=0}^n a_j*Psi_j(x) with Psi_j being the j-th Hermite function.
 * For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to Hermite functions is used. */
func Hermite_func_series_e(n int, x float64, a []float64, result *Result) err.GSLError {
	if n < 0 {
		return DomainError(result)
	} else if n == 0 {
		result.val = a[0] * math.Exp(-0.5*x*x) / math.Sqrt(math.SqrtPi)
		result.err = gsl.Float64Eps * math.Abs(result.val)
		return nil
	} else if n == 1 {
		result.val = (a[0] + a[1]*math.Sqrt2*x) * math.Exp(-0.5*x*x) / math.Sqrt(math.SqrtPi)
		result.err = 2. * gsl.Float64Eps * (math.Abs(a[0]) + math.Abs(a[1]*math.Sqrt2*x)) * math.Exp(-0.5*x*x) / math.Sqrt(math.SqrtPi)
		return nil
	} else {
		/* downward recurrence: b_n = a_n + math.Sqrt(2/(n+1))*x b_{n+1} - math.Sqrt((n+1)/(n+2)) b_{n+2} */

		b0 := 0.
		b1 := 0.
		btmp := 0.

		e0 := 0.
		e1 := 0.
		etmp := e1

		var j int

		for j = n; j >= 0; j-- {
			fj := float64(j)
			btmp = b0
			b0 = a[j] + math.Sqrt(2./(fj+1))*x*b0 - math.Sqrt((fj+1.)/(fj+2.))*b1
			b1 = btmp

			etmp = e0
			e0 = (gsl.Float64Eps*math.Abs(a[j]) + math.Sqrt(2./(fj+1))*math.Abs(x)*e0 + math.Sqrt((fj+1.)/(fj+2.))*e1)
			e1 = etmp
		}

		result.val = b0 * math.Exp(-0.5*x*x) / math.Sqrt(math.SqrtPi)
		result.err = e0 + math.Abs(result.val)*gsl.Float64Eps

		return nil
	}
}

func Hermite_func_series(n int, x float64, a []float64) float64 {
	result := new(Result)
	status := Hermite_func_series_e(n, x, a, result)
	return EvalResult(result, status)
}

/* Evaluates the m-th derivative of the Hermite function of order n at position x.
 * A summation including upward recurrences is used. */
func Hermite_func_der_e(m, n int, x float64, result *Result) err.GSLError {
	if m < 0 || n < 0 {
		return DomainError(result)
	} else if m == 0 {
		return Hermite_func_e(n, x, result)
	} else if m == 1 {
		hi2 := 1.0 / math.Sqrt(math.SqrtPi)
		hi1 := math.Sqrt2 * x * hi2
		hi := 0.0
		sum_log_scale := 0.0
		var abshi float64
		var i int
		fn := float64(n)
		for i = 2; i <= n; i++ {
			fi := float64(i)
			hi = math.Sqrt(2.0/fi)*x*hi1 - math.Sqrt((fi-1.0)/fi)*hi2
			hi2 = hi1
			hi1 = hi

			abshi = math.Abs(hi)
			if abshi > 1.0 {
				log_scale := rnd(math.Log(abshi))
				scale := math.Exp(-log_scale)

				hi *= scale
				hi1 *= scale
				hi2 *= scale
				sum_log_scale += log_scale
			}
		}

		/* psi'_n(x) = math.Sqrt(2 n) psi_{n-1} - x psi_n */
		result.val = (math.Sqrt(2.0*fn)*hi2 - x*hi) * math.Exp(-0.5*x*x+sum_log_scale)
		result.err = fn * gsl.Float64Eps * math.Abs(result.val)

		return nil
	} else {
		var j int
		var r, er, b float64
		h0 := 1.
		h1 := x
		eh0 := gsl.Float64Eps
		eh1 := gsl.Float64Eps
		p0 := 1.
		p1 := math.Sqrt2 * x
		ep0 := gsl.Float64Eps
		ep1 := math.Sqrt2 * gsl.Float64Eps
		f := 1.

		for j = int(gsl.Max(1, n-m+1)); j <= n; j++ {
			f *= math.Sqrt(2. * float64(j))
		}

		if m > n {
			if gsl.IsOdd(m - n) {
				f = -f
			}

			for j = 0; j < int(gsl.Min(n, m-n)); j++ {
				f *= float64(m-j) / float64(j+1)
			}
		}

		for j = 1; j <= m-n; j++ {
			b = x*h1 - float64(j)*h0
			h0 = h1
			h1 = b

			b = (math.Abs(x)*eh1 + float64(j)*eh0)
			eh0 = eh1
			eh1 = b
		}

		b = 0.
		for j = 1; j <= n-m; j++ {
			fj := float64(j)
			b = (math.Sqrt2*x*p1 - math.Sqrt(fj)*p0) / math.Sqrt(fj+1.)
			p0 = p1
			p1 = b

			b = (math.Sqrt2*math.Abs(x)*ep1 + math.Sqrt(fj)*ep0) / math.Sqrt(fj+1.)
			ep0 = ep1
			ep1 = b
		}

		b = 0.
		r = 0.
		er = 0.
		for j = int(gsl.Max(0, m-n)); j <= m; j++ {
			r += f * h0 * p0
			er += eh0*math.Abs(f*p0) + ep0*math.Abs(f*h0) + gsl.Float64Eps*math.Abs(f*h0*p0)

			b = x*h1 - float64(j+1)*h0
			h0 = h1
			h1 = b

			b = 0.5 * (math.Abs(x)*eh1 + float64(j+1)*eh0)
			eh0 = eh1
			eh1 = b

			b = (math.Sqrt2*x*p1 - math.Sqrt(float64(n-m+j+1))*p0) / math.Sqrt(float64(n-m+j+2))
			p0 = p1
			p1 = b

			b = 0.5 * (math.Sqrt2*math.Abs(x)*ep1 + math.Sqrt(float64(n-m+j+1))*ep0) / math.Sqrt(float64(n-m+j+2))
			ep0 = ep1
			ep1 = b

			f *= -float64(m-j) / float64(j+1) / math.Sqrt(float64(n-m+j+1)) * gsl.SqrtHalf
		}

		result.val = r * math.Exp(-0.5*x*x) / math.Sqrt(math.SqrtPi)
		result.err = er*math.Abs(math.Exp(-0.5*x*x)/math.Sqrt(math.SqrtPi)) + gsl.Float64Eps*math.Abs(result.val)

		return nil
	}
}

func Hermite_func_der(m, n int, x float64) float64 {
	result := new(Result)
	status := Hermite_func_der_e(m, n, x, result)
	return EvalResult(result, status)
}

func h_zero_init(n, k int) float64 {
	p := 1.
	x := 1.
	y := 1.
	if k == 1 && n > 50 {
		if gsl.IsOdd(n) {
			x = 1. / math.Sqrt(float64(n-1)/6.)
		} else {
			x = 1. / math.Sqrt(0.5*float64(n))
		}
	} else {
		p = -0.7937005259840997373758528196 * Airy_zero_Ai(n/2-k+1)
		x = math.Sqrt(2*float64(n) + 1.)
		y = math.Pow(2*float64(n)+1., 1/6.)
		x = x - p/y - 0.1*p*p/(x*y*y) + (9/280.-p*p*p*11/350.)/(x*x*x) + (p*277/12600.-Pow_int(p, 4)*823/63000.)/Pow_int(x, 4)/y
	}
	p = math.Acos(x / math.Sqrt(2*float64(n)+1.))
	y = gsl.Pi * (-2*(float64(n)/2-float64(k)) - 1.5) / (float64(n) + 0.5)
	if sys.Fcmp(y, math.Sin(2.*p)-2*p, gsl.SqrtFloat64Eps) == 0 { // verify
		return x /* initial approx sufficiently accurate */
	}
	if y > -gsl.Float64Eps {
		return math.Sqrt(2*float64(n) + 1.)
	}
	if p < gsl.Float64Eps {
		p = gsl.Float64Eps
	}
	if p > gsl.PiOver2 {
		p = gsl.PiOver2
	}
	if math.Sin(2.*p)-2*p > y {
		x = math.Max((math.Sin(2.*p)-2*p-y)/4., gsl.SqrtFloat64Eps)
		for ok := true; ok; ok = math.Sin(2.*p)-2*p > y {
			x *= 2.
			p += x
		}
	}
	for ok := true; ok; ok = (sys.Fcmp(x, p, 100*gsl.Float64Eps) != 0) {
		x = p
		p -= (math.Sin(2.*p) - 2.*p - y) / (2.*math.Cos(2.*p) - 2.)
		if p < 0. || p > gsl.PiOver2 {
			p = gsl.PiOver2
		}
	}
	return math.Sqrt(2*float64(n)+1.) * math.Cos(p)
}

/* lookup table for the positive zeros of the probabilists' Hermite polynomials of order 3 through 20 */
var (
	he_zero_tab = []float64{
		1.73205080756887729352744634151,
		0.741963784302725857648513596726,
		2.33441421833897723931751226721,
		1.35562617997426586583052129087,
		2.85697001387280565416230426401,
		0.616706590192594152193686099399,
		1.88917587775371067550566789858,
		3.32425743355211895236183546247,
		1.154405394739968127239597758838,
		2.36675941073454128861885646856,
		3.75043971772574225630392202571,
		0.539079811351375108072461918694,
		1.63651904243510799922544657297,
		2.80248586128754169911301080618,
		4.14454718612589433206019783917,
		1.023255663789132524828148225810,
		2.07684797867783010652215614374,
		3.20542900285646994336567590292,
		4.51274586339978266756667884317,
		0.484935707515497653046233483105,
		1.46598909439115818325066466416,
		2.48432584163895458087625118368,
		3.58182348355192692277623675546,
		4.85946282833231215015516494660,
		0.928868997381063940144111999584,
		1.87603502015484584534137013967,
		2.86512316064364499771968407254,
		3.93616660712997692868589612142,
		5.18800122437487094818666404539,
		0.444403001944138945299732445510,
		1.34037519715161672153112945211,
		2.25946445100079912386492979448,
		3.22370982877009747166319001956,
		4.27182584793228172295999293076,
		5.50090170446774760081221630899,
		0.856679493519450033897376121795,
		1.72541837958823916151095838741,
		2.62068997343221478063807762201,
		3.56344438028163409162493844661,
		4.59139844893652062705231872720,
		5.80016725238650030586450565322,
		0.412590457954601838167454145167,
		1.24268895548546417895063983219,
		2.08834474570194417097139675101,
		2.96303657983866750254927123447,
		3.88692457505976938384755016476,
		4.89693639734556468372449782879,
		6.08740954690129132226890147034,
		0.799129068324547999424888414207,
		1.60671006902872973652322479373,
		2.43243682700975804116311571682,
		3.28908242439876638890856229770,
		4.19620771126901565957404160583,
		5.19009359130478119946445431715,
		6.36394788882983831771116094427,
		0.386760604500557347721047189801,
		1.16382910055496477419336819907,
		1.95198034571633346449212362880,
		2.76024504763070161684598142269,
		3.60087362417154828824902745506,
		4.49295530252001124266582263095,
		5.47222570594934308841242925805,
		6.63087819839312848022981922233,
		0.751842600703896170737870774614,
		1.50988330779674075905491513417,
		2.28101944025298889535537879396,
		3.07379717532819355851658337833,
		3.90006571719800990903311840097,
		4.77853158962998382710540812497,
		5.74446007865940618125547815768,
		6.88912243989533223256205432938,
		0.365245755507697595916901619097,
		1.09839551809150122773848360538,
		1.83977992150864548966395498992,
		2.59583368891124032910545091458,
		3.37473653577809099529779309480,
		4.18802023162940370448450911428,
		5.05407268544273984538327527397,
		6.00774591135959752029303858752,
		7.13946484914647887560975631213,
		0.712085044042379940413609979021,
		1.42887667607837287134157901452,
		2.15550276131693514033871248449,
		2.89805127651575312007902775275,
		3.66441654745063847665304033851,
		4.46587262683103133615452574019,
		5.32053637733603803162823765939,
		6.26289115651325170419416064557,
		7.38257902403043186766326977122,
		0.346964157081355927973322447164,
		1.04294534880275103146136681143,
		1.74524732081412671493067861704,
		2.45866361117236775131735057433,
		3.18901481655338941485371744116,
		3.94396735065731626033176813604,
		4.73458133404605534390170946748,
		5.57873880589320115268040332802,
		6.51059015701365448636289263918,
		7.61904854167975829138128156060,
	}
)

/*
 * Computes the s-th zero the probabilists' Hermite polynomial of order n.
 * A Newton iteration using a continued fraction representation adapted from:
 *
 * [E.T. Whittaker (1914), On the continued fractions which represent the
 * functions of Hermite and other functions defined by differential equations,
 * Proceedings of the Edinburgh Mathematical Society, 32, 65-74]
 *
 * is performed with the initial approximation from
 *
 * [Arpad Elbert and Martin E. Muldoon, Approximations for zeros of Hermite
 * functions, pp. 117-126 in D. Dominici and R. S. Maier, eds, "Special Functions
 * and Orthogonal Polynomials", Contemporary Mathematics, vol 471 (2008)]
 *
 * refined via the bisection method.
 */

func Hermite_prob_zero_e(n, s int, result *Result) err.GSLError {
	if n <= 0 || s < 0 || s > n/2 {
		return DomainError(result)
	} else if s == 0 {
		if gsl.IsOdd(n) {
			result.val = 0.
			result.err = 0.
			return nil
		} else {
			return DomainError(result)
		}
	} else if n == 2 {
		result.val = 1.
		result.err = 0.
		return nil
	} else if n < 21 {
		var idx int
		if gsl.IsOdd(n) {
			idx = n / 2
		}

		idx += ((n / 2) * (n/2 - 1)) + s - 2

		result.val = he_zero_tab[idx]
		result.err = gsl.Float64Eps * (result.val)
		return nil
	} else {
		d := 1.0
		x := 1.0
		x0 := 1.0
		var j int
		x = h_zero_init(n, s) * math.Sqrt2
		for ok := true; ok; ok = sys.Fcmp(x, x0, 10*gsl.Float64Eps) != 0 {
			x0 = x
			d = 0.0
			for j = 1; j < n; j++ {
				d = float64(j) / (x - d)
			}
			x -= (x - d) / float64(n)
			/* gsl_fcmp can be used since the smallest zero approaches 1/sqrt(n) or 1/sqrt((n-1)/3.) for large n and thus all zeros are non-zero (except for the trivial case handled above) */
		}
		result.val = x
		result.err = 2*gsl.Float64Eps*x + math.Abs(x-x0)
		return nil
	}
}
func Hermite_prob_zero(n, s int) float64 {
	result := new(Result)
	status := Hermite_prob_zero_e(n, s, result)
	return EvalResult(result, status)
}

/* lookup table for the positive zeros of the physicists' Hermite polynomials of order 3 through 20 */
var (
	h_zero_tab = []float64{
		1.22474487139158904909864203735,
		0.524647623275290317884060253835,
		1.65068012388578455588334111112,
		0.958572464613818507112770593893,
		2.02018287045608563292872408814,
		0.436077411927616508679215948251,
		1.335849074013696949714895282970,
		2.35060497367449222283392198706,
		0.816287882858964663038710959027,
		1.67355162876747144503180139830,
		2.65196135683523349244708200652,
		0.381186990207322116854718885584,
		1.157193712446780194720765779063,
		1.98165675669584292585463063977,
		2.93063742025724401922350270524,
		0.723551018752837573322639864579,
		1.46855328921666793166701573925,
		2.26658058453184311180209693284,
		3.19099320178152760723004779538,
		0.342901327223704608789165025557,
		1.03661082978951365417749191676,
		1.75668364929988177345140122011,
		2.53273167423278979640896079775,
		3.43615911883773760332672549432,
		0.656809566882099765024611575383,
		1.32655708449493285594973473558,
		2.02594801582575533516591283121,
		2.78329009978165177083671870152,
		3.66847084655958251845837146485,
		0.314240376254359111276611634095,
		0.947788391240163743704578131060,
		1.59768263515260479670966277090,
		2.27950708050105990018772856942,
		3.02063702512088977171067937518,
		3.88972489786978191927164274724,
		0.605763879171060113080537108602,
		1.22005503659074842622205526637,
		1.85310765160151214200350644316,
		2.51973568567823788343040913628,
		3.24660897837240998812205115236,
		4.10133759617863964117891508007,
		0.291745510672562078446113075799,
		0.878713787329399416114679311861,
		1.47668273114114087058350654421,
		2.09518325850771681573497272630,
		2.74847072498540256862499852415,
		3.46265693360227055020891736115,
		4.30444857047363181262129810037,
		0.565069583255575748526020337198,
		1.13611558521092066631913490556,
		1.71999257518648893241583152515,
		2.32573248617385774545404479449,
		2.96716692790560324848896036355,
		3.66995037340445253472922383312,
		4.49999070730939155366438053053,
		0.273481046138152452158280401965,
		0.822951449144655892582454496734,
		1.38025853919888079637208966969,
		1.95178799091625397743465541496,
		2.54620215784748136215932870545,
		3.17699916197995602681399455926,
		3.86944790486012269871942409801,
		4.68873893930581836468849864875,
		0.531633001342654731349086553718,
		1.06764872574345055363045773799,
		1.61292431422123133311288254454,
		2.17350282666662081927537907149,
		2.75776291570388873092640349574,
		3.37893209114149408338327069289,
		4.06194667587547430689245559698,
		4.87134519367440308834927655662,
		0.258267750519096759258116098711,
		0.776682919267411661316659462284,
		1.30092085838961736566626555439,
		1.83553160426162889225383944409,
		2.38629908916668600026459301424,
		2.96137750553160684477863254906,
		3.57376906848626607950067599377,
		4.24811787356812646302342016090,
		5.04836400887446676837203757885,
		0.503520163423888209373811765050,
		1.01036838713431135136859873726,
		1.52417061939353303183354859367,
		2.04923170985061937575050838669,
		2.59113378979454256492128084112,
		3.15784881834760228184318034120,
		3.76218735196402009751489394104,
		4.42853280660377943723498532226,
		5.22027169053748216460967142500,
		0.245340708300901249903836530634,
		0.737473728545394358705605144252,
		1.23407621539532300788581834696,
		1.73853771211658620678086566214,
		2.25497400208927552308233334473,
		2.78880605842813048052503375640,
		3.34785456738321632691492452300,
		3.94476404011562521037562880052,
		4.60368244955074427307767524898,
		5.38748089001123286201690041068,
	}
)

/*
 * Computes the s-th zero the physicists' Hermite polynomial of order n, thus also
 * the s-th zero of the Hermite function of order n. A Newton iteration using a continued
 * fraction representation adapted from:
 *
 * [E.T. Whittaker (1914), On the continued fractions which represent the functions of Hermite
 * and other functions defined by differential equations, Proceedings of the Edinburgh Mathematical
 * Society, 32, 65-74]
 *
 * An initial approximation is used from:
 *
 * [Arpad Elbert and Martin E. Muldoon, Approximations for zeros of Hermite functions,
 * pp. 117-126 in D. Dominici and R. S. Maier, eds, "Special Functions and Orthogonal Polynomials",
 * Contemporary Mathematics, vol 471 (2008)]
 *
 * which is refined via the bisection method.
 */

func Hermite_zero_e(n, s int, result *Result) err.GSLError {
	if n <= 0 || s < 0 || s > n/2 {
		return DomainError(result)
	} else if s == 0 {
		if gsl.IsOdd(n) {
			result.val = 0.
			result.err = 0.
			return nil
		} else {
			return DomainError(result)
		}
	} else if n == 2 {
		result.val = gsl.SqrtHalf
		result.err = 0.
		return nil
	} else if n < 21 {
		var idx int
		if gsl.IsOdd(n) {
			idx = n / 2
		}
		idx += ((n / 2) * (n/2 - 1)) + s - 2
		result.val = h_zero_tab[idx]
		result.err = gsl.Float64Eps * (result.val)
		return nil
	} else {
		d := 1.
		x := 1.
		x0 := 1.
		var j int

		x = h_zero_init(n, s)
		for ok := true; ok; ok = sys.Fcmp(x, x0, 10*gsl.Float64Eps) != 0 {
			x0 = x
			d = 0.

			for j = 1; j < n; j++ {
				d = 2 * float64(j) / (2.*x - d)
			}

			x -= (2*x - d) * 0.5 / float64(n)

			/* gsl_fcmp can be used since the smallest zero approaches 1/math.Sqrt(n) or 1/math.Sqrt((n-1)/3.)
			 * for large n and thus all zeros are non-zero (except for the trivial case handled above) */
		}

		result.val = x
		result.err = 2*gsl.Float64Eps*x + math.Abs(x-x0)

		return nil
	}
}

func Hermite_zero(n, s int) float64 {
	result := new(Result)
	status := Hermite_zero_e(n, s, result)
	return EvalResult(result, status)
}

func Hermite_func_zero_e(n, s int, result *Result) err.GSLError {
	return Hermite_zero_e(n, s, result)
}

func Hermite_func_zero(n, s int) float64 {
	result := new(Result)
	status := Hermite_func_zero_e(n, s, result)
	return EvalResult(result, status)
}

func Hermite_phys_e(n int, x float64, result *Result) err.GSLError {
	return Hermite_e(n, x, result)
}

func Hermite_phys(n int, x float64) float64 {
	result := new(Result)
	status := Hermite_phys_e(n, x, result)
	return EvalResult(result, status)
}

func Hermite_phys_der_e(m, n int, x float64, result *Result) err.GSLError {
	return Hermite_deriv_e(m, n, x, result)
}

func Hermite_phys_der(m, n int, x float64) float64 {
	result := new(Result)
	status := Hermite_phys_der_e(m, n, x, result)
	return EvalResult(result, status)
}

func Hermite_phys_array(nmax int, x float64, result_array []float64) err.GSLError {
	return Hermite_array(nmax, x, result_array)
}

func Hermite_phys_series_e(n int, x float64, a []float64, result *Result) err.GSLError {
	return Hermite_series_e(n, x, a, result)
}

func Hermite_phys_series(n int, x float64, a []float64) float64 {
	result := new(Result)
	status := Hermite_phys_series_e(n, x, a, result)
	return EvalResult(result, status)
}

func Hermite_phys_array_der(m, nmax int, x float64, result_array []float64) err.GSLError {
	return Hermite_array_deriv(m, nmax, x, result_array)
}

func Hermite_phys_der_array(mmax, n int, x float64, result_array []float64) err.GSLError {
	return Hermite_deriv_array(mmax, n, x, result_array)
}

func Hermite_phys_zero_e(n, s int, result *Result) err.GSLError {
	return Hermite_zero_e(n, s, result)
}

func Hermite_phys_zero(n, s int) float64 {
	result := new(Result)
	status := Hermite_phys_zero_e(n, s, result)
	return EvalResult(result, status)
}

func Hermite_prob_array_der(m, nmax int, x float64, result_array []float64) err.GSLError {
	return Hermite_prob_array_deriv(m, nmax, x, result_array)
}

func Hermite_prob_der_array(mmax, n int, x float64, result_array []float64) err.GSLError {
	return Hermite_prob_deriv_array(mmax, n, x, result_array)
}

func Hermite_prob_der_e(m, n int, x float64, result *Result) err.GSLError {
	return Hermite_prob_deriv_e(m, n, x, result)
}

func Hermite_prob_der(m, n int, x float64) float64 {
	result := new(Result)
	status := Hermite_prob_deriv_e(m, n, x, result)
	return EvalResult(result, status)
}
