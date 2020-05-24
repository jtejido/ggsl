package err

import (
	"runtime"
)

const (
	CONTINUE = iota - 2
	FAILURE
	SUCCESS // sanity purpose
	EDOM
	ERANGE
	EFAULT
	EINVAL
	EFAILED
	EFACTOR
	ESANITY
	ENOMEM
	EBADFUNC
	ERUNAWAY
	EMAXITER
	EZERODIV
	EBADTOL
	ETOL
	EUNDRFLW
	EOVRFLW
	ELOSS
	EROUND
	EBADLEN
	ENOTSQR
	ESING
	EDIVERGE
	EUNSUP
	EUNIMPL
	ECACHE
	ETABLE
	ENOPROG
	ENOPROGJ
	ETOLF
	ETOLX
	ETOLG
	EOF
)

const (
	continue_msg = "Iteration has not converged"
	failure_msg  = "General failure"
	edom_msg     = "Input domain error, e.g sqrt(-1)"
	erange_msg   = "Output range error, e.g. exp(1e100)"
	efault_msg   = "Invalid pointer"
	einval_msg   = "Invalid argument supplied by user"
	efailed_msg  = "Generic failure"
	efactor_msg  = "Factorization failed"
	esanity_msg  = "Sanity check failed - shouldn't happen"
	enomem_msg   = "Malloc failed"
	ebadfunc_msg = "Problem with user-supplied function"
	erunaway_msg = "Iterative process is out of control"
	emaxiter_msg = "Exceeded max number of iterations"
	ezerodiv_msg = "Tried to divide by zero"
	ebadtol_msg  = "User specified an invalid tolerance"
	etol_msg     = "Failed to reach the specified tolerance"
	eundrflw_msg = "Underflow"
	eovrflw_msg  = "Overflow"
	eloss_msg    = "Loss of accuracy"
	eround_msg   = "Failed because of roundoff error"
	ebadlen_msg  = "Matrix, vector lengths are not conformant"
	enotsqr_msg  = "Matrix not square"
	esing_msg    = "Apparent singularity detected"
	ediverge_msg = "Integral or series is divergent"
	eunsup_msg   = "Requested feature is not supported by the hardware"
	eunimpl_msg  = "Requested feature not (yet) implemented"
	ecache_msg   = "Cache limit exceeded"
	etable_msg   = "Table limit exceeded"
	enoprog_msg  = "Iteration is not making progress towards solution"
	enoprogj_msg = "Jacobian evaluations are not improving the solution"
	etolf_msg    = "Cannot reach the specified tolerance in F"
	etolx_msg    = "Cannot reach the specified tolerance in X"
	etolg_msg    = "Cannot reach the specified tolerance in gradient"
	eof_msg      = "End of file"
)

/* GSL_ERROR: call the error handler, and return the error */
func ERROR(reason string, gsl_errno int) GSLError {
	_, file, line, _ := runtime.Caller(1)
	Error(reason, file, line, gsl_errno)
	return New(gsl_errno, reason)
}

/* GSL_ERROR_VAL: call the error handler, and return the given value */
func ERROR_VAL(reason string, gsl_errno int, value float64) float64 {
	_, file, line, _ := runtime.Caller(1)
	Error(reason, file, line, gsl_errno)
	return value
}

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */
func ERROR_VOID(reason string, gsl_errno int) {
	_, file, line, _ := runtime.Caller(1)
	Error(reason, file, line, gsl_errno)
	return
}

/* GSL_ERROR_NULL suitable for out-of-memory conditions */
func ERROR_NULL(reason string, gsl_errno int) {
	ERROR_VAL(reason, gsl_errno, 0)
}

// If multiple errors, select first non-Success one
func ErrorSelect(err ...GSLError) GSLError {
	for _, e := range err {
		if e != nil {
			return e
		}
	}

	return nil
}

// helpers
func Continue() GSLError {
	return New(CONTINUE, continue_msg)
}

func Failure() GSLError {
	return New(FAILURE, failure_msg)
}

func Success() GSLError {
	return GSLError(nil)
}

func Domain() GSLError {
	return New(EDOM, edom_msg)
}

func Range() GSLError {
	return New(ERANGE, erange_msg)
}

func Fault() GSLError {
	return New(EFAULT, efault_msg)
}

func Invalid() GSLError {
	return New(EINVAL, einval_msg)
}

func Generic() GSLError {
	return New(EFAILED, efailed_msg)
}

func Factor() GSLError {
	return New(EFACTOR, efactor_msg)
}

func Sanity() GSLError {
	return New(ESANITY, esanity_msg)
}

func NoMemory() GSLError {
	return New(ENOMEM, enomem_msg)
}

func BadFunc() GSLError {
	return New(EBADFUNC, ebadfunc_msg)
}

func RunAway() GSLError {
	return New(ERUNAWAY, erunaway_msg)
}

func MaxIteration() GSLError {
	return New(EMAXITER, emaxiter_msg)
}

func ZeroDiv() GSLError {
	return New(EZERODIV, ezerodiv_msg)
}

func BadTolerance() GSLError {
	return New(EBADTOL, ebadtol_msg)
}

func Tolerance() GSLError {
	return New(ETOL, etol_msg)
}

func Underflow() GSLError {
	return New(EUNDRFLW, eundrflw_msg)
}

func Overflow() GSLError {
	return New(EOVRFLW, eovrflw_msg)
}

func Loss() GSLError {
	return New(ELOSS, eloss_msg)
}

func Round() GSLError {
	return New(EROUND, eround_msg)
}

func BadLength() GSLError {
	return New(EBADLEN, ebadlen_msg)
}

func NotSquare() GSLError {
	return New(ENOTSQR, enotsqr_msg)
}

func Singularity() GSLError {
	return New(ESING, esing_msg)
}

func Divergent() GSLError {
	return New(EDIVERGE, ediverge_msg)
}

func NotSupported() GSLError {
	return New(EUNSUP, eunsup_msg)
}

func NotImplemented() GSLError {
	return New(EUNIMPL, eunimpl_msg)
}

func CacheExceeded() GSLError {
	return New(ECACHE, ecache_msg)
}

func TableExceeded() GSLError {
	return New(ETABLE, etable_msg)
}

func NoProgress() GSLError {
	return New(ENOPROG, enoprog_msg)
}

func NoProgressJacobian() GSLError {
	return New(ENOPROGJ, enoprogj_msg)
}

func ToleranceF() GSLError {
	return New(ETOLF, etolf_msg)
}

func ToleranceX() GSLError {
	return New(ETOLX, etolx_msg)
}

func ToleranceGradient() GSLError {
	return New(ETOLG, etolg_msg)
}

func EndOfFile() GSLError {
	return New(EOF, eof_msg)
}
