package err

import (
	"fmt"
	"github.com/lucky-se7en/ggsl/test"
	"strconv"
	"testing"
)

const (
	MAX_ERRS = 64
)

type error struct {
	number int
	name   string
}

var (
	i, j, n int
	errors  = make([]error, MAX_ERRS)
)

func check(x int) {
	errors[n].number = x
	errors[n].name = strconv.Itoa(x)
	n++
}

func TestError(t *testing.T) {
	check(SUCCESS)
	check(FAILURE)
	check(CONTINUE)
	check(EDOM)
	check(ERANGE)
	check(EFAULT)
	check(EINVAL)
	check(EFAILED)
	check(EFACTOR)
	check(ESANITY)
	check(ENOMEM)
	check(EBADFUNC)
	check(ERUNAWAY)
	check(EMAXITER)
	check(EZERODIV)
	check(EBADTOL)
	check(ETOL)
	check(EUNDRFLW)
	check(EOVRFLW)
	check(ELOSS)
	check(EROUND)
	check(EBADLEN)
	check(ENOTSQR)
	check(ESING)
	check(EDIVERGE)
	check(EUNSUP)
	check(EUNIMPL)
	check(ECACHE)
	check(ETABLE)
	check(ENOPROG)
	check(ENOPROGJ)
	check(ETOLF)
	check(ETOLX)
	check(ETOLG)
	check(EOF)

	for i := 0; i < n; i++ {
		fmt.Printf("%s = %d\n", errors[i].name, errors[i].number)
	}

	for i := 0; i < n; i++ {
		status := 0
		for j := 0; j < n; j++ {
			if j != i {
				if errors[i].number == errors[j].number {
					status |= 1
				}
			}
		}

		test.Test(t, status, "%s is distinct from other error values")
	}

}
