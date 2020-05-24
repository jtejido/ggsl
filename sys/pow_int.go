package sys

func PowInt(x float64, n int) float64 {
	var un uint
	if n < 0 {
		x = 1.0 / x
		un = uint(-n)
	} else {
		un = uint(n)
	}

	return PowUint(x, un)
}

func PowUint(x float64, n uint) float64 {
	value := 1.0

	/* repeated squaring method
	 * returns 0.0^0 = 1.0, so continuous in x
	 */
	for n > 0 {
		if n&1 == 1 {
			value *= x /* for n odd */
		}
		n >>= 1
		x *= x
	}

	return value
}
