package poly

func PolyEval(c []float64, len int, x float64) float64 {
	ans := c[len-1]
	for i := len - 1; i > 0; i-- {
		ans = c[i-1] + x*ans
	}

	return ans
}
