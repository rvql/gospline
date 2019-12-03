package gospline

// solves diagonal matrix using Thomas algorithm
// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
// will change input slices
// the matrix should have full rank
func triThomas(a, b, c, d []float64) []float64 {
	n := len(b)
	if len(a)+1 != n || n != len(c)+1 || n != len(d) {
		panic("invalid input slices")
	}

	d[0] /= b[0]
	if n == 1 {
		return d
	}

	c[0] /= b[0]
	for i := 1; i < n-1; i++ {
		div := b[i] - a[i-1]*c[i-1]
		c[i] /= div
		d[i] = (d[i] - a[i-1]*d[i-1]) / div
	}
	d[n-1] = (d[n-1] - a[n-2]*d[n-2]) / (b[n-1] - a[n-2]*c[n-2])
	for i := n - 2; i >= 0; i-- {
		d[i] -= c[i] * d[i+1]
	}
	return d
}
