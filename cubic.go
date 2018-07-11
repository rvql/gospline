package gospline

type boundary uint

// bounday types
const (
	CubicFirstDeriv boundary = iota
	CubicSecondDeriv
	CubicThirdDeriv
	CubicPeriodic
)

type cubic struct {
	x []float64
	y []float64
	n int

	boundary
	f0 float64
	fn float64

	m []float64
}

// NewCubicSpline returns cubic spline with natural boundary
// the values in x must be in [0.0, 1.0], and the whole array must be in ascending order
func NewCubicSpline(x, y []float64) Spline {
	return newSpline(x, y, CubicSecondDeriv, 0, 0)
}

func (c *cubic) At(x float64) float64 {
	return 0
}

func (c *cubic) Range(start, end, step float64) []float64 {
	return make([]float64, 0)
}

func newSpline(x, y []float64, b boundary, f0, fn float64) Spline {
	if len(x) != len(y) {
		panic("array length mismatch")
	}
	n := len(x)
	for i := 0; i < n; i++ {
		if x[i] < 0 || x[i] > 1 {
			panic("values in x must be in [0.0, 1.0]")
		}
		if i < n-1 && x[i] >= x[i+1] {
			panic("values in x must be in ascending order")
		}
	}
	c := cubic{
		x:        x,
		y:        y,
		n:        n,
		boundary: b,
		f0:       f0,
		fn:       fn,
	}
	return &c
}

func (c *cubic) calculateM() {
	h := make([]float64, c.n)
	for i := 1; i < c.n; i++ {
		h[i] = c.x[i] - c.x[i-1]
	}

	mu := make([]float64, c.n)
	lambda := make([]float64, c.n)
	diag := make([]float64, c.n)
	d := make([]float64, c.n)
	for i := 1; i < c.n-1; i++ {
		mu[i] = h[i] / (h[i] + h[i+1])
		lambda[i] = 1 - mu[i]
		diag[i] = 2
		d[i] = 6 * (c.y[i-1]/h[i]/(h[i]+h[i+1]) - c.y[i]/h[i]/h[i+1] + c.y[i]/(h[i]+h[i+1])/h[i+1])
	}
	diag[0] = 2
	diag[c.n-1] = 2

	// boundary
	switch c.boundary {
	case CubicFirstDeriv:
		mu[c.n-1] = 1
		lambda[0] = 1
		d[0] = 6 * ((c.y[1]-c.y[0])/h[1] - c.f0) / h[1]
		d[c.n-1] = 6 * (c.fn - (c.y[c.n-1]-c.y[c.n-2])/h[c.n-1]) / h[c.n-1]
	case CubicSecondDeriv:
		// lambda[0] == mu[n-1] == 0
		d[0] = 2 * c.f0
		d[c.n-1] = 2 * c.fn
	case CubicThirdDeriv:
		fallthrough
	case CubicPeriodic:
		panic("not yet implemented")
	}

	c.m = triThomas(mu[1:], diag, lambda[:c.n], d)
}
