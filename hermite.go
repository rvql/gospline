package gospline

// See https://en.wikipedia.org/wiki/Monotone_cubic_interpolation

import (
	"math"
	"sort"
)

type hermite struct {
	x []float64
	p []float64
	m []float64
	n int

	segs []*hermiteSegment
}

// p(x) = a(x - xk)^3 + b(x - xk)^2 + c(x - xk) + d
type hermiteSegment struct {
	a float64
	b float64
	c float64
	d float64
}

// NewMonotoneSpline returns cubic monotone hermite spline
// the whole array must be in ascending order
func NewMonotoneSpline(x, y []float64) Spline {
	if len(x) != len(y) {
		panic("array length mismatch")
	}
	n := len(x)
	if !sort.Float64sAreSorted(x) {
		panic("values in x must be in ascending order")
	}

	// 1. Compute the slopes of the secant lines between successive points.
	dk := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		dk[i] = (y[i+1] - y[i]) / (x[i+1] - x[i])
	}

	// 2. Initialize interior tangents as the average of the secants.
	m := make([]float64, n)
	for i := 1; i < n-1; i++ {
		if dk[i-1] < 0 && dk[i] > 0 || dk[i-1] > 0 && dk[i] < 0 {
			m[i] = 0
		} else {
			m[i] = (dk[i-1] + dk[i]) / 2
		}
	}
	m[0] = dk[0]
	m[n-1] = dk[n-2]

	for i := 0; i < n-1; i++ {
		// 3. If dk[i] == 0, spline connecting these points
		// must be flat to preserve monotonicity.
		if y[i] == y[i+1] {
			m[i] = 0
			m[i+1] = 0
			continue
		}
		// 4. Check for strict monotonicity.
		alpha := m[i] / dk[i]
		if i >= 1 {
			prevBeta := m[i] / dk[i-1]
			if alpha < 0 && prevBeta < 0 {
				m[i] = 0
			}
		}
		beta := m[i+1] / dk[i]
		// 5. Prevent overshoot.
		if alpha*alpha+beta*beta > 9 {
			tao := 3 / math.Hypot(alpha, beta)
			m[i] = tao * alpha * dk[i]
			m[i+1] = tao * beta * dk[i]
		}
	}

	return newHermite(x, y, m)
}

func (hm *hermite) At(x float64) float64 {
	if x <= hm.x[0] {
		return hm.p[0]
	}
	if x >= hm.x[hm.n-1] {
		return hm.p[hm.n-1]
	}

	// Find the greatest numbered segment with hm.x > x, then
	// choose the one before it.  If x is one of the input points,
	// this selects the segment such that hm.x[seg] == x.
	seg := sort.SearchFloat64s(hm.x, x)
	if hm.x[seg] != x {
		seg--
	}

	if hm.segs == nil {
		hm.segs = make([]*hermiteSegment, hm.n-1)
	}

	s := hm.segs[seg]
	if s == nil {
		h := hm.x[seg+1] - hm.x[seg]
		s = &hermiteSegment{
			a: ((2*hm.p[seg]-2*hm.p[seg+1])/h + hm.m[seg] + hm.m[seg+1]) / h / h,
			b: ((-3*hm.p[seg]+3*hm.p[seg+1])/h - 2*hm.m[seg] - hm.m[seg+1]) / h,
			c: hm.m[seg],
			d: hm.p[seg],
		}
		hm.segs[seg] = s
	}

	dx := x - hm.x[seg]

	return dx*(dx*(dx*s.a+s.b)+s.c) + s.d
}

func (hm *hermite) Range(start, end, step float64) []float64 {
	if start > end {
		panic("start must be smaller than end")
	}
	n := int((end-start)/step) + 1
	v := make([]float64, n)
	x := start
	for i := 0; i < n; i++ {
		v[i] = hm.At(x)
		x += step
	}
	return v
}

func newHermite(x, p, m []float64) Spline {
	if len(x) != len(p) || len(p) != len(m) {
		panic("array length mismatch")
	}
	n := len(x)
	xx := make([]float64, n)
	copy(xx, x)
	pp := make([]float64, n)
	copy(pp, p)
	mm := make([]float64, n)
	copy(mm, m)
	return &hermite{
		x: xx,
		p: pp,
		m: mm,
		n: n,
	}
}
