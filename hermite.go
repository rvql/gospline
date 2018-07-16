package gospline

import (
	"fmt"
	"math"
)

type hermite struct {
	x []float64
	p []float64
	m []float64
	n int

	segs []hermiteSegment
}

// p(x) = a(x - xk)^3 + b(x - xk)^2 + c(x - xk) + d
type hermiteSegment struct {
	ready bool
	a     float64
	b     float64
	c     float64
	d     float64
}

// NewMonotoneSpline returns cubic monotone hermite spline
// the whole array must be in ascending order
func NewMonotoneSpline(x, y []float64) Spline {
	if len(x) != len(y) {
		panic("array length mismatch")
	}
	n := len(x)
	for i := 0; i < n; i++ {
		if i < n-1 && x[i] >= x[i+1] {
			panic("values in x must be in ascending order")
		}
	}

	dk := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		dk[i] = (y[i+1] - y[i]) / (x[i+1] - x[i])
	}

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
		// TODO: maybe fix this equality test of float64
		if y[i] == y[i+1] {
			m[i] = 0
			m[i+1] = 0
			continue
		}
		alpha := m[i] / dk[i]
		if i >= 1 {
			prevBeta := m[i] / dk[i-1]
			if alpha < 0 && prevBeta < 0 {
				m[i] = 0
			}
		}
		beta := m[i+1] / dk[i]
		if alpha*alpha+beta*beta > 9 {
			tao := 3 / math.Hypot(alpha, beta)
			m[i] = tao * alpha * dk[i]
			m[i+1] = tao * beta * dk[i]
		}
	}

	return newHermite(x, y, m)
}

func (hm *hermite) At(x float64) float64 {
	if x < hm.x[0] || x > hm.x[hm.n-1] {
		panic(fmt.Sprintf("x %g should be between %g and %g", x, hm.x[0], hm.x[hm.n-1]))
	}

	var seg int
	if x == hm.x[hm.n-1] {
		seg = hm.n - 1
	} else {
		for seg = 0; seg < hm.n-1; seg++ {
			if hm.x[seg] <= x && x < hm.x[seg+1] {
				break
			}
		}
	}

	if hm.segs == nil {
		hm.segs = make([]hermiteSegment, hm.n-1)
	}

	s := &hm.segs[seg]
	if !s.ready {
		h := hm.x[seg+1] - hm.x[seg]
		s.a = ((2*hm.p[seg]-2*hm.p[seg+1])/h + hm.m[seg] + hm.m[seg+1]) / h / h
		s.b = ((-3*hm.p[seg]+3*hm.p[seg+1])/h - 2*hm.m[seg] - hm.m[seg+1]) / h
		s.c = hm.m[seg]
		s.d = hm.p[seg] / h
		s.ready = true
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
	fmt.Println(x, p, m)
	n := len(x)
	for i := 0; i < n; i++ {
		if i < n-1 && x[i] >= x[i+1] {
			panic("values in x must be in ascending order")
		}
	}
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
