package gospline

import (
	"testing"
)

func TestHermiteSpline(t *testing.T) {
	x := []float64{0, 0.16, 0.42, 0.6425, 0.8575}
	y := []float64{0, 32, 237, 255, 0}
	s := NewMonotoneSpline(x, y)

	// Test the input points are mapped exactly.
	for i := range x {
		if at := s.At(x[i]); !floatEquals(at, y[i]) {
			t.Errorf("interpolate incorrect at %f: %f != %f", x[i], at, y[i])
		}
	}

	const intermediate = 10000

	// Test that intermediate points are monotonic.
	for i := range x[:len(x)-1] {
		diff := y[i+1] - y[i]
		ypos := y[i]
		for j := 0; j < intermediate; j++ {
			xval := x[i] + (x[i+1]-x[i])*float64(j)/intermediate
			yval := s.At(xval)
			jdiff := yval - ypos

			if jdiff*diff < 0 {
				t.Errorf("not monotone at x=%f y=%f", xval, yval)
			}
			ypos = yval
		}
	}
}
