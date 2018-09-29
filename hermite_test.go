package gospline

import (
	"testing"
)

func TestHermiteSpline(t *testing.T) {
	const epsilon = 1e-9
	x := [...]float64{0, 0.16, 0.42, 0.6425, 0.8575}
	y := [...]float64{0, 32, 237, 255, 0}
	s := NewMonotoneSpline(x[:], y[:])

	for i := range x {
		if s.At(x[i]) != y[i] {
			t.Errorf("interpolate incorrect at %f", x[i])
		}
	}
}
