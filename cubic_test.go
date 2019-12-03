package gospline

import (
	"testing"
)

func TestSimpleCubicSpline(t *testing.T) {
	x := [...]float64{0, 1, 2, 3}
	y := [...]float64{0, 0.5, 2, 1.5}
	s := newSpline(x[:], y[:], CubicSecondDeriv, 0, 0)
	for i := range x {
		if !floatEquals(s.At(x[i]), y[i]) {
			t.Errorf("expected f(%g) = %g, but the result is %g", x[i], y[i], s.At(x[i]))
		}
	}
	if len(s.Range(0, 1, 0.25)) != 5 {
		t.Error("s.Range(0, 1, 0.25) should have 5 elements, but is", s.Range(0, 1, 0.25))
	}
	if len(s.Range(0, 1, 0.3)) != 4 {
		t.Error("s.Range(0, 1, 0.3) should have 4 elements, but is", s.Range(0, 1, 0.3))
	}
}
