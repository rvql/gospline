package gospline

import (
	"math"
	"testing"
)

func TestSolveTridiagonal(t *testing.T) {
	const epsilon = 1e-9
	a := [4]float64{2, 3, 4, 5}
	b := [5]float64{1, 2, 3, 4, 5}
	c := [4]float64{2, 3, 4, 5}
	d := [5]float64{5, 15, 31, 53, 45}
	x := triThomas(a[:], b[:], c[:], d[:])
	if math.Abs(x[0]-1) > epsilon {
		t.Error("x[0] should be 1")
	}
	if math.Abs(x[1]-2) > epsilon {
		t.Error("x[1] should be 2")
	}
	if math.Abs(x[4]-5) > epsilon {
		t.Error("x[4] should be 5")
	}
}
