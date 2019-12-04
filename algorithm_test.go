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

func TestFindSegment(t *testing.T) {
	testcases := []struct {
		xs     []float64
		x      float64
		result int
	}{
		// in segment
		{
			xs:     []float64{1, 2, 3, 4, 5},
			x:      4.5,
			result: 3,
		},
		// on endpoint
		{
			xs:     []float64{1, 2, 3, 4, 5},
			x:      3,
			result: 2,
		},
		// negative endpoints
		{
			xs:     []float64{-1.2, -1.0, 0.6, 1.3, 100},
			x:      0,
			result: 1,
		},
		// below any
		{
			xs:     []float64{-1.2, -1.0, 0.6, 1.3, 100},
			x:      -100,
			result: 0,
		},
		// above any
		{
			xs:     []float64{-1.2, -1.0, 0.6, 1.3, 100},
			x:      101,
			result: 3,
		},
		// left endpoint
		{
			xs:     []float64{-1.2, -1.0, 0.6, 1.3, 100},
			x:      -1.2,
			result: 0,
		},
		// right endpoint
		{
			xs:     []float64{-1.2, -1.0, 0.6, 1.3, 100},
			x:      100,
			result: 3,
		},
	}
	for _, tc := range testcases {
		if result := findSegment(tc.xs, tc.x); result != tc.result {
			t.Errorf("testcase %v failed, expected %d, result %d", tc, tc.result, result)
		}
	}
}
