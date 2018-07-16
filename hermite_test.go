package gospline

import (
	"fmt"
	"testing"
)

func TestHermiteSpline(t *testing.T) {
	const epsilon = 1e-9
	x := [...]float64{0, 0.16, 0.42, 0.6425, 0.8575, 1.0}
	y := [...]float64{0, 32, 237, 255, 0, 0}
	s := NewMonotoneSpline(x[:], y[:])
	fmt.Println(s.At(0))
	fmt.Println(s.At(0.16))
	fmt.Println(s.At(0.42))
	fmt.Println(s.At(0.6425))
	fmt.Println(s.At(0.8575))
}
