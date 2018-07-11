package gospline

import (
	"fmt"
	"testing"
)

func TestSimpleCubicSpline(t *testing.T) {
	x := [...]float64{0, 1, 2, 3}
	y := [...]float64{0, 0.5, 2, 1.5}
	s := newSpline(x[:], y[:], CubicSecondDeriv, -0.3, 3.3)
	fmt.Println(s.At(0))
	fmt.Println(s.At(1))
	fmt.Println(s.At(2))
	fmt.Println(s.At(3))
}
