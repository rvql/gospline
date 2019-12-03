package gospline

import "math"

const epsilon = 1e-9

func floatEquals(a, b float64) bool {
	return math.Abs(a-b) < epsilon
}
