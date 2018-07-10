package gospline

// Spline interface for interpolation functions
type Spline interface {
	// At returns interpolated value at x
	// 0 <= x <= 1
	At(x float64) float64

	// Range returns interpolated values in [start, end] with step
	Range(start, end, step float64) []float64
}
