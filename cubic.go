package gospline

type boundary uint

// bounday types
const (
	CubicNatual boundary = iota
	CubicClamped
	CubicPeriodic
)

type cubic struct {
	x []float64
	y []float64
	boundary
	// donates to f''0 and f''n if boundary == CubicNatual
	// or f'0 and f'n if boundary == CubicClamped
	f0 float64
	fn float64
}
