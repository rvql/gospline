# gospline
Golang cubic spline library

This library generates a cubic spline for given points.

## Usage

### Create a cubic spline
```golang
s := gospline.NewCubicSpline([]float64{0, 1, 2, 3}, []float64{0, 0.5, 2, 1.5})
```

### Get an interpolated value
```golang
s.At(3.5)
```

### Get an array of interpolated values
```
s.Range(0, 3, 0.25)
```

### Supported boundaries
First derivation boundary: gospline.NewClampedCubicSpline
Second derivation boundary: gospline.NewNaturalCubicSpline

## Installation
Just
`go get github.com/cnkei/gospline`
