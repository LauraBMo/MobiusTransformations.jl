```@meta
CurrentModule = MobiusTransformations
```

# MobiusTransformations.jl

A Möbius transformation is a rational function of the form

```math
f(z) = \frac{a z + b}{c z + d}, \qquad ad - bc \neq 0,
```

mapping the extended complex plane (the Riemann sphere) to itself. The package
provides a lightweight `MöbiusTransformation` type with constructors,
evaluation, composition, and projective equality.

## Installation

```julia
julia> ] add MobiusTransformations
```

## Usage

```julia
using MobiusTransformations

m = Möbius(1, 2, 3, 4)             # z -> (z + 2) / (3z + 4)
m(0), m(1), m(Inf)                 # evaluation, Inf included

f = Möbius(0, 1, Inf, 1+im, 2+2im, 3+3im)   # 0 -> 1+i, 1 -> 2+2i, Inf -> 3+3i
```

## Sections

- [Construction](construction.md) — build a `MöbiusTransformation`
- [Evaluation](evaluation.md) — apply a transformation to a point
- [Operations](operations.md) — invert, compose, compare, matrix, determinant, normalization
- [Set Infinity](infinity.md) — value at a pole

## License

MIT. Derived from the Möbius-transform code in
[ComplexRegions.jl](https://github.com/complexvariables/ComplexRegions.jl)
(© Tobin Driscoll, MIT).
