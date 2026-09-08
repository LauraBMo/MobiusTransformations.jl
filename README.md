# MobiusTransformations.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://LauraBMo.github.io/MobiusTransformations.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://LauraBMo.github.io/MobiusTransformations.jl/dev/) [![Build Status](https://github.com/LauraBMo/MobiusTransformations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/LauraBMo/MobiusTransformations.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/LauraBMo/MobiusTransformations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/LauraBMo/MobiusTransformations.jl)

A Möbius transformation is a rational function of the form

$$f(z) = \frac{a z + b}{c z + d}, \qquad ad - bc \neq 0,$$

mapping the extended complex plane (the Riemann sphere) to itself.
This package provides a lightweight `MöbiusTransformation` type with constructors, evaluation, composition, and projective equality.

## Installation

```julia
julia> ] add MobiusTransformations
```

## Usage

```julia
using MobiusTransformations
```

A transformation is built from its four coefficients:

```julia
julia> m = Möbius(1, 2, 3, 4)          # z -> (z + 2) / (3z + 4)
Möbius map z --> (1*z + 2) / (3*z + 4)

julia> m(0), m(1), m(Inf)              # evaluation, Inf included
(0.5, 0.42857142857142855, 0.3333333333333333)
```

Map three points to three points (using `Inf` to fix the point at infinity):

```julia
julia> f = Möbius(0, 1, Inf, 1+im, 2+2im, 3+3im)   # 0 -> 1+i, 1 -> 2+2i, Inf -> 3+3i

julia> f(0), f(1), f(Inf)
(1.0 + 1.0im, 2.0 + 2.0im, 3.0 + 3.0im)
```

Compose, invert, and compare up to scaling:

```julia
julia> g = Möbius(0, 1, 1, 0)          # z -> 1/z
julia> isone(g * inv(g))               # g and its inverse cancel
true

julia> Möbius(1, 2, 3, 4) == Möbius(2, 4, 6, 8)   # scalar multiples are equal
true
```

Möbius maps broadcast over collections, and expose their matrix:

```julia
julia> g.([1, 2, 4])                   # elementwise 1/z
3-element Vector{Float64}:
 1.0
 0.5
 0.25

julia> Matrix(m)                       # [a b; c d]
2×2 Matrix{Int64}:
 1  2
 3  4
```

## API

| Function | Description |
| --- | --- |
| `Möbius(a, b, c, d)` | Transformation `z -> (a z + b) / (c z + d)` |
| `Möbius(M)` | Transformation from the `2×2` matrix `[a b; c d]` |
| `Möbius(x, y, z)` | Maps `(0, 1, Inf)` to `(x, y, z)` |
| `Möbius(x, y, z, X, Y, Z)` | Maps `(x, y, z)` to `(X, Y, Z)` |
| `Möbius([...])`, `Möbius(src, tgt)` | Vector forms of the above |
| `Möbius()` / `Möbius(T)` | Identity map (defaults to `ComplexF64`) |
| `m(z)` | Evaluate the map at `z` |
| `m * n`, `m ∘ n` | Composition |
| `inv(m)` | Inverse map |
| `LinearAlgebra.det(m)`, `LinearAlgebra.normalize(m)` | Determinant and SL-normalization |
| `Matrix(m)` | Coefficients as a `2×2` matrix |
| `isone(m)`, `m == n`, `hash(m)` | Projective identity/equality/hash |
| `set_infinity(x)` | Value a pole evaluates to (default `complex(Inf)`) |

`Mobius` is an ASCII alias for `Möbius`.

## Exact fields and Nemo

The three-point construction is adapted from `mobius.jl` in
[ComplexRegions.jl](https://github.com/complexvariables/ComplexRegions.jl).
Unlike that code, which branches on `isinf` and returns a hardcoded
`complex(abs(num)/abs(den))` at a pole, this package stores coefficients
parametrically and uses a local `_isinf`, so maps work over any field — including
exact fields where `isinf` isn't defined. Over a Nemo number field, point the
package at the field's own infinity:

```julia
using MobiusTransformations, Nemo
R, x = polynomial_ring(QQ, "x")
K, a = number_field(x^2 + 1, "a")     # ℚ(i)
m = Möbius(K(1), K(2), K(3), K(4))    # z ↦ (z + 2)/(3z + 4)
set_infinity(Nemo.inf)
m(-K(4)//K(3))                        # pole → Nemo.inf, not complex(Inf)
```

If you currently use ComplexRegions.jl's `mobius.jl` for Möbius maps, this package
is a standalone drop-in with exact-field support.

## Attribution

This package is derived from the Möbius-transform code in
[ComplexRegions.jl](https://github.com/complexvariables/ComplexRegions.jl)
(© Tobin Driscoll, MIT). See `LICENSE` for details.

## License

MIT.
