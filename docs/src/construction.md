```@meta
CurrentModule = MobiusTransformations
```

# Construction

A Möbius transformation is built from its four coefficients, from a `2×2`
matrix, or from the images of three points. `Mobius` is an ASCII alias for
`Möbius`.

```julia
m = Möbius(1, 2, 3, 4)          # z -> (z + 2) / (3z + 4)
n = Möbius([1 2; 3 4])          # same map, from its matrix
f = Möbius(0, 1, Inf, 1, 2, 3)  # 0 -> 1, 1 -> 2, Inf -> 3
```

The constructor overloads:

| Form | Result |
| --- | --- |
| `Möbius(a, b, c, d)` | build a `MöbiusTransformation` (generic) |
| `Möbius(M::AbstractMatrix)` | from the `2×2` matrix `[a b; c d]` |
| `Möbius(x, y, z)` | `[0, 1, Inf]` → `[x, y, z]` |
| `Möbius(x, y, z, X, Y, Z)` | `[x, y, z]` → `[X, Y, Z]` |
| `Möbius(target)` | vector → 3-arg form |
| `Möbius(source, target)` | vectors → 6-arg form |
| `Möbius(::Type{T}=ComplexF64)` | identity map |

```@docs
MöbiusTransformation
Möbius
Mobius
```
