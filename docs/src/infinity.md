```@meta
CurrentModule = MobiusTransformations
```

# Set Infinity

A Möbius transformation is a map on the Riemann sphere, so the package needs a
concrete value for the point at infinity — in particular, to return when a map
evaluates at a pole (where the denominator vanishes). That value is stored in
`INF` and set with `set_infinity`.

It plays two roles:

- **returned** — a pole evaluates to `INF[]`;
- **accepted** — `m(INF[])` or `Mobius([0, 1, INF[]])` are recognized.

The default is `complex(Inf)`, so plain `Inf` works out of the box.

!!! note
    `set_infinity` is strictly needed only to choose the value a pole
    **returns**.  When `isinf` is defined for a type `T`, the infinity values of that type are **accepted** on input, regardless of `INF[]`.

## Over a Nemo number field

```julia
using Nemo
R, x = polynomial_ring(QQ, "x")
K, a = number_field(x^2 + 1, "a")   # QQ(i)

m = Möbius(K(1), K(2), K(3), K(4))  # z -> (z + 2) / (3z + 4)

m(-K(4)//K(3))   # pole: default INF[] = complex(Inf), a ComplexF64 (not in K)
m(Nemo.inf)      # accepted on input via isinf → a/c = 1/3

set_infinity(Nemo.inf)              # poles now return Nemo's infinity
m(-K(4)//K(3))   # denominator 3·(-4/3) + 4 = 0 → Nemo.inf
```

`Nemo.inf` is AbstractAlgebra's positive infinity. Storing it in `INF` makes
poles return it; on input, `m(Nemo.inf)` is recognized as infinity regardless,
since `isinf(Nemo.inf)` is `true`.

```@docs
INF
set_infinity
```
