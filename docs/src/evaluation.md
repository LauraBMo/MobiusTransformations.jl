```@meta
CurrentModule = MobiusTransformations
```

# Evaluation

A `MöbiusTransformation` is callable, mapping a point `z` to

```math
m(z) = \frac{a z + b}{c z + d}.
```

Evaluation at `Inf` is handled, and a pole (vanishing denominator) returns the
value set by [`set_infinity`](@ref). Maps broadcast elementwise over
collections.

```julia
m = Möbius(1, 2, 3, 4)
m(0), m(1), m(Inf)
m.([0, 1, Inf])
```

```@docs
MöbiusTransformation(z)
Base.broadcastable
```
