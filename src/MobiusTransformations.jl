module MobiusTransformations

export Mobius, Möbius, set_infinity

const INF = Ref{Any}(complex(Inf))

"""
    set_infinity(infinity)

Sets the representation of infinity used by the package.
"""
function set_infinity(infinity)
    INF[] = infinity
end

# Möbius transformation
"""
    MöbiusTransformation{T}

Represents a Möbius transformation of the form:
f(z) = (a*z + b) / (c*z + d)

# Fields
- `a::T, b::T, c::T, d::T`: Coefficients of the transformation
"""
struct MöbiusTransformation{T} <: Function
    a::T
    b::T
    c::T
    d::T
end

# Dealing with diferent types
MöbiusTransformation(a, b, c, d) = MöbiusTransformation(promote(a, b, c, d)...)

# Möbius transformations via Möbius(...) calls:
# 4 args: Build MöbiusTransformation (generic)
# 3 args: Standard transformation: [0, 1, Inf] -> target Number-like
# 6 args: Transformation: source Number-like -> target Number-like
# 2 args: Source, Target vectors -> 6 args
# 1 args: Target vector -> 3 args
# 1 args::Type{T} (arg is a Type) returns Identity (default ComplexF64)
# 1 args::Matrix{T} MöbiusTransformation(args')

"""
    Möbius(a, b, c, d)

Creates a [`MöbiusTransformation`](@ref) with coefficients `a`, `b`, `c`, `d`.
"""
Möbius(a, b, c, d) = MöbiusTransformation(a, b, c, d)

const Mobius = Möbius   # for us lazy Americans

"""
    Möbius(M::AbstractMatrix)

Creates a Möbius transformation with coefficients `a`, `b`, `c`, `d`, where `M=[a b; c d]`.
"""
Möbius(M::AbstractMatrix) = Möbius((transpose(M))...)

# The three-point construction below (mapping 0, 1, ∞ to given points via the
# isinf branches) is adapted from ComplexRegions.jl by Tobin Driscoll (MIT):
# https://github.com/complexvariables/ComplexRegions.jl/blob/main/src/mobius.jl
"""
    Möbius(x, y, z)

Returns the Möbius transformation that maps `[0, 1, Inf]` to points `[x, y, z]`.
Values of `Inf` are permitted.
"""
function Möbius(x, y, z)
    if isinf(x)
        return Möbius(z, y - z, one(x), zero(x))
    elseif isinf(y)
        return Möbius(-z, x, -one(x), one(x))
    elseif isinf(z)
        return Möbius(y - x, x, zero(x), one(x))
    else
        xy, yz = y - x, z - y
        return Möbius(z * xy, x * yz, xy, yz)
    end
end

"""
    Möbius(x, y, z, X, Y, Z)

Returns the Möbius transformation that maps `[x, y, z]` to `[X, Y, Z]`.
Values of `Inf` are permitted.
"""
function Möbius(x, y, z, X, Y, Z)
    m_source = Möbius(x, y, z) # (0, 1, Inf) -> (x, y, z)
    m_image = Möbius(X, Y, Z) # (0, 1, Inf) -> (X, Y, Z)
    # Obs. No 'inv' performed.
    return m_image * inv(m_source) # (x, y, z) -> (X, Y, Z)
end

"""
    Möbius(target)

Returns the Möbius transformation that maps `[0, 1, Inf]` to `target = [x, y z]`.
Values of `Inf` are permitted.
"""
Möbius(target) = Möbius(target...)

"""
    Möbius(source, target)

Returns the Möbius transformation that maps `source` to `target`.
Values of `Inf` are permitted.
"""
Möbius(source, target) = Möbius(source..., target...)

"""
    Möbius(::Type{T}=ComplexF64)

Returns the identity Möbius transformation of type `T`.
"""
Möbius(::Type{T}=ComplexF64) where T = Möbius(one(T), zero(T), zero(T), one(T))

end # of module MobiusTransformations.
