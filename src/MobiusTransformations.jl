module MobiusTransformations

using UnPack: @unpack
import LinearAlgebra: det, normalize

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

"""
    Mobius

ASCII alias for `Möbius`. Provided for convenience when typing `ö`
is cumbersome.
"""
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

"""
    isone(m::MöbiusTransformation)

Return `true` if `m` is the identity Mobius transformation and `false` otherwise.
"""
function Base.isone(m::MöbiusTransformation)
    @unpack a, b, c, d = m
    return iszero(b) && iszero(c) && (a == d)
end

"""
    ==(m::MöbiusTransformation, n::MöbiusTransformation)

Return `true` if `m` and `n` are equal as projective transformations, i.e.
`m * inv(n)` is the identity.
"""
Base.:(==)(m::MöbiusTransformation, n::MöbiusTransformation) = isone(m * inv(n))

# TODO
# ≈(m::MöbiusTransformation, n::MöbiusTransformation) = isapproxone(m*inv(n))

"""
    eltype(m::MöbiusTransformation{T}) where {T}

Return the coefficient type `T` of `m`.
"""
Base.eltype(_::MöbiusTransformation{T}) where {T} = T

"""
    hash(m::MöbiusTransformation, h::UInt64)

Return a hash for `m`, based on its values at `0`, `1`, and `Inf`.
"""
function Base.hash(m::MöbiusTransformation, h::UInt64=UInt64(0))
    z = 0.0 + 0.0 * im # kludge to make -0.0 and -0.0im into +versions
    a = m(0) + z
    b = m(1) + z
    c = m(Inf) + z
    return hash(a, hash(b, hash(c, h)))
end

# Vectorized operations
"""
    broadcastable(m::MöbiusTransformation)

Return `Ref(m)` so that `m.(z)` broadcasts elementwise over `z` rather than
treating `m` as a collection.
"""
Base.broadcastable(m::MöbiusTransformation) = Ref(m)

"""
    Matrix(m::MöbiusTransformation)

Return the `2×2` coefficient matrix `[a b; c d]` of `m`.
"""
function Base.Matrix(m::MöbiusTransformation)
    @unpack a, b, c, d = m
    return [a b; c d]
end

"""
    det(m::MöbiusTransformation)

Return the determinant `a*d - b*c` of the coefficient matrix of `m=Möbius(a, b, c, d)`.
"""
function det(m::MöbiusTransformation)
    @unpack a, b, c, d = m
    return a * d - b * c
end

"""
    normalize(m::MöbiusTransformation)

Returns a Möbius transformation `m2` such that `m2 == m` and `det(m2) = 1`.
Requires `sqrt(det(m))` to be defined in the coefficient field.
"""
normalize(m::MöbiusTransformation) = m * inv(sqrt(det(m)))

# Inverse Möbius transformation
"""
    inv(m::MöbiusTransformation)

Return the inverse Möbius transformation. For `m = Möbius(a, b, c, d)`,
`inv(m)` has coefficients `(d, -b, -c, a)`.
"""
function Base.inv(m::MöbiusTransformation)
    @unpack a, b, c, d = m
    MöbiusTransformation(d, -b, -c, a)
end

"""
    *(λ, m::MöbiusTransformation)
    *(m::MöbiusTransformation, λ)

Scale every coefficient of `m` by the scalar `λ`.
"""
Base.:(*)(λ, m::MöbiusTransformation) = Möbius(λ.*Matrix(m))
Base.:(*)(m::MöbiusTransformation, λ) = *(λ, m)

"""
    *(m::MöbiusTransformation, n::MöbiusTransformation)

Compose two Möbius transformations.
"""
function Base.:(*)(m::MöbiusTransformation, n::MöbiusTransformation)
    @unpack a, b, c, d = n
    e, f, g, h = a, b, c, d
    @unpack a, b, c, d = m
    MöbiusTransformation(a * e + b * g, a * f + b * h,
                         c * e + d * g, c * f + d * h)
end

"""
    ∘(m::MöbiusTransformation, n::MöbiusTransformation)

Compose two Möbius transformations (same as *).
"""
Base.:∘(m::MöbiusTransformation, n::MöbiusTransformation) = m * n

#
# Eval
# 
"""
    (m::MöbiusTransformation)(z)

Apply to a number z the Möbius transformation
`m(z) = (a*z + b) / (c*z + d)`, where `m = Möbius([a b; c d])`.

Values of `Inf` are permitted.
"""
function (m::MöbiusTransformation)(z)
    @unpack a, b, c, d = m
    if isinf(z)
        numer, denom = a, c
    else
        numer, denom = a * z + b, c * z + d
    end

    if abs(denom) == 0
        return INF[]
    else
        return numer * inv(denom)
    end
end

#
# Display
#
function Base.show(io::IO, m::MöbiusTransformation)
    @unpack a, b, c, d = m
    print(IOContext(io, :compact => true), "Möbius map z --> ($a*z + $b) / ($c*z + $d)")
end

function Base.show(io::IO, ::MIME"text/plain", m::MöbiusTransformation)
    string_linear((X, Y)) = "(" * X * ")*z + " * Y
    @unpack a, b, c, d = m
    # if abs(c) < NUM_TOL
    #     A, B = [repr("text/plain", x) for x in [a*inv(d), b*inv(d)]]
    #     numer = string_linear((A, B))
    #     print(io, "Möbius:\n   ", numer)
    # end
    A, B, C, D = [repr("text/plain", x) for x in [a, b, c, d]]
    numer, denom = string_linear.([(A, B), (C, D)])
    newline = "\n   "
    hline = reduce(*, fill("–", maximum(length, [numer, denom])))

    ## Print message
    print(io, "Möbius: ", eltype(m), newline,
        numer, newline,
        hline, newline,
        denom)
end

end # of module MobiusTransformations.
