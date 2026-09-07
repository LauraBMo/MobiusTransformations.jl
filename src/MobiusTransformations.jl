module MobiusTransformations

export set_infinity

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

end # of module MobiusTransformations.
