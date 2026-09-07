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

end # of module MobiusTransformations.
