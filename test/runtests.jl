using MobiusTransformations
using Test
using Aqua
using LinearAlgebra: det, normalize

const MT = MobiusTransformations

@testset "MobiusTransformations.jl" begin

    @testset "constructors" begin
        # Möbius / Mobius alias
        @test Mobius === Möbius

        # 4-arg constructor promotes mixed types
        @test eltype(Möbius(1, 2, 3, 4)) == Int
        @test Möbius(1, 2, 3, 4) isa MT.MöbiusTransformation{Int}
        @test eltype(Möbius(1, 2, 3.0, 4)) == Float64
        @test eltype(Möbius(1.0 + 0.0im, 2, 3, 4)) == ComplexF64

        # fields are stored as given
        m = Möbius(1, 2, 3, 4)
        @test m.a == 1 && m.b == 2 && m.c == 3 && m.d == 4

        # identity maps
        @test Möbius() isa MT.MöbiusTransformation{ComplexF64}
        @test isone(Möbius())
        @test Möbius(ComplexF64) isa MT.MöbiusTransformation{ComplexF64}
        @test Möbius(Int) isa MT.MöbiusTransformation{Int}
        @test Möbius(Float64) isa MT.MöbiusTransformation{Float64}
        @test isone(Möbius(Int))
        @test isone(Möbius(0, 1, Inf))

        # matrix constructor (row-wise, not conjugated)
        @test Möbius([1 2; 3 4]) == Möbius(1, 2, 3, 4)
        @test Möbius([1 + 2im 3; 4 5]) == Möbius(1 + 2im, 3, 4, 5)
    end

    @testset "three-point Möbius(x, y, z)" begin
        # all finite: maps 0 -> x, 1 -> y, Inf -> z
        x, y, z = 1 + 2im, 3 - 1im, 2 + 5im
        m = Möbius(x, y, z)
        @test m(0) ≈ x
        @test m(1) ≈ y
        @test m(Inf) ≈ z

        # Inf in each slot
        m = Möbius(Inf, 1 + 1im, 2 + 2im)
        @test isinf(m(0))
        @test m(1) ≈ 1 + 1im
        @test m(Inf) ≈ 2 + 2im

        m = Möbius(1 + 1im, Inf, 2 + 2im)
        @test m(0) ≈ 1 + 1im
        @test isinf(m(1))
        @test m(Inf) ≈ 2 + 2im

        m = Möbius(1 + 1im, 2 + 2im, Inf)
        @test m(0) ≈ 1 + 1im
        @test m(1) ≈ 2 + 2im
        @test isinf(m(Inf))
    end

    @testset "six-point Möbius(x, y, z, X, Y, Z)" begin
        m = Möbius(0, 1, Inf, 1 + 1im, 2 + 2im, 3 + 3im)
        @test m(0) ≈ 1 + 1im
        @test m(1) ≈ 2 + 2im
        @test m(Inf) ≈ 3 + 3im

        m = Möbius(1, 2, 3, 2 + 1im, 0, 4 - 1im)
        @test m(1) ≈ 2 + 1im
        @test m(2) ≈ 0
        @test m(3) ≈ 4 - 1im
    end

    @testset "vector forms" begin
        m = Möbius([0, 1, Inf])
        @test m(0) == 0
        @test m(1) == 1
        @test isinf(m(Inf))

        m = Möbius([1, 2, 3], [10, 20, 30])
        @test m(1) ≈ 10
        @test m(2) ≈ 20
        @test m(3) ≈ 30
    end

    @testset "evaluation" begin
        m = Möbius(1, 2, 3, 4)
        @test m(0) ≈ 0.5
        @test m(1) ≈ 3 / 7
        @test m(-2) ≈ 0
        @test m(1im) ≈ (1im + 2) / (3im + 4)
        @test m(Inf) ≈ 1 / 3
        @test m(-Inf) ≈ 1 / 3             # -Inf is the same point as Inf
        @test m(complex(0, Inf)) ≈ 1 / 3  # imaginary infinity too

        # z -> 1/z has a pole at 0
        invmap = Möbius(0, 1, 1, 0)
        @test isinf(invmap(0))
        @test invmap(Inf) ≈ 0
        @test invmap(2) ≈ 0.5
    end

    @testset "composition and inverse" begin
        m = Möbius(1, 2, 3, 4)
        n = Möbius(5, 6, 7, 8)

        # m * n applies n then m; agrees with ∘ and matrix product
        mn = m * n
        @test mn == (m ∘ n)
        @test Matrix(mn) ≈ Matrix(m) * Matrix(n)
        for z in (-1.5, 0.0, 0.7, 2 + 3im)
            @test mn(z) ≈ m(n(z))
        end

        # inverse round-trips to the identity
        g = Möbius(2 + 1im, 3, 4, 5 - 2im)
        @test isone(g * inv(g))
        @test isone(inv(g) * g)
        for z in (0, 1, 2 - 3im)
            @test inv(g)(g(z)) ≈ z
        end
        @test isinf(inv(g)(g(Inf)))
    end

    @testset "det and normalize" begin
        m = Möbius(1, 2, 3, 4)
        @test det(m) == 1 * 4 - 2 * 3

        # SL-normalization: det -> 1, projectively equal
        g = Möbius(1 + 1im, 2, 3, 4 - 1im)
        @test det(normalize(g)) ≈ 1
        @test normalize(g) == g

        # normalize requires a real square root of a negative det
        @test_throws DomainError normalize(Möbius(1, 2, 3, 4))
    end

    @testset "scalar multiplication" begin
        m = Möbius(1, 2, 3, 4)

        # regression: coefficients must scale in (a, b, c, d) order, not
        # the column-major order a splatted matrix would produce
        @test Matrix(2 * m) == [2 4; 6 8]
        @test 2 * m == Möbius(2, 4, 6, 8)
        @test m * 2 == 2 * m
        @test (1 + 2im) * m == Möbius(1 + 2im, 2 + 4im, 3 + 6im, 4 + 8im)

        # det scales quadratically
        @test det(2 * m) == 4 * det(m)
    end

    @testset "isone and equality" begin
        @test isone(Möbius(1, 0, 0, 1))
        @test isone(Möbius(2, 0, 0, 2))          # scalar multiple of identity
        @test !isone(Möbius(1, 1, 0, 1))

        # projective equality: scalar multiples are equal
        @test Möbius(1, 2, 3, 4) == Möbius(2, 4, 6, 8)
        @test Möbius(1, 2, 3, 4) == Möbius(1 + 0im, 2, 3, 4)
        @test Möbius(1, 2, 3, 4) != Möbius(1, 2, 3, 5)

        # hash is consistent with equality
        @test hash(Möbius(1, 2, 3, 4)) == hash(Möbius(2, 4, 6, 8))
        # -0.0 normalizes to 0.0 (negative scalar multiple of the identity)
        @test hash(Möbius(1, 0, 0, 1)) == hash(Möbius(-1, 0, 0, -1))
    end

    @testset "Matrix and broadcasting" begin
        @test Matrix(Möbius(1, 2, 3, 4)) == [1 2; 3 4]

        m = Möbius(0, 1, 1, 0)  # z -> 1/z
        @test m.([1, 2, 4]) ≈ [1, 0.5, 0.25]
    end

    @testset "display" begin
        m = Möbius(1, 2, 3, 4)
        @test occursin("Möbius", sprint(show, m))
        @test occursin("Möbius", sprint(show, MIME("text/plain"), m))
    end

    @testset "Aqua quality checks" begin
        Aqua.test_all(MobiusTransformations)
    end

end
