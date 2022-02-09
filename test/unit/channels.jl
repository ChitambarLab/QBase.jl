using Test
using LinearAlgebra
using QBase

@testset "./src/channels.jl" begin

@testset "replacer_channel()" begin
    @testset "simple qubit examples" begin
        ρ = State([1 0;0 0])
        σ = State([0 0;0 1])

        r = replacer_channel(ρ, σ, 0.5)

        @test r isa State
        @test r == [0.5 0;0 0.5]
    end
    @testset "errors" begin
        ρ = State([1 0;0 0])
        σ = State([0 0;0 1])

        @test_throws DomainError replacer_channel(ρ, σ, -0.1)
        @test_throws DomainError replacer_channel(ρ, σ, 1.1)
        @test_throws DomainError replacer_channel(
            State([1 0;0 0]), State([0 0 0;0 1 0;0 0 0]), 0.5
        )
    end
end

@testset "depolarizing_channel()" begin
    @testset "simple State examples" begin
        ρ0 = State([1 0;0 0])

        ρ0_out1 = depolarizing_channel(ρ0, 1)

        @test ρ0_out1 isa State
        @test ρ0_out1 == ρ0

        ρ0_out2 = depolarizing_channel(ρ0, 0)

        @test ρ0_out2 isa State
        @test ρ0_out2 == [1/2 0;0 1/2]

        ρ0_out3 = depolarizing_channel(ρ0, 0.5)

        @test ρ0_out3 isa State
        @test ρ0_out3 == [3/4 0;0 1/4]
    end

    @testset "simple State examples" begin
        ρ = State([0.5 0 0.5;0 0 0;0.5 0 0.5])

        ρ_out1 = depolarizing_channel(ρ, 1)

        @test ρ_out1 isa State
        @test ρ_out1 == ρ

        ρ_out2 = depolarizing_channel(ρ, 0.5)

        @test ρ_out2 isa State
        @test ρ_out2 == [(1/4+1/6) 0 1/4;0 1/6 0;1/4 0 (1/4+1/6)]
    end

    @testset "errors" begin
        ρ = State([1 0;0 0])

        @test_throws DomainError depolarizing_channel(ρ, 1.1)
        @test_throws DomainError depolarizing_channel(ρ, -0.1)
    end
end

@testset "erasure_channel()" begin
    @testset "simple State examples" begin
        ρ = State([1 0;0 0])

        ρ_out1 = erasure_channel(ρ, 1)

        @test ρ_out1 isa State
        @test ρ_out1 == [1 0 0;0 0 0;0 0 0]

        ρ_out2 = erasure_channel(ρ, 0)

        @test ρ_out2 isa State
        @test ρ_out2 == [0 0 0;0 0 0;0 0 1]

        ρ_out3 = erasure_channel(ρ, 0.5)

        @test ρ_out3 isa State
        @test ρ_out3 == [0.5 0 0;0 0 0;0 0 0.5]
    end

    @testset "simple State examples" begin
        ρ = State([0.5 0 0.5;0 0 0;0.5 0 0.5])

        ρ_out1 = erasure_channel(ρ, 1)

        @test ρ_out1 isa State
        @test ρ_out1 == [0.5 0 0.5 0;0 0 0 0;0.5 0 0.5 0;0 0 0 0]

        ρ_out2 = erasure_channel(ρ, 0.5)

        @test ρ_out2 isa State
        @test ρ_out2 == [1/4 0 1/4 0;0 0 0 0;1/4 0 1/4 0;0 0 0 1/2]
    end

    @testset "errors" begin
        ρ = State([1 0;0 0])

        @test_throws DomainError erasure_channel(ρ, 1.1)
        @test_throws DomainError erasure_channel(ρ, -0.1)
    end
end

end
