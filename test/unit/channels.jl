using Test
using LinearAlgebra
using Suppressor
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

@testset "is_choi_matrix()" begin

    @testset "valid cases" begin
        @test is_choi_matrix([1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1], [2,2])
        @test is_choi_matrix([1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 0.9], [2,2], atol=0.2)
    end

    @testset "invalid cases" begin
        warn_msg = @capture_err begin
            @test !is_choi_matrix([1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 0.9], [2,2])
        end

        @test occursin("The Choi matrix `Λ` is not trace-preserving.", warn_msg)
    end
end

@testset "choi_matrix()" begin
    𝒩_depol(X) = 1/2*[1 0 ; 0 1]
    Λ_depol = choi_matrix(𝒩_depol, [2,2])
    @test Λ_depol == [1 0 1 0;0 1 0 1;1 0 1 0;0 1 0 1] / 2
    @test is_choi_matrix(Λ_depol, [2,2])

    𝒩_id(X) = X
    Λ_id = choi_matrix(𝒩_id, [2,2])
    @test Λ_id == [1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1]
    @test is_choi_matrix(Λ_id, [2,2])
end

@testset "ChoiOp()" begin
        @testset "function instantiation" begin
            Λ_depol = ChoiOp(x -> 1/2*[1 0; 0 1], [2,2])

            @test Λ_depol isa Operator
            @test Λ_depol isa ChoiOp{ComplexF64}
            @test Λ_depol.M isa Matrix{ComplexF64}
            @test Λ_depol.M == [1 0 1 0;0 1 0 1;1 0 1 0;0 1 0 1]/2
            @test Λ_depol.dims == [2,2]
        end

        @testset "matrix instantiation" begin
            Λ = [1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1]
            Λ_id = ChoiOp(Λ, [2,2])

            @test Λ_id isa Operator
            @test Λ_id isa ChoiOp{Int}
            @test Λ_id.M isa Matrix{Int}
            @test Λ_id.M == Λ
            @test Λ_id.dims == [2,2]
        end

        #
        # @testset "kraus operator instantiation" begin
        #     kraus_ops = [σI, σx, σy, σz]/2
        #     choi_channel = Choi(kraus_ops)
        #
        #     @test choi_channel isa Choi{ComplexF64}
        #     @test choi_channel.JN isa Matrix{ComplexF64}
        #     @test choi_channel.JN == I/2
        #     @test choi_channel.in_dim == 2
        #     @test choi_channel.out_dim == 2
        # end
        #

        @testset "DomainError" begin
            @suppress_err @test_throws DomainError ChoiOp(Matrix(I, 4, 4), [2,2])
        end
end

@testset "choi_evolve()" begin
    Λ_depol = choi_matrix(x -> [1 0;0 1]/2, [2,2])
    ρ_in = [1 0;0 0]

    ρ_out = choi_evolve(Λ_depol, ρ_in)
    @test ρ_out == [1 0;0 1]/2
end

end
