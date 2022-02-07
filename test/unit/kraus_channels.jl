using Test
using Suppressor
using QBase

@testset "is_kraus_channel()" begin
    @testset "valid cases" begin
        @test is_kraus_channel([[1 0;0 0],[0 0;0 1]])
        @test is_kraus_channel([[1 0;0 0],[0 0;0 0.9]], atol=0.2)
        @test is_kraus_channel([
            sqrt(1/2)*σI, sqrt(1/6)*σx, sqrt(1/6)*σy, sqrt(1/6)*σz,
        ])
    end

    @testset "invalid cases" begin
        warn_msg = @capture_err begin
            @test !is_kraus_channel([[1 0;0 0],[0 0;0 0.9]])
        end

        @test occursin("`kraus_ops` do not sum to identity.", warn_msg)
    end
end

@testset "KrausChannel()" begin
    @testset "KrausChannel(::Vector{Matrix})" begin
        kraus_ops = [[1 0;0 0],[0 0;0 1]]
        𝒩 = KrausChannel(kraus_ops)

        @test 𝒩 isa KrausChannel{Int}
        @test 𝒩.kraus_ops == kraus_ops
    end

    @testset "DomainErrors" begin
        @suppress_err @test_throws DomainError KrausChannel([[1 0;0 0],[1 0;0 0]])
    end

    @testset "show(::KrausChannel)" begin
        show_msg = @capture_out begin
            show(stdout, MIME("text/plain"), KrausChannel([
                [1 0;0 0],[0 0;0 1]
            ]))
        end

        @test show_msg == """KrausChannel{Int64}
        atol : 1.0e-7
        kraus_ops : 2-element Array{Array{Int64,2},1}:
         [1 0; 0 0]
         [0 0; 0 1]"""
    end
end

@testset "kraus_evolve()" begin
    kraus_ops = [
        sqrt(1/2)*σI,
        sqrt(1/6)*σx,
        sqrt(1/6)*σy,
        sqrt(1/6)*σz,
    ]

    ρ_out = kraus_evolve(kraus_ops, [1 0;0 0])
    @test ρ_out ≈ [2/3 0;0 1/3]
end
