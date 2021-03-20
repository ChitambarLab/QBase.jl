using Test, LinearAlgebra
using QBase

@testset "./src/types/unitaries.jl" begin

@testset "is_unitary()" begin
    @testset "valid unitaries" begin
        valid_unitaries = [
            [0 1;1 0],
            [0 -im 0;im 0 0;0 0 1],
            [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0],
        ]

        for U in valid_unitaries
            @test is_unitary(U)
        end
    end

    @testset "invalid unitaries" begin
        invalid_unitaries = [
            [0.5 0;0 0.5],
            [0 im 0;-im 0 0;0 0 0],
            [1 0 0;0 1 0],
        ]

        for U in invalid_unitaries
            @test !is_unitary(U)
        end
    end
end

@testset "Unitary()" begin
    @testset "valid cases" begin
        U = Unitary([0 1;1 0])
        @test U == [0 1;1 0]
        @test U isa Unitary{Int64}
        @test U isa AbstractUnitary

        U = Unitary([1/sqrt(2) 1/sqrt(2) 0;-1im/sqrt(2) 1im/sqrt(2) 0;0 0 im])
        @test U == [1/sqrt(2) 1/sqrt(2) 0;-1im/sqrt(2) 1im/sqrt(2) 0;0 0 im]
        @test U isa Unitary{Complex{Float64}}
        @test U isa AbstractUnitary
    end

    @testset "invalid unitaries" begin
        @test_throws DomainError Unitary([0 1;1 0;0 0])
        @test_throws DomainError Unitary([0.5 0.5;0.5 0.5])
    end

    @testset "atol" begin
        @test_throws DomainError Unitary([1 0.1;0 1])
        @test Unitary([1 0.1;0 1], atol=0.2) isa Unitary{Float64}
    end
end

@testset "unitary multiplication" begin
    Id = σx*σx
    @test Id == [1 0;0 1]
    @test Id isa Unitary{Int64}

    U_prod = *(random_unitary.([3,3,3,3,3])...)
    @test U_prod isa Unitary{Complex{Float64}}
end

@testset "unitary kron" begin
    σ_xz = kron(σx, σz)
    @test σ_xz == [0 0 1 0;0 0 0 -1;1 0 0 0;0 -1 0 0]
    @test σ_xz isa Unitary{Int64}

    σ_yII = kron(σy,σI,σI)
    @test σ_yII == [zeros(4,4) -im*Matrix(I,4,4);im*Matrix(I,4,4) zeros(4,4)]
    @test σ_yII isa Unitary{Complex{Int64}}
end

end
