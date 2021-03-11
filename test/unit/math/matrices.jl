using Test, LinearAlgebra
using QBase

@testset "./src/math/matrices.jl" begin

@testset "partial_trace()" begin
    @testset "invalid inputs" begin
        @test_throws DomainError partial_trace([1 0;0 1], [2,2], 1)
        @test_throws DomainError partial_trace([1 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0], [2,2], 1)
    end

    @testset "diagonal density operator" begin
        ρ = (1/6)*Matrix(1I,6,6)

        @test partial_trace(ρ, [2,3],2) == [0.5 0; 0 0.5]
        @test partial_trace(ρ, [2,3],1) == (1/3)*I
        @test partial_trace(ρ, [3,2],2) == (1/3)*I
    end

    @testset "maximally entangled bell state" begin
        ρ = 0.5*[1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]

        @test partial_trace(ρ,[2,2],2) == [0.5 0; 0 0.5]
        @test partial_trace(ρ,[2,2],1) == [0.5 0; 0 0.5]
    end
end

@testset "computational_basis_vectors()" begin
    @test computational_basis_vectors(2) == [[1,0],[0,1]]
    @test computational_basis_vectors(5) == [
        [1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]
    ]
end

@testset "commutes()" begin
    @test_throws DomainError commutes([1 0;0 1;1 1],[0 1;1 0])

    @test commutes([1 0;0 0],[0 0;0 1])
    @test !commutes([0 1;1 0],[1 0;0 -1])

    @test commutes.([[0 1;1 0],[1 0;0 -1]],[[0 im;-im 0],[0 im;-im 0]]) == [false,false]
end

@testset "is_hermitian()" begin
    @testset "invalid cases" begin
        @test !is_hermitian([1 0;1 1])
        @test !is_hermitian([1 1;-1 1])
    end
    @testset "valid cases" begin
        @test is_hermitian([1 0;0 1])
        @test is_hermitian([1 0;0 -1])
        @test is_hermitian([1 im; -im 1])
    end
    @testset "using atol" begin
        @test !is_hermitian([1 1;1+1e-6 1])
        @test is_hermitian([1 1;1+1e-6 1], atol=1e-5)
    end
end

@testset "is_positive_semidefinite()" begin
    @testset "invalid cases" begin
        @test !is_positive_semidefinite([1 0 ; 0 -1])
        @test !is_positive_semidefinite([0 1;1 0])
        @test !is_positive_semidefinite([1 1;-1 1])
        @test !is_positive_semidefinite([1 0 0;0 1 0])
    end

    @testset "valid cases" begin
        @test is_positive_semidefinite([1 0; 0 1])
    end

    @testset "using atol" begin
        @test !is_positive_semidefinite([1 0;0 -1e-6])
        @test is_positive_semidefinite([1 0;0 -1e-6], atol=1e-6)
    end
end

@testset "n_product_id()" begin
    m0 = [1 0;0 0]
    m1 = [0 0;1 0]
    m2 = [0 0 0;0 0 0;0 1 0;0 0 0]
    m3 = [0 0 0 0;0 0 0 0;0 0 0 1]

    col_id = n_product_id([1,1],[2,2])
    row_id = n_product_id([1,2],[2,2])
    @test col_id == 1
    @test row_id == 2
    @test kron(m0,m1)[row_id,col_id] == 1

    col_id = n_product_id([1,1,2,4],[2,2,3,4])
    row_id = n_product_id([1,2,3,3],[2,2,4,3])
    @test col_id == 8
    @test row_id == 21
    @test kron(m0,m1,m2,m3)[row_id,col_id] == 1

    @testset "errors" begin
        @test_throws DomainError n_product_id([1,1,1],[2,2])
        @test_throws DomainError n_product_id([1,3],[2,2])
    end
end

end
