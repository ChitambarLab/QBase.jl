using Test, LinearAlgebra

@testset "QMath/matrices.jl" begin

using QBase: QMath

@testset "block_diagonal()" begin
    I2 = [2 0; 0 2]
    P2 = [0 2; 2 0]
    I5 = Matrix(1I, 5, 5)
    four_two = [1 2;3 4;5 6;7 8]
    two_four = [1 2 3 4;5 6 7 8]

    @test QMath.block_diagonal([I5]) == I5
    @test QMath.block_diagonal([four_two]) == four_two
    @test QMath.block_diagonal([two_four]) == two_four

    @test QMath.block_diagonal([I2,P2,I2]) == [
        2 0 0 0 0 0;
        0 2 0 0 0 0;
        0 0 0 2 0 0;
        0 0 2 0 0 0;
        0 0 0 0 2 0;
        0 0 0 0 0 2
    ]

    @test QMath.block_diagonal([four_two,two_four]) == [
        1 2 0 0 0 0;
        3 4 0 0 0 0;
        5 6 0 0 0 0;
        7 8 0 0 0 0;
        0 0 1 2 3 4;
        0 0 5 6 7 8;
    ]

    @test QMath.block_diagonal([[1][:,:], I2]) == [
        1 0 0;
        0 2 0;
        0 0 2;
    ]
end

@testset "partial_trace()" begin
    @testset "invalid inputs" begin
        @test_throws DomainError QMath.partial_trace([1 0;0 1], [2,2], 1)
        @test_throws DomainError QMath.partial_trace([1 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0], [2,2], 1)
    end

    @testset "diagonal density operator" begin
        ρ = (1/6)*Matrix(1I,6,6)

        @test QMath.partial_trace(ρ, [2,3],2) == [0.5 0; 0 0.5]
        @test QMath.partial_trace(ρ, [2,3],1) == (1/3)*I
        @test QMath.partial_trace(ρ, [3,2],2) == (1/3)*I
    end

    @testset "maximally entangled bell state" begin
        ρ = 0.5*[1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]

        @test QMath.partial_trace(ρ,[2,2],2) == [0.5 0; 0 0.5]
        @test QMath.partial_trace(ρ,[2,2],1) == [0.5 0; 0 0.5]
    end

    @testset "DensityMatrix types work with partial trace" begin
        ρ = QBase.States.DensityMatrix([1 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0])

        reduced_ρ = QMath.partial_trace(ρ, [2,2], 1)

        @test !isa(reduced_ρ, QBase.States.DensityMatrix)
        @test reduced_ρ isa Array{Complex{Float64}, 2}
        @test reduced_ρ == [1 0;0 0]
    end
end

@testset "computational_basis_vectors()" begin
    @test QMath.computational_basis_vectors(2) == [[1,0],[0,1]]
    @test QMath.computational_basis_vectors(5) == [
        [1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]
    ]
end

@testset "commutes()" begin
    @test_throws ArgumentError QMath.commutes([1 0;0 1;1 1],[0 1;1 0])

    @test QMath.commutes([1 0;0 0],[0 0;0 1])
    @test !QMath.commutes([0 1;1 0],[1 0;0 -1])

    @test QMath.commutes.([[0 1;1 0],[1 0;0 -1]],[[0 im;-im 0],[0 im;-im 0]]) == [false,false]
end

end
