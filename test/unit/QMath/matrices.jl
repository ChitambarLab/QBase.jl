using Test, LinearAlgebra

@testset "QMath/matrices.jl" begin

using QBase: QMath

@testset "block_diagonal()" begin
    I2 = [2 0; 0 2]
    P2 = [0 2; 2 0]
    I5 = diagm(0 => fill(5,5))
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

@testset "id()" begin
    @test QMath.id(2) == [1 0;0 1]
    @test QMath.id(1) == fill(1,(1,1))
    @test QMath.id(0) == fill(0,(0,0))
    @test QMath.id(3) == [1 0 0; 0 1 0; 0 0 1]

    @testset "dim = 31" begin
        id_31 = QMath.id(31)

        @test size(id_31) == (31,31)
        @test tr(id_31) == 31
    end
end

@testset "partial_trace()" begin

    @testset "diagonal density operator" begin
        ρ = (1/6)*QMath.id(6)

        @test QMath.partial_trace(ρ, (2,3),2) == [0.5 0; 0 0.5]
        @test QMath.partial_trace(ρ, (2,3),1) == (1/3)*QMath.id(3)
        @test QMath.partial_trace(ρ, (3,2),2) == (1/3)*QMath.id(3)
    end

    @testset "maximally entangled bell state" begin
        ρ = 0.5*[1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]

        @test QMath.partial_trace(ρ,(2,2),2) == [0.5 0; 0 0.5]
        @test QMath.partial_trace(ρ,(2,2),1) == [0.5 0; 0 0.5]

    end
end

@testset "commutes()" begin
    @test_throws ArgumentError QMath.commutes([1 0;0 1;1 1],[0 1;1 0])

    @test QMath.commutes([1 0;0 0],[0 0;0 1])
    @test !QMath.commutes([0 1;1 0],[1 0;0 -1])

    @test QMath.commutes.([[0 1;1 0],[1 0;0 -1]],[[0 im;-im 0],[0 im;-im 0]]) == [false,false]
end

end
