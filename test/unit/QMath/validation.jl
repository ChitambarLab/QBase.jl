using Test

@testset "QMath/validation.jl" begin

using QBase: QMath

@testset "is_hermitian()" begin
    @test QMath.is_hermitian([1 0;0 1])
    @test QMath.is_hermitian([1 0;0 -1])
    @test !( QMath.is_hermitian([1 0;1 1]) )
    @test !( QMath.is_hermitian([1 1;-1 1]) )
    @test QMath.is_hermitian([1 im; -im 1])
end

@testset "is_positive_semidefinite()" begin
    @test !( QMath.is_positive_semidefinite([1 0 ; 0 -1]) )
    @test !( QMath.is_positive_semidefinite([0 1;1 0]) )
    @test !( QMath.is_positive_semidefinite([1 1;-1 1]) )
    @test QMath.is_positive_semidefinite([1 0; 0 1])

    @test_throws ArgumentError QMath.is_positive_semidefinite([1 0 0;0 1 0])
end

@testset "is_square_matrix()" begin
    @test QMath.is_square_matrix(fill(1,(3,3)))
    @test !( QMath.is_square_matrix(fill(1,(3,2))) )
    @test !( QMath.is_square_matrix(fill(1,(2,3))) )
    @test !( QMath.is_square_matrix(fill(1,(3,3,3))) )
end

end
