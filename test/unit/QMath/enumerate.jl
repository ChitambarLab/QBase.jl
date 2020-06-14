using Test, LinearAlgebra

@testset "QMath/enumerate.jl" begin

using QBase: QMath

@testset "permutation_maps()" begin
    @testset "dim = 2" begin
        @test QMath.permutation_maps(2) == [[1 0; 0 1], [0 1; 1 0]]
    end

    @testset "dim = 3" begin
        @test issetequal(QMath.permutation_maps(3), [
            [1 0 0; 0 1 0; 0 0 1],
            [1 0 0; 0 0 1; 0 1 0],
            [0 1 0; 1 0 0; 0 0 1],
            [0 1 0; 0 0 1; 1 0 0],
            [0 0 1; 1 0 0; 0 1 0],
            [0 0 1; 0 1 0; 1 0 0],
        ])
    end

    @testset "permtutation map validity" for N in 2:6
        perms = QMath.permutation_maps(N)

        @testset "maps are unique" begin
            @test size(perms) == (factorial(N),)
            @test size(perms) == size(unique(perms))
        end

        @testset "matrix dimension" begin
            for perm in perms
                @test size(perm) == (N,N)
            end
        end

        @testset "determinant is (+/-) 1" begin
            for perm in perms
                parity = det(perm)
                @test (parity == -1) | (parity == 1)
            end
        end

        @testset "invertibility" begin
            for perm in perms
                @test perm*perm' == diagm(0 => fill(1,N))
            end
        end
    end
end

@testset "base_n_val()" begin
    @testset "little endian" begin
        @test 1 == QMath.base_n_val(digits(1,base=2,pad=2), 2, endian="little")
        @test 57 == QMath.base_n_val(digits(57,base=3,pad=5), 3, endian="little")
    end

    @testset "big endian" begin
        @test 12 == QMath.base_n_val([1,1,0,0],2)
        @test 974 == QMath.base_n_val([1,2,3,4,4],5)
    end

    @testset "argument errors" begin
        @test_throws ArgumentError QMath.base_n_val([1,2,0,5],5)
        @test_throws ArgumentError QMath.base_n_val([1,2,0,5],4)
        @test_throws ArgumentError QMath.base_n_val([-1,0,1],2)
    end
end

@testset "stirling2()" begin
    @test QMath.stirling2(1,1) == 1

    @testset "n = 2" begin
        @test QMath.stirling2(2,1) == 1
        @test QMath.stirling2(2,2) == 1
    end

    @testset "n = 3" begin
        @test QMath.stirling2(3,1) == 1
        @test QMath.stirling2(3,2) == 3
        @test QMath.stirling2(3,3) == 1
    end

    @testset "n = 4" begin
        @test QMath.stirling2(4,1) == 1
        @test QMath.stirling2(4,2) == 7
        @test QMath.stirling2(4,3) == 6
        @test QMath.stirling2(4,4) == 1
    end

    @testset "n = 5" begin
        @test QMath.stirling2(5,1) == 1
        @test QMath.stirling2(5,2) == 15
        @test QMath.stirling2(5,3) == 25
        @test QMath.stirling2(5,4) == 10
        @test QMath.stirling2(5,5) == 1
    end

    @testset "n = 6" begin
        @test QMath.stirling2(6,1) == 1
        @test QMath.stirling2(6,2) == 31
        @test QMath.stirling2(6,3) == 90
        @test QMath.stirling2(6,4) == 65
        @test QMath.stirling2(6,5) == 15
        @test QMath.stirling2(6,6) == 1
    end

    @testset "n = 10" begin
        @test QMath.stirling2(10,1) == 1
        @test QMath.stirling2(10,2) == 511
        @test QMath.stirling2(10,3) == 9330
        @test QMath.stirling2(10,4) == 34105
        @test QMath.stirling2(10,5) == 42525
        @test QMath.stirling2(10,6) == 22827
        @test QMath.stirling2(10,7) == 5880
        @test QMath.stirling2(10,8) == 750
        @test QMath.stirling2(10,9) == 45
        @test QMath.stirling2(10,10) == 1
    end
end
end
