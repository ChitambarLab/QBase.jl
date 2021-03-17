using Test, LinearAlgebra
using QBase

@testset "./src/types/measurements.jl" begin

@testset "is_povm()" begin
    @testset "positive test case from Nielsen and Chuang" begin
        E1 = sqrt(2)/(1 + sqrt(2)) * [0 0;0 1]
        E2 = sqrt(2)/(1 + sqrt(2))*(1/2) * [1 -1; -1 1]
        E3 = [1 0;0 1] - E1 - E2

        @test is_povm([E1,E2,E3])

        @test is_povm([
            [1/3 0; 0 1/3],
            2/3*[-cos(2*π/3) 0.5*sin(2*π/3); 0.5*sin(2*π/3) -cos(2*π/3)],
            2/3*[-cos(2*π/3) -0.5*sin(2*π/3); -0.5*sin(2*π/3) -cos(2*π/3)]
        ])
    end

    @testset "negative test case" begin
        @test !is_povm([[1 0;0 0],[0 0;0 1],[0 0;0 0.1]])
    end

    @testset "atol" begin
        Π = [[0.5 0.6;0.6 0.5],[0.5 -0.5;-0.5 0.5]]
        @test !is_povm(Π)
        @test is_povm(Π,atol=0.2)
    end

    @test is_povm(POVM([[1 0;0 0],[0 0;0 1]]))
end

@testset "POVM()" begin
    @testset "valid povms:" begin
        Π = [[1 0;0 0],[0 0;0 1]]
        povm = POVM(Π)
        @test povm isa POVM{Int64}
        @test povm isa AbstractMeasurement
        @test povm.Π == Π

        Π = [[0.5 (0.1-0.1*im);(0.1+0.1*im) 0.5],[0.5 -(0.1-0.1*im);-(0.1+0.1*im) 0.5]]
        povm = POVM(Π)
        @test povm isa POVM{Complex{Float64}}
        @test povm isa AbstractMeasurement
        @test povm.Π == Π

        Π = [[2/3 0;0 0],[1/6 1/sqrt(12);1/sqrt(12) 1/2],[1/6 -1/sqrt(12);-1/sqrt(12) 1/2]]
        povm = POVM(Π)
        @test povm isa POVM{Float64}
        @test povm isa AbstractMeasurement
        @test povm.Π == Π
    end

    @testset "invalid povms" begin
        invalid_povms = [
            [[1 0;0 0],[0 0;0 1],[0.5 0;0 0.5]],
            [[0.5 (0.1+0.1*im);(0.1+0.1*im) 0.5],[0.5 -(0.1+0.1*im);-(0.1+0.1*im) 0.5]],
            [[0.5 0.5;0.5 0.5],[0.5 -0.5*im;-0.5*im 0.5]]
        ]

        for Π in invalid_povms
            @test_throws DomainError POVM(Π)
        end
    end

    @testset "atol" begin
        Π = [[1 0;0 0],[0 0.1;0 1]]
        @test_throws DomainError POVM(Π)
        @test POVM(Π, atol=0.11) isa POVM{Float64}
    end
end

end
