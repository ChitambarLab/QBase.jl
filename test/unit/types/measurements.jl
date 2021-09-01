using Test, LinearAlgebra
using QBase

@testset "./src/types/measurements.jl" begin

@testset "is_povm_element()" begin
    @testset "valid cases" begin
        @test is_povm_element([1 0;0 0])
        @test is_povm_element([0.5 -0.5im;0.5im 0.5])
        @test is_povm_element([0.5 0.2;0.3 0.5], atol=0.11)
    end

    @testset "invalid cases" begin
        @test !is_povm_element([0.5 0.2;0.3 0.5])
        @test !is_povm_element([1 0.1;0.1 0])
    end
end

@testset "is_povm()" begin
    @testset "valid cases" begin
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

@testset "POVMel()" begin
    @testset "valid povm elements" begin
        M = POVMel(State([1 0 0;0 0 0;0 0 0]))
        @test M isa POVMel{Int64}
        @test M == [1 0 0;0 0 0;0 0 0]

        M = POVMel([0.5 0.5im;-0.5im 0.5])
        @test M isa POVMel{Complex{Float64}}
        @test M == [0.5 0.5im;-0.5im 0.5]

        M = POVMel([0.5 0.45im;-0.5im 0.5], atol=0.06)
        @test M isa POVMel{Complex{Float64}}
        @test M == [0.5 0.45im;-0.5im 0.5]
    end

    @testset "errors" begin
        @test_throws DomainError POVMel([0.5 0.5im;0.5im 0.5])
    end
end

@testset "kron for POVMel" begin
    M1 = POVMel([1 0;0 0])
    M2 = POVMel([0 0;0 1.])

    M1_M2 =  kron(M1,M2)
    @test M1_M2 isa POVMel{Float64}
    @test M1_M2 == [0 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0]
end

@testset "POVM()" begin
    @testset "valid povms:" begin
        Π = [[1 0;0 0],[0 0;0 1]]
        povm = POVM(Π)
        @test povm isa POVM{Int64}
        @test povm isa Measurement
        @test povm == Π
        @test length(povm) == 2

        Π = [[0.5 (0.1-0.1*im);(0.1+0.1*im) 0.5],[0.5 -(0.1-0.1*im);-(0.1+0.1*im) 0.5]]
        povm = POVM(Π)
        @test povm isa POVM{Complex{Float64}}
        @test povm isa Measurement
        @test povm == Π
        @test length(povm) == 2

        Π = [[2/3 0;0 0],[1/6 1/sqrt(12);1/sqrt(12) 1/2],[1/6 -1/sqrt(12);-1/sqrt(12) 1/2]]
        povm = POVM(Π)
        @test povm isa POVM{Float64}
        @test povm isa Measurement
        @test povm == Π
        @test length(povm) == 3
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

@testset "PVM()" begin
    @testset "valid cases" begin
        Π = [[1,0],[0,1]]
        pvm = PVM(Π)
        @test pvm isa PVM{Int64}
        @test pvm isa Measurement
        @test pvm == Π
        @test length(pvm) == 2

        pvm = PVM(Ket.(Π))
        @test pvm isa PVM{Int64}
        @test pvm isa Measurement
        @test pvm == Π

        Π = [[1,-1,-im,-im],[1,-1,im,im],[1,1,im,-im],[1,1,-im,im]]/2
        pvm = PVM(Π)
        @test pvm isa PVM{Complex{Float64}}
        @test pvm isa Measurement
        @test pvm == Π
        @test length(pvm) == 4
    end

    @testset "invalid cases" begin
        @test_throws MethodError PVM([[1 0;0 0],[0 0;0 1]])
        @test_throws DomainError PVM([[1,0,0],[0,1,0]])
        @test_throws DomainError PVM([[0.5,0.5],[0.5,-0.5]])
    end

    @testset "atol" begin
        @test_throws DomainError PVM([[1,1],[1,1]])
        @test PVM([[1,1],[1,1]], atol=2) isa PVM{Int64}
    end
end

end
