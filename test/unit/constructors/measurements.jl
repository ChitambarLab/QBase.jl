using Test, LinearAlgebra
using QBase

@testset "./src/constructors/measurements.jl" begin

@testset "mirror_symmetric_qubit_3povm" begin
    @testset "valid input θ = $θ" for θ in π/4:0.1:π/2
        Π = mirror_symmetric_qubit_3povm(θ)
        @test is_povm(Π)
        @test Π isa POVM{Float64}
    end

    @testset "invalid inputs" begin
        @test_throws DomainError mirror_symmetric_qubit_3povm(π/4 - 0.1)
        @test_throws DomainError mirror_symmetric_qubit_3povm(π/2 + 0.1)
    end
end

@testset "asymmetric_qubit_3povm" begin
    @testset "valid range of angles" begin
        @test asymmetric_qubit_3povm(π/2, -0.0001) isa POVM{Float64}
        @test asymmetric_qubit_3povm(0.00001, -π/2) isa POVM{Float64}

        for θ1 in π/24:π/24:π/2
            for θ2 in -π/2:π/24:(θ1-π/2 - 0.00001) # slight offset because θ2 cannot be 0
                Π1 = asymmetric_qubit_3povm(θ1, θ2)
                Π2 = asymmetric_qubit_3povm(θ2, θ1)

                @test Π1 isa POVM{Float64}
                @test Π2 isa POVM{Float64}
            end
        end

        for θ1 in -π/24:-π/24:-π/2
            for θ2 in π/2:-π/24:(π/2+θ1+0.00001)
                Π1 = asymmetric_qubit_3povm(θ1, θ2)
                Π2 = asymmetric_qubit_3povm(θ2, θ1)

                @test Π1 isa POVM{Float64}
                @test Π2 isa POVM{Float64}
            end
        end
    end

    @testset "invalid arguments" begin
        @test_throws DomainError asymmetric_qubit_3povm(π/2, 0)
        @test_throws DomainError asymmetric_qubit_3povm(π/2-0.0001, 0)
        @test_throws DomainError asymmetric_qubit_3povm(0, -π/2)
        @test_throws DomainError asymmetric_qubit_3povm(0, -π/2+0.001)
        @test_throws DomainError asymmetric_qubit_3povm(π/3, π/4)
    end
end

@testset "planar_symmetric_qubit_povm()" begin
    @testset "scanning over simple cases" begin
        for n in 2:50
            Π= planar_symmetric_qubit_povm(n)

            @test length(Π) == n
            @test Π isa POVM{Float64}
        end
    end
end

@testset "trine_qubit_povm()" begin
    trine_Π = trine_qubit_povm()
    @test trine_Π isa POVM{Float64}
    @test isapprox(trine_Π, mirror_symmetric_qubit_3povm(π/3), atol=1e-7)
end

@testset "sic_qubit_povm()" begin
    @test sic_qubit_povm() isa POVM{Complex{Float64}}
end

@testset "sqrt_povm()" begin
    priors = Probabilities([1/3;1/3;1/3])

    states = trine_qubit_states()
    Π = sqrt_povm(priors, states)

    @test Π isa POVM{Float64}
    @test isapprox(Π, trine_qubit_povm(), atol=1e-7)
end

@testset "_naimark_kraus_operators()" begin
    @testset "trine qubit povm" begin
        k = QBase._naimark_kraus_operators(trine_qubit_povm())
        trine_Π = trine_qubit_povm()

        @test k[1]'*k[1] ≈ trine_Π[1]
        @test k[2]'*k[2] ≈ trine_Π[2]
        @test k[3]'*k[3] ≈ trine_Π[3]
    end

    @testset "pentagram qubit povm" begin
        Π = planar_symmetric_qubit_povm(5)
        k = QBase._naimark_kraus_operators(Π)

        @test sum(k_i -> k_i'*k_i, k) ≈ [1 0;0 1]

        @test k[1]'*k[1] ≈ Π[1]
        @test k[2]'*k[2] ≈ Π[2]
        @test k[3]'*k[3] ≈ Π[3]
        @test k[4]'*k[4] ≈ Π[4]
        @test k[5]'*k[5] ≈ Π[5]
    end
end

@testset "naimark_dilation()" begin
    @testset "reconstructs trine probabilites" begin
        projector_dict = naimark_dilation(trine_qubit_povm())

        ρ_set = trine_qubit_states()

        @test 2/3 ≈ tr(Matrix{Complex{Float64}}(
            projector_dict["projectors"][1] * kron(ρ_set[1], projector_dict["ancilla"])
        ))
        @test 2/3 ≈ tr(Matrix{Complex{Float64}}(
            projector_dict["projectors"][2] * kron(ρ_set[2], projector_dict["ancilla"])
        ))
        @test 2/3 ≈ tr(Matrix{Complex{Float64}}(
            projector_dict["projectors"][3] * kron(ρ_set[3], projector_dict["ancilla"])
        ))

        @test 1/6 ≈ tr(Matrix{Complex{Float64}}(
            projector_dict["projectors"][1] * kron(ρ_set[2], projector_dict["ancilla"])
        ))
        @test 1/6 ≈ tr(Matrix{Complex{Float64}}(
            projector_dict["projectors"][2] * kron(ρ_set[1], projector_dict["ancilla"])
        ))
        @test 1/6 ≈ tr(Matrix{Complex{Float64}}(
            projector_dict["projectors"][2] * kron(ρ_set[3], projector_dict["ancilla"])
        ))
    end
end

end
