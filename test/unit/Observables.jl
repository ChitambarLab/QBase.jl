using Test, LinearAlgebra

@testset "QBase/Observables.jl" begin

using QBase

@testset "Observables.is_povm()" begin
    @testset "positive test case from Nielsen and Chuang" begin
        E1 = sqrt(2)/(1 + sqrt(2)) * [0 0;0 1]
        E2 = sqrt(2)/(1 + sqrt(2))*(1/2) * [1 -1; -1 1]
        E3 = [1 0;0 1] - E1 - E2

        @test Observables.is_povm([E1,E2,E3])

        @test Observables.is_povm([
            [1/3 0; 0 1/3],
            2/3*[-cos(2*π/3) 0.5*sin(2*π/3); 0.5*sin(2*π/3) -cos(2*π/3)],
            2/3*[-cos(2*π/3) -0.5*sin(2*π/3); -0.5*sin(2*π/3) -cos(2*π/3)]
        ])
    end

    @testset "negative test case" begin
        @test !Observables.is_povm([[1 0;0 0],[0 0;0 1],[0 0;0 0.1]])
    end
end


@testset "Observables.POVM()" begin
    valid_povms = [
        [[1 0;0 0],[0 0;0 1]],
        [[0.5 (0.1-0.1*im);(0.1+0.1*im) 0.5],[0.5 -(0.1-0.1*im);-(0.1+0.1*im) 0.5]],
        [[2/3 0;0 0],[1/6 1/sqrt(12);1/sqrt(12) 1/2],[1/6 -1/sqrt(12);-1/sqrt(12) 1/2]]
    ]

    @testset "valid povms: $id/$(length(valid_povms))" for id in 1:length(valid_povms)
        Π = Observables.POVM(valid_povms[id])
        @test isa(Π, Observables.POVM)
        @test Π == valid_povms[id]
    end

    invalid_povms = [
        [[1 0;0 0],[0 0;0 1],[0.5 0;0 0.5]],
        [[0.5 (0.1+0.1*im);(0.1+0.1*im) 0.5],[0.5 -(0.1+0.1*im);-(0.1+0.1*im) 0.5]],
        [[0.5 0.5;0.5 0.5],[0.5 -0.5*im;-0.5*im 0.5]]
    ]

    @testset "invalid povms: $id/$(length(invalid_povms))" for id in 1:length(invalid_povms)
        @test_throws DomainError Observables.POVM(invalid_povms[id])
    end
end

@testset "Observables.QubitPOVM()" begin

    valid_cases = [
        [[0.5 0.5;0.5 0.5],[0.5 -0.5;-0.5 0.5]]
    ]

    num_valid_cases = length(valid_cases)

    @testset "valid qubit povms: $i/$num_valid_cases" for i in 1:num_valid_cases
        Π = Observables.QubitPOVM(valid_cases[i])
        @test Π isa Observables.QubitPOVM
        @test Π == valid_cases[i]
    end

    @testset "invalid povms" begin
        @test_throws DomainError Observables.QubitPOVM([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 1]])
        @test_throws DomainError Observables.QubitPOVM([[1 0;0 0],[0 0;0 0.9]])
    end
end

@testset "mirror_symmetric_qubit_3povm" begin
    @testset "valid input θ = $θ" for θ in π/4:0.1:π/2
        Π = Observables.mirror_symmetric_qubit_3povm(θ)
        @test Observables.is_povm(Π)
        @test isa(Π, Observables.QubitPOVM)
    end

    @testset "invalid inputs" begin
        @test_throws DomainError Observables.mirror_symmetric_qubit_3povm(π/4 - 0.1)
        @test_throws DomainError Observables.mirror_symmetric_qubit_3povm(π/2 + 0.1)
    end
end

@testset "asymmetric_qubit_3povm" begin
    @testset "valid range of angles" begin
        @test isa(Observables.asymmetric_qubit_3povm(π/2, -0.0001), Observables.QubitPOVM)
        @test isa(Observables.asymmetric_qubit_3povm(0.00001, -π/2), Observables.QubitPOVM)

        for θ1 in π/24:π/24:π/2
            for θ2 in -π/2:π/24:(θ1-π/2 - 0.00001) # slight offset because θ2 cannot be 0
                Π1 = Observables.asymmetric_qubit_3povm(θ1, θ2)
                Π2 = Observables.asymmetric_qubit_3povm(θ2, θ1)

                @test isa(Π1, Observables.QubitPOVM)
                @test isa(Π2, Observables.QubitPOVM)
            end
        end

        for θ1 in -π/24:-π/24:-π/2
            for θ2 in π/2:-π/24:(π/2+θ1+0.00001)
                Π1 = Observables.asymmetric_qubit_3povm(θ1, θ2)
                Π2 = Observables.asymmetric_qubit_3povm(θ2, θ1)

                @test isa(Π1, Observables.QubitPOVM)
                @test isa(Π2, Observables.QubitPOVM)
            end
        end
    end

    @testset "invalid arguments" begin
        @test_throws DomainError Observables.asymmetric_qubit_3povm(π/2, 0)
        @test_throws DomainError Observables.asymmetric_qubit_3povm(π/2-0.0001, 0)
        @test_throws DomainError Observables.asymmetric_qubit_3povm(0, -π/2)
        @test_throws DomainError Observables.asymmetric_qubit_3povm(0, -π/2+0.001)
        @test_throws DomainError Observables.asymmetric_qubit_3povm(π/3, π/4)
    end
end

@testset "sqrt_povm()" begin
    priors = QMath.Marginals([1/3;1/3;1/3])

    states = States.trine_qubits
    Π = Observables.sqrt_povm(priors, states)

    @test isa(Π, Observables.POVM)
end

@testset "sic_qubit_povm" begin
    @test isa(Observables.sic_qubit_povm, Observables.QubitPOVM)
end

@testset "trine_qubit_povm" begin
    @test isa(Observables.trine_qubit_povm, Observables.QubitPOVM)
    @test Observables.trine_qubit_povm ≈ Observables.mirror_symmetric_qubit_3povm(π/3)
end

@testset "kraus_operators()" begin

    @testset "trine qubit povm" begin
        k = Observables.kraus_operators(Observables.trine_qubit_povm)

        @test k[1]'*k[1] ≈ Observables.trine_qubit_povm[1]
        @test k[2]'*k[2] ≈ Observables.trine_qubit_povm[2]
        @test k[3]'*k[3] ≈ Observables.trine_qubit_povm[3]
    end

    @testset "pentagram qubit povm" begin
        Π = Observables.planar_symmetric_qubit_povm(5)
        k = Observables.kraus_operators(Π)

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
        projector_dict = Observables.naimark_dilation(Observables.trine_qubit_povm)

        ρ_set = States.trine_qubits

        @test 2/3 ≈ tr(projector_dict["projectors"][1] * kron(ρ_set[1], projector_dict["ancilla"]))
        @test 2/3 ≈ tr(projector_dict["projectors"][2] * kron(ρ_set[2], projector_dict["ancilla"]))
        @test 2/3 ≈ tr(projector_dict["projectors"][3] * kron(ρ_set[3], projector_dict["ancilla"]))

        @test 1/6 ≈ tr(projector_dict["projectors"][1] * kron(ρ_set[2], projector_dict["ancilla"]))
        @test 1/6 ≈ tr(projector_dict["projectors"][2] * kron(ρ_set[1], projector_dict["ancilla"]))
        @test 1/6 ≈ tr(projector_dict["projectors"][2] * kron(ρ_set[3], projector_dict["ancilla"]))
    end
end

@testset "planar_symmetric_qubit_povm()" begin
    @testset "scanning over simple cases" begin
        for n in 2:100
            Π = Observables.planar_symmetric_qubit_povm(n)

            @test length(Π) == n
            @test Π isa Observables.QubitPOVM
        end
    end
end

end
