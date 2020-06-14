using Test, LinearAlgebra

@testset "QBase/Observables.jl" begin

using QBase

@testset "QBase.is_povm()" begin
    @testset "positive test case from Nielsen and Chuang" begin
        E1 = sqrt(2)/(1 + sqrt(2)) * [0 0;0 1]
        E2 = sqrt(2)/(1 + sqrt(2))*(1/2) * [1 -1; -1 1]
        E3 = [1 0;0 1] - E1 - E2

        @test QBase.is_povm([E1,E2,E3])

        @test QBase.is_povm([
            [1/3 0; 0 1/3],
            2/3*[-cos(2*π/3) 0.5*sin(2*π/3); 0.5*sin(2*π/3) -cos(2*π/3)],
            2/3*[-cos(2*π/3) -0.5*sin(2*π/3); -0.5*sin(2*π/3) -cos(2*π/3)]
        ])
    end

    @testset "negative test case" begin
        @test !QBase.is_povm([[1 0;0 0],[0 0;0 1],[0 0;0 0.1]])
    end
end


@testset "QBase.POVM()" begin
    valid_povms = [
        [[1 0;0 0],[0 0;0 1]],
        [[0.5 (0.1-0.1*im);(0.1+0.1*im) 0.5],[0.5 -(0.1-0.1*im);-(0.1+0.1*im) 0.5]],
        [[2/3 0;0 0],[1/6 1/sqrt(12);1/sqrt(12) 1/2],[1/6 -1/sqrt(12);-1/sqrt(12) 1/2]]
    ]

    @testset "valid povms: $id/$(length(valid_povms))" for id in 1:length(valid_povms)
        Π = QBase.POVM(valid_povms[id])
        @test isa(Π, QBase.POVM)
        @test Π == valid_povms[id]
    end

    invalid_povms = [
        [[1 0;0 0],[0 0;0 1],[0.5 0;0 0.5]],
        [[0.5 (0.1+0.1*im);(0.1+0.1*im) 0.5],[0.5 -(0.1+0.1*im);-(0.1+0.1*im) 0.5]],
        [[0.5 0.5;0.5 0.5],[0.5 -0.5*im;-0.5*im 0.5]]
    ]

    @testset "invalid povms: $id/$(length(invalid_povms))" for id in 1:length(invalid_povms)
        @test_throws DomainError QBase.POVM(invalid_povms[id])
    end
end

@testset "QBase.QubitPOVM()" begin

    valid_cases = [
        [[0.5 0.5;0.5 0.5],[0.5 -0.5;-0.5 0.5]]
    ]

    num_valid_cases = length(valid_cases)

    @testset "valid qubit povms: $i/$num_valid_cases" for i in 1:num_valid_cases
        Π = QBase.QubitPOVM(valid_cases[i])
        @test Π isa QBase.QubitPOVM
        @test Π == valid_cases[i]
    end

    @testset "invalid povms" begin
        @test_throws DomainError QBase.QubitPOVM([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 1]])
        @test_throws DomainError QBase.QubitPOVM([[1 0;0 0],[0 0;0 0.9]])
    end
end

@testset "mirror_symmetric_qubit_3povm" begin
    @testset "valid input θ = $θ" for θ in π/4:0.1:π/2
        Π = QBase.mirror_symmetric_qubit_3povm(θ)
        @test QBase.is_povm(Π)
        @test isa(Π, QBase.QubitPOVM)
    end

    @testset "invalid inputs" begin
        @test_throws DomainError QBase.mirror_symmetric_qubit_3povm(π/4 - 0.1)
        @test_throws DomainError QBase.mirror_symmetric_qubit_3povm(π/2 + 0.1)
    end
end

@testset "asymmetric_qubit_3povm" begin
    @testset "valid range of angles" begin
        @test isa(QBase.asymmetric_qubit_3povm(π/2, -0.0001), QBase.QubitPOVM)
        @test isa(QBase.asymmetric_qubit_3povm(0.00001, -π/2), QBase.QubitPOVM)

        for θ1 in π/24:π/24:π/2
            for θ2 in -π/2:π/24:(θ1-π/2 - 0.00001) # slight offset because θ2 cannot be 0
                Π1 = QBase.asymmetric_qubit_3povm(θ1, θ2)
                Π2 = QBase.asymmetric_qubit_3povm(θ2, θ1)

                @test isa(Π1, QBase.QubitPOVM)
                @test isa(Π2, QBase.QubitPOVM)
            end
        end

        for θ1 in -π/24:-π/24:-π/2
            for θ2 in π/2:-π/24:(π/2+θ1+0.00001)
                Π1 = QBase.asymmetric_qubit_3povm(θ1, θ2)
                Π2 = QBase.asymmetric_qubit_3povm(θ2, θ1)

                @test isa(Π1, QBase.QubitPOVM)
                @test isa(Π2, QBase.QubitPOVM)
            end
        end
    end

    @testset "invalid arguments" begin
        @test_throws DomainError QBase.asymmetric_qubit_3povm(π/2, 0)
        @test_throws DomainError QBase.asymmetric_qubit_3povm(π/2-0.0001, 0)
        @test_throws DomainError QBase.asymmetric_qubit_3povm(0, -π/2)
        @test_throws DomainError QBase.asymmetric_qubit_3povm(0, -π/2+0.001)
        @test_throws DomainError QBase.asymmetric_qubit_3povm(π/3, π/4)
    end
end

@testset "sqrt_povm()" begin
    priors = QBase.QMath.Marginals([1/3;1/3;1/3])

    states = QBase.trine_qubits
    Π = QBase.sqrt_povm(priors, states)

    @test isa(Π, QBase.POVM)
end

@testset "sic_qubit_povm" begin
    @test isa(QBase.sic_qubit_povm, QBase.QubitPOVM)
end

@testset "trine_qubit_povm" begin
    @test isa(QBase.trine_qubit_povm, QBase.QubitPOVM)
    @test QBase.trine_qubit_povm ≈ QBase.mirror_symmetric_qubit_3povm(π/3)
end

@testset "measurement_probabilites()" begin
    @testset "trine_measurements" begin
        conditionals = QBase.measurement_probabilities(
            QBase.mirror_symmetric_qubit_3povm(π/3), QBase.trine_qubits
        )
        @test conditionals ≈ [2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3]
        @test conditionals isa QBase.QMath.Conditionals
    end

    @testset "classical_measurements" begin
        Π = QBase.QubitPOVM([[1 0;0 0],[0 0;0 1],[0 0;0 0]])
        ρ_set = QBase.Qubit.([[1 0;0 0],[1 0;0 0],[0 0;0 1]])

        @test QBase.measurement_probabilities(Π,ρ_set) == [1 1 0;0 0 1;0 0 0]

        Π = QBase.POVM([[1 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]])
        ρ_set = QBase.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]])

        @test QBase.measurement_probabilities(Π,ρ_set) == [1 1 0;0 0 1]
    end
end

@testset "kraus_operators()" begin
    k = QBase.kraus_operators(QBase.trine_qubit_povm)

    @test k[1]'*k[1] ≈ QBase.trine_qubit_povm[1]
    @test k[2]'*k[2] ≈ QBase.trine_qubit_povm[2]
    @test k[3]'*k[3] ≈ QBase.trine_qubit_povm[3]
end

@testset "naimark_dilation()" begin
    @testset "reconstructs trine probabilites" begin
        projector_dict = QBase.naimark_dilation(QBase.trine_qubit_povm)

        ρ_set = QBase.trine_qubits

        @test 2/3 ≈ tr(projector_dict["projectors"][1] * kron(ρ_set[1], projector_dict["ancilla"]))
        @test 2/3 ≈ tr(projector_dict["projectors"][2] * kron(ρ_set[2], projector_dict["ancilla"]))
        @test 2/3 ≈ tr(projector_dict["projectors"][3] * kron(ρ_set[3], projector_dict["ancilla"]))

        @test 1/6 ≈ tr(projector_dict["projectors"][1] * kron(ρ_set[2], projector_dict["ancilla"]))
        @test 1/6 ≈ tr(projector_dict["projectors"][2] * kron(ρ_set[1], projector_dict["ancilla"]))
        @test 1/6 ≈ tr(projector_dict["projectors"][2] * kron(ρ_set[3], projector_dict["ancilla"]))
    end
end

end
