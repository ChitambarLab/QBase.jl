using Test, LinearAlgebra

@testset "QBase/States.jl" begin

using QBase

@testset "is_ket()" begin
    @testset "valid wavefunctions" begin
        @test QBase.is_ket([1;0])
        @test QBase.is_ket([1])
        @test QBase.is_ket(1/sqrt(3)*[1;1;-1])
        @test QBase.is_ket(1/sqrt(2)*[1;im])
    end

    @testset "invalid wavefunctions" begin
        @test !( QBase.is_ket([0;0]) )
        @test !( QBase.is_ket([1;1]) )
    end

    @testset "valid inputs" begin
        @test QBase.is_ket([1,0]::Vector{Int64})
        @test QBase.is_ket([1.0,0.0]::Vector{Float64})
        @test QBase.is_ket([1.0+im*0.0,0.0+im*0.0 ]::Vector{Complex{Float64}})
    end

    @testset "invalid inputs" begin
        @test_throws MethodError QBase.is_ket([1 0;0 0])
        @test_throws MethodError QBase.is_ket(1/sqrt(2)*[1 1])
    end
end

@testset "QBase.Ket()" begin
    valid_cases = [
        [1;0;0],
        [1,0,0,0],
        [1,1,im]./sqrt(3),
        [0.5im,0,-sqrt(3)/2]
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ψ = QBase.Ket(valid_cases[id])
        @test ψ isa QBase.Ket
        @test ψ isa QBase.States.AbstractKet
        @test ψ == valid_cases[id]
    end

    @testset "invalid inputs" begin
        @test_throws DomainError QBase.Ket([1,1,1])
        @test_throws MethodError QBase.Ket([1 0])
    end
end

@testset "QBase.QubitKet()" begin
    valid_cases = [
        [1,0],
        [1,im]/sqrt(2),
        [0.5;sqrt(3)/2]
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ψ = QBase.QubitKet(valid_cases[id])
        @test ψ isa QBase.QubitKet
        @test ψ isa QBase.States.AbstractKet
        @test ψ == valid_cases[id]
    end

    @testset "invalid inputs" begin
        @test_throws DomainError QBase.QubitKet([1 0])
        @test_throws DomainError QBase.QubitKet([1,1])
        @test_throws DomainError QBase.QubitKet([1,0,0])
    end
end

@testset "is_density_matrix()" begin
    @testset "invalid matrices" begin
        @test !(QBase.is_density_matrix([1 0;0 1]))
        @test !(QBase.is_density_matrix([0.5 1;-1 0.5]))
        @test !(QBase.is_density_matrix([0.5 1;1 0.5]))
        @test !(QBase.is_density_matrix([1.5 0;0 -0.5]))
        @test !(QBase.is_density_matrix([0.5 im;-im 0.5]))
        @test !(QBase.is_density_matrix([1 4;4 1]))

        @test_throws DomainError QBase.is_density_matrix([0.5 0;0 0.5; 0 0])
    end
    @testset "valid matrices" begin
        @test QBase.is_density_matrix([0.5 0;0 0.5]::Matrix{Float64})
        @test QBase.is_density_matrix([1 0;0 0]::Matrix{Int64})
        @test QBase.is_density_matrix([0 0;0 1])
        @test QBase.is_density_matrix([0.5 0.5;0.5 0.5])
        @test QBase.is_density_matrix([0.5 0.5im;-0.5im 0.5]::Matrix{Complex{Float64}})
    end
end

@testset "QBase.DensityMatrix()" begin
    valid_cases = [
        [1 0 0;0 0 0;0 0 0]::Matrix{Int64},
        [0.5 0.5;0.5 0.5]::Matrix{Float64},
        [0.2 0 0;0 0.6 -0.1;0 -0.1 0.2],
        [0.5 0.5im;-0.5im 0.5]::Matrix{Complex{Float64}}
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ρ = QBase.DensityMatrix(valid_cases[id])
        @test ρ isa QBase.DensityMatrix
        @test ρ isa QBase.States.AbstractDensityMatrix
        @test ρ == valid_cases[id]
    end

    invalid_cases = [
        [1 0;0 1],
        [0.7 0.5;0.5 0.3],
        [0.5 0.5;-0.5 0.5],
        [0.5 0.5im;0.5im 0.5]
    ]

    @testset "invalid input: $id/$(length(invalid_cases))" for id in 1:length(invalid_cases)
        ρ = invalid_cases[id]
        @test_throws DomainError QBase.DensityMatrix(ρ)
    end
end

@testset "QBase.Qubit()" begin
    valid_cases = [
        [0.5 0.5;0.5 0.5]::Matrix{Float64},
        [1 0;0 0],
        [0.5 0.5im;-0.5im 0.5]::Matrix{Complex{Float64}}
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ρ = QBase.Qubit(valid_cases[id])
        @test ρ isa QBase.Qubit
        @test ρ isa QBase.States.AbstractDensityMatrix
        @test ρ == valid_cases[id]
    end

    invalid_cases = [
        [1 0;0 1],
        [1 0 0;0 0 0;0 0 0],
        [0.5 0.5;-0.5 0.5],
        [0.5 0.5im;0.5im 0.5]
    ]

    @testset "invalid input: $id/$(length(invalid_cases))" for id in 1:length(invalid_cases)
        ρ = invalid_cases[id]
        @test_throws DomainError QBase.Qubit(ρ)
    end
end

@testset "pure_state()" begin
    valid_cases = [
        (QBase.Ket([1,0,0]), [1 0 0;0 0 0;0 0 0]),
        ([0;0;1], [0 0 0;0 0 0;0 0 1]),
        ([1;1]/sqrt(2), [0.5 0.5;0.5 0.5])
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ρ = QBase.pure_state(valid_cases[id][1])
        @test ρ isa QBase.DensityMatrix
        @test rank(ρ) == 1
        @test ρ ≈ valid_cases[id][2]
    end

    @test_throws DomainError QBase.pure_state([1;1])
end

@testset "pure_qubit()" begin
    valid_cases = [
        (QBase.Ket([1,0]), [1 0;0 0]),
        ([0;1], [0 0;0 1]),
        ([1;1]/sqrt(2), [0.5 0.5;0.5 0.5])
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ρ = QBase.pure_qubit(valid_cases[id][1])
        @test ρ isa QBase.Qubit
        @test rank(ρ) == 1
        @test ρ ≈ valid_cases[id][2]
    end

    @test_throws DomainError QBase.pure_qubit([1;1])
end

@testset "mixed_state()" begin
    @test QBase.mixed_state(QBase.QMath.Marginals([0.7,0.2,0.1]), QBase.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 0 0;0 0 1],[0 0 0;0 1 0;0 0 0]])) == [0.7 0 0;0 0.1 0;0 0 0.2]
    @test QBase.mixed_state(QBase.QMath.Marginals([0.7,0.2,0.1]), QBase.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 0 0;0 0 1],[0 0 0;0 1 0;0 0 0]])) isa QBase.DensityMatrix
    @test rank(QBase.mixed_state(QBase.QMath.Marginals([0.7,0.2,0.1]), QBase.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 0 0;0 0 1],[0 0 0;0 1 0;0 0 0]]))) == 3

    @test QBase.mixed_state(QBase.QMath.Marginals([0.7,0.3]), QBase.Qubit.([[1 0;0 0],[0 0;0 1]])) == [0.7 0;0 0.3]
    @test QBase.mixed_state(QBase.QMath.Marginals([0.7,0.3]), QBase.Qubit.([[1 0;0 0],[0 0;0 1]])) isa QBase.DensityMatrix
    @test rank(QBase.mixed_state(QBase.QMath.Marginals([0.7,0.3]), QBase.Qubit.([[1 0;0 0],[0 0;0 1]]))) == 2

    @test_throws DomainError QBase.mixed_state(QBase.QMath.Marginals([0.7,0.2]), QBase.DensityMatrix.([[1 0;0 0],[0 0;0 1]]))
end

@testset "mixed_qubit()" begin
    @test QBase.mixed_qubit(QBase.QMath.Marginals([0.7,0.3]), QBase.Qubit.([[1 0;0 0],[0 0;0 1]])) == [0.7 0;0 0.3]
    @test QBase.mixed_qubit(QBase.QMath.Marginals([0.7,0.3]), QBase.Qubit.([[1 0;0 0],[0 0;0 1]])) isa QBase.Qubit
    @test rank(QBase.mixed_qubit(QBase.QMath.Marginals([0.7,0.3]), QBase.Qubit.([[1 0;0 0],[0 0;0 1]]))) == 2

    @test QBase.mixed_qubit(QBase.QMath.Marginals([0.7,0.3]), QBase.Qubit.([[1 0;0 0],[0 0;0 1]])) == [0.7 0;0 0.3]
    @test QBase.mixed_qubit(QBase.QMath.Marginals([0.7,0.3]), QBase.Qubit.([[1 0;0 0],[0 0;0 1]])) isa QBase.Qubit
    @test rank(QBase.mixed_qubit(QBase.QMath.Marginals([0.7,0.3]), QBase.Qubit.([[1 0;0 0],[0 0;0 1]]))) == 2

    @test_throws DomainError QBase.mixed_qubit(QBase.QMath.Marginals([0.7,0.2]), QBase.DensityMatrix.([[1 0;0 0],[0 0;0 1]]))
end

@testset "trine_qubit_kets" begin
    @test isapprox(QBase.trine_qubit_kets, QBase.mirror_symmetric_qubit_kets(π/3))
    @test QBase.trine_qubit_kets isa Vector{QBase.QubitKet}
end

@testset "mirror_symmetric_qubit_kets()" begin
    @test QBase.mirror_symmetric_qubit_kets(0) isa Vector{QBase.QubitKet}

    @test isapprox(QBase.mirror_symmetric_qubit_kets(0), [[1.0;0],[1.0;0],[1.0;0]])
    @test isapprox(QBase.mirror_symmetric_qubit_kets(π/4), [[1.0;0],[1;1]/sqrt(2),[1;-1]/sqrt(2)])
    @test isapprox(QBase.mirror_symmetric_qubit_kets(π/3),[[1.0 0]',[0.5 (sqrt(3)/2)]',[0.5 -sqrt(3)/2]'] )
    @test isapprox(QBase.mirror_symmetric_qubit_kets(π/2), [[1;0],[0;1],[0;-1]])
end

@testset "trine_qubits" begin
    @test isapprox(QBase.trine_qubits, QBase.mirror_symmetric_qubits(π/3))
    @test QBase.trine_qubits isa Vector{QBase.Qubit}
end

@testset "mirror_symmetric_qubits()" begin
    @testset "is density matrix: θ = $θ" for θ in 0:0.1π:π/2
        @test QBase.mirror_symmetric_qubits(θ) isa Vector{QBase.Qubit}
    end

    # if we pass in a bloch angle rather than hilbert angle we get a DomainError
    @test_throws DomainError QBase.mirror_symmetric_qubits(π)

    @test QBase.mirror_symmetric_qubits(π/2) ≈ [[1 0; 0 0], [0 0; 0 1], [0 0; 0 1]]
end

@testset "bloch_qubit_ket()" begin
    ψ = QBase.bloch_qubit_ket(π/2,π/2)

    @test ψ isa QBase.QubitKet
    @test norm(ψ) == 1
    @test length(ψ) == 2

    @test isapprox(ψ[1], 1/sqrt(2))
    @test isapprox(ψ[2], im/sqrt(2))

    @testset "valid wavefunction" for θ in 0:0.2:π
        @test all(map(
            (ϕ) -> QBase.bloch_qubit_ket(θ,ϕ) isa QBase.QubitKet,
            0:0.4:2*π
        ))
    end
end

@testset "bloch_qubit()" begin
    @test [1 0;0 0] == QBase.bloch_qubit(0,0,1)
    @test 0.5*[1 -im;im 1] == QBase.bloch_qubit(0,1,0)
    @test 0.5*[1 1;1 1] == QBase.bloch_qubit(1,0,0)
    @test 0.5*[1 0;0 1] == QBase.bloch_qubit(0,0,0)

    @test QBase.bloch_qubit(0,0,1) == QBase.bloch_qubit(0,0)
    @test QBase.bloch_qubit(0,1,0) ≈ QBase.bloch_qubit(π/2,π/2)
    @test QBase.bloch_qubit(1,0,0) ≈ QBase.bloch_qubit(π/2,0)

    @test QBase.bloch_qubit(0,0,1) isa QBase.Qubit
end

@testset "sic_qubits" begin
    @test QBase.sic_qubits isa Vector{QBase.Qubit}
end

@testset "bb84_qubits" begin
     @test QBase.bb84_qubits isa Vector{QBase.Qubit}
end

end
