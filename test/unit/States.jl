using Test, LinearAlgebra

@testset "QBase/States.jl" begin

using QBase

@testset "States.is_ket()" begin
    @testset "valid wavefunctions" begin
        @test States.is_ket([1;0])
        @test States.is_ket([1])
        @test States.is_ket(1/sqrt(3)*[1;1;-1])
        @test States.is_ket(1/sqrt(2)*[1;im])
    end

    @testset "invalid wavefunctions" begin
        @test !( States.is_ket([0;0]) )
        @test !( States.is_ket([1;1]) )
    end

    @testset "valid inputs" begin
        @test States.is_ket([1,0]::Vector{Int64})
        @test States.is_ket([1.0,0.0]::Vector{Float64})
        @test States.is_ket([1.0+im*0.0,0.0+im*0.0 ]::Vector{Complex{Float64}})
    end

    @testset "invalid inputs" begin
        @test_throws MethodError States.is_ket([1 0;0 0])
        @test_throws MethodError States.is_ket(1/sqrt(2)*[1 1])
    end

    @testset "unitary evolution" begin
        ψ = Unitaries.σx*States.Ket([1,0])

        @test ψ == [0,1]
        @test ψ isa States.Ket
    end
end

@testset "States.Ket()" begin
    valid_cases = [
        [1;0;0],
        [1,0,0,0],
        [1,1,im]./sqrt(3),
        [0.5im,0,-sqrt(3)/2]
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ψ = States.Ket(valid_cases[id])
        @test ψ isa States.Ket
        @test ψ isa States.States.AbstractKet
        @test ψ == valid_cases[id]
    end

    @testset "invalid inputs" begin
        @test_throws DomainError States.Ket([1,1,1])
        @test_throws MethodError States.Ket([1 0])
    end
end

@testset "States.QubitKet()" begin
    valid_cases = [
        [1,0],
        [1,im]/sqrt(2),
        [0.5;sqrt(3)/2]
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ψ = States.QubitKet(valid_cases[id])
        @test ψ isa States.QubitKet
        @test ψ isa States.States.AbstractKet
        @test ψ == valid_cases[id]
    end

    @testset "invalid inputs" begin
        @test_throws DomainError States.QubitKet([1 0])
        @test_throws DomainError States.QubitKet([1,1])
        @test_throws DomainError States.QubitKet([1,0,0])
    end
end

@testset "is_density_matrix()" begin
    @testset "invalid matrices" begin
        @test !(States.is_density_matrix([1 0;0 1]))
        @test !(States.is_density_matrix([0.5 1;-1 0.5]))
        @test !(States.is_density_matrix([0.5 1;1 0.5]))
        @test !(States.is_density_matrix([1.5 0;0 -0.5]))
        @test !(States.is_density_matrix([0.5 im;-im 0.5]))
        @test !(States.is_density_matrix([1 4;4 1]))

        @test_throws DomainError States.is_density_matrix([0.5 0;0 0.5; 0 0])
    end
    @testset "valid matrices" begin
        @test States.is_density_matrix([0.5 0;0 0.5]::Matrix{Float64})
        @test States.is_density_matrix([1 0;0 0]::Matrix{Int64})
        @test States.is_density_matrix([0 0;0 1])
        @test States.is_density_matrix([0.5 0.5;0.5 0.5])
        @test States.is_density_matrix([0.5 0.5im;-0.5im 0.5]::Matrix{Complex{Float64}})
    end
end

@testset "States.DensityMatrix()" begin
    valid_cases = [
        [1 0 0;0 0 0;0 0 0]::Matrix{Int64},
        [0.5 0.5;0.5 0.5]::Matrix{Float64},
        [0.2 0 0;0 0.6 -0.1;0 -0.1 0.2],
        [0.5 0.5im;-0.5im 0.5]::Matrix{Complex{Float64}}
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ρ = States.DensityMatrix(valid_cases[id])
        @test ρ isa States.DensityMatrix
        @test ρ isa States.States.AbstractDensityMatrix
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
        @test_throws DomainError States.DensityMatrix(ρ)
    end

    @testset "kron" begin
        q0 = States.DensityMatrix([1 0;0 0])
        q_plus = States.DensityMatrix(1/2*[1 1;1 1])

        @test q0 isa States.DensityMatrix
        @test q_plus isa States.DensityMatrix

        @test kron(q0,q_plus) isa States.DensityMatrix
    end

    @testset "partial_trace" begin
        ρ = States.DensityMatrix(0.5*[1 0 0 0;0 0 0 0;0 0 1 0;0 0 0 0])

        ρ_sub1 = QMath.partial_trace(ρ, [2,2], 1)

        @test ρ_sub1 isa States.DensityMatrix
        @test ρ_sub1 == [1 0;0 0]

        ρ_sub2 = QMath.partial_trace(ρ, [2,2], 2)

        @test ρ_sub2 isa States.DensityMatrix
        @test ρ_sub2 == 0.5*[1 0;0 1]
    end
end

@testset "States.Qubit()" begin
    valid_cases = [
        [0.5 0.5;0.5 0.5]::Matrix{Float64},
        [1 0;0 0],
        [0.5 0.5im;-0.5im 0.5]::Matrix{Complex{Float64}}
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ρ = States.Qubit(valid_cases[id])
        @test ρ isa States.Qubit
        @test ρ isa States.States.AbstractDensityMatrix
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
        @test_throws DomainError States.Qubit(ρ)
    end

    @testset "kron" begin
        q0 = States.Qubit([1 0;0 0])
        q_plus = States.Qubit(1/2*[1 1;1 1])

        @test q0 isa States.Qubit
        @test q_plus isa States.Qubit

        @test kron(q0,q_plus) isa States.DensityMatrix
    end
end

@testset "pure_state()" begin
    valid_cases = [
        (States.Ket([1,0,0]), [1 0 0;0 0 0;0 0 0]),
        ([0;0;1], [0 0 0;0 0 0;0 0 1]),
        ([1;1]/sqrt(2), [0.5 0.5;0.5 0.5])
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ρ = States.pure_state(valid_cases[id][1])
        @test ρ isa States.DensityMatrix
        @test rank(ρ) == 1
        @test ρ ≈ valid_cases[id][2]
    end

    @test_throws DomainError States.pure_state([1;1])
end

@testset "pure_qubit()" begin
    valid_cases = [
        (States.Ket([1,0]), [1 0;0 0]),
        ([0;1], [0 0;0 1]),
        ([1;1]/sqrt(2), [0.5 0.5;0.5 0.5])
    ]

    @testset "valid input: $id/$(length(valid_cases))" for id in 1:length(valid_cases)
        ρ = States.pure_qubit(valid_cases[id][1])
        @test ρ isa States.Qubit
        @test rank(ρ) == 1
        @test ρ ≈ valid_cases[id][2]
    end

    @test_throws DomainError States.pure_qubit([1;1])
end

@testset "mixed_state()" begin
    @test States.mixed_state(States.QMath.Marginals([0.7,0.2,0.1]), States.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 0 0;0 0 1],[0 0 0;0 1 0;0 0 0]])) == [0.7 0 0;0 0.1 0;0 0 0.2]
    @test States.mixed_state(States.QMath.Marginals([0.7,0.2,0.1]), States.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 0 0;0 0 1],[0 0 0;0 1 0;0 0 0]])) isa States.DensityMatrix
    @test rank(States.mixed_state(States.QMath.Marginals([0.7,0.2,0.1]), States.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 0 0;0 0 1],[0 0 0;0 1 0;0 0 0]]))) == 3

    @test States.mixed_state(States.QMath.Marginals([0.7,0.3]), States.Qubit.([[1 0;0 0],[0 0;0 1]])) == [0.7 0;0 0.3]
    @test States.mixed_state(States.QMath.Marginals([0.7,0.3]), States.Qubit.([[1 0;0 0],[0 0;0 1]])) isa States.DensityMatrix
    @test rank(States.mixed_state(States.QMath.Marginals([0.7,0.3]), States.Qubit.([[1 0;0 0],[0 0;0 1]]))) == 2

    @test_throws DomainError States.mixed_state(States.QMath.Marginals([0.7,0.2]), States.DensityMatrix.([[1 0;0 0],[0 0;0 1]]))
end

@testset "mixed_qubit()" begin
    @test States.mixed_qubit(States.QMath.Marginals([0.7,0.3]), States.Qubit.([[1 0;0 0],[0 0;0 1]])) == [0.7 0;0 0.3]
    @test States.mixed_qubit(States.QMath.Marginals([0.7,0.3]), States.Qubit.([[1 0;0 0],[0 0;0 1]])) isa States.Qubit
    @test rank(States.mixed_qubit(States.QMath.Marginals([0.7,0.3]), States.Qubit.([[1 0;0 0],[0 0;0 1]]))) == 2

    @test States.mixed_qubit(States.QMath.Marginals([0.7,0.3]), States.Qubit.([[1 0;0 0],[0 0;0 1]])) == [0.7 0;0 0.3]
    @test States.mixed_qubit(States.QMath.Marginals([0.7,0.3]), States.Qubit.([[1 0;0 0],[0 0;0 1]])) isa States.Qubit
    @test rank(States.mixed_qubit(States.QMath.Marginals([0.7,0.3]), States.Qubit.([[1 0;0 0],[0 0;0 1]]))) == 2

    @test_throws DomainError States.mixed_qubit(States.QMath.Marginals([0.7,0.2]), States.DensityMatrix.([[1 0;0 0],[0 0;0 1]]))
end

@testset "trine_qubit_kets" begin
    @test isapprox(States.trine_qubit_kets, States.mirror_symmetric_qubit_kets(π/3))
    @test States.trine_qubit_kets isa Vector{States.QubitKet}
end

@testset "mirror_symmetric_qubit_kets()" begin
    @test States.mirror_symmetric_qubit_kets(0) isa Vector{States.QubitKet}

    @test isapprox(States.mirror_symmetric_qubit_kets(0), [[1.0;0],[1.0;0],[1.0;0]])
    @test isapprox(States.mirror_symmetric_qubit_kets(π/4), [[1.0;0],[1;1]/sqrt(2),[1;-1]/sqrt(2)])
    @test isapprox(States.mirror_symmetric_qubit_kets(π/3),[[1.0 0]',[0.5 (sqrt(3)/2)]',[0.5 -sqrt(3)/2]'] )
    @test isapprox(States.mirror_symmetric_qubit_kets(π/2), [[1;0],[0;1],[0;-1]])
end

@testset "trine_qubits" begin
    @test isapprox(States.trine_qubits, States.mirror_symmetric_qubits(π/3))
    @test States.trine_qubits isa Vector{States.Qubit}
end

@testset "mirror_symmetric_qubits()" begin
    @testset "is density matrix: θ = $θ" for θ in 0:0.1π:π/2
        @test States.mirror_symmetric_qubits(θ) isa Vector{States.Qubit}
    end

    # if we pass in a bloch angle rather than hilbert angle we get a DomainError
    @test_throws DomainError States.mirror_symmetric_qubits(π)

    @test States.mirror_symmetric_qubits(π/2) ≈ [[1 0; 0 0], [0 0; 0 1], [0 0; 0 1]]
end

@testset "bloch_qubit_ket()" begin
    ψ = States.bloch_qubit_ket(π/2,π/2)

    @test ψ isa States.QubitKet
    @test norm(ψ) == 1
    @test length(ψ) == 2

    @test isapprox(ψ[1], 1/sqrt(2))
    @test isapprox(ψ[2], im/sqrt(2))

    @testset "valid wavefunction" for θ in 0:0.2:π
        @test all(map(
            (ϕ) -> States.bloch_qubit_ket(θ,ϕ) isa States.QubitKet,
            0:0.4:2*π
        ))
    end
end

@testset "bloch_qubit()" begin
    @test [1 0;0 0] == States.bloch_qubit(0,0,1)
    @test 0.5*[1 -im;im 1] == States.bloch_qubit(0,1,0)
    @test 0.5*[1 1;1 1] == States.bloch_qubit(1,0,0)
    @test 0.5*[1 0;0 1] == States.bloch_qubit(0,0,0)

    @test States.bloch_qubit(0,0,1) == States.bloch_qubit(0,0)
    @test States.bloch_qubit(0,1,0) ≈ States.bloch_qubit(π/2,π/2)
    @test States.bloch_qubit(1,0,0) ≈ States.bloch_qubit(π/2,0)

    @test States.bloch_qubit(0,0,1) isa States.Qubit
end

@testset "sic_qubits" begin
    @test States.sic_qubits isa Vector{States.Qubit}
end

@testset "bb84_qubits" begin
     @test States.bb84_qubits isa Vector{States.Qubit}
end

@testset "States.basis_kets()" begin
    @test States.basis_kets(3) == [[1,0,0],[0,1,0],[0,0,1]]
    @test States.basis_kets(8) == [
        [1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]
    ]
end

@testset "States.basis_states()" begin
    @test States.basis_states(3) == [[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]]
    @test States.basis_states(2) == [[1 0;0 0],[0 0;0 1]]
end

@testset "States.bell_kets" begin
    kets = States.bell_kets
    @test length(kets) == 4
    @test kets isa Array{States.Ket,1}
    @test kets[1] == 1/sqrt(2)*(kron([1,0],[1,0])+kron([0,1],[0,1]))
    @test kets[2] == 1/sqrt(2)*(kron([1,0],[1,0])-kron([0,1],[0,1]))
    @test kets[3] == 1/sqrt(2)*(kron([1,0],[0,1])+kron([0,1],[1,0]))
    @test kets[4] == 1/sqrt(2)*(kron([1,0],[0,1])-kron([0,1],[1,0]))
end

@testset "States.bell_states" begin
    states = States.bell_states

    @test states isa Vector{States.DensityMatrix}
    @test states[1] ≈ 1/2*[1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1]
    @test states[2] ≈ 1/2*[1 0 0 -1;0 0 0 0;0 0 0 0;-1 0 0 1]
    @test states[3] ≈ 1/2*[0 0 0 0;0 1 1 0;0 1 1 0;0 0 0 0]
    @test states[4] ≈ 1/2*[0 0 0 0;0 1 -1 0;0 -1 1 0;0 0 0 0]
end

@testset "States.generalized_bell_kets()" begin
    @testset "qubit bell states" begin
        kets = States.generalized_bell_kets(2)

        @test kets isa Vector{States.Ket}
        @test kets ≈ States.bell_kets
    end

    @testset "qutrit bell states" begin
        kets = States.generalized_bell_kets(3)

        basis = States.basis_kets(3)
        ρ1 = (kron(basis[1],basis[1]) + kron(basis[2],basis[2]) + kron(basis[3],basis[3]))/sqrt(3)
        ρ2 = (kron(basis[1],basis[1]) + exp(im*2*π/3)*kron(basis[2],basis[2]) + exp(im*4*π/3)*kron(basis[3],basis[3]))/sqrt(3)
        ρ3 = (kron(basis[1],basis[1]) + exp(im*4*π/3)*kron(basis[2],basis[2]) + exp(im*2*π/3)*kron(basis[3],basis[3]))/sqrt(3)

        ρ4 = (kron(basis[1],basis[2]) + kron(basis[2],basis[3]) + kron(basis[3],basis[1]))/sqrt(3)
        ρ5 = (kron(basis[1],basis[2]) + exp(im*2*π/3)*kron(basis[2],basis[3]) + exp(im*4*π/3)*kron(basis[3],basis[1]))/sqrt(3)
        ρ6 = (kron(basis[1],basis[2]) + exp(im*4*π/3)*kron(basis[2],basis[3]) + exp(im*2*π/3)*kron(basis[3],basis[1]))/sqrt(3)

        ρ7 = (kron(basis[1],basis[3]) + kron(basis[2],basis[1]) + kron(basis[3],basis[2]))/sqrt(3)
        ρ8 = (kron(basis[1],basis[3]) + exp(im*2*π/3)*kron(basis[2],basis[1]) + exp(im*4*π/3)*kron(basis[3],basis[2]))/sqrt(3)
        ρ9 = (kron(basis[1],basis[3]) + exp(im*4*π/3)*kron(basis[2],basis[1]) + exp(im*2*π/3)*kron(basis[3],basis[2]))/sqrt(3)

        ψ_match = [ρ1,ρ2,ρ3,ρ4,ρ5,ρ6,ρ7,ρ8,ρ9]

        @test kets isa Vector{States.Ket}
        @test kets ≈ ψ_match
    end
end

@testset "generalized_bell_states()" begin
    @testset "qubit bell states" begin
        states = States.generalized_bell_states(2)

        @test states isa Vector{States.DensityMatrix}
        @test states ≈ States.bell_states
    end
end

@testset "planar_symmetric_qubit_kets()" begin
    @testset "simple cases" begin
        @testset "orthogonal qubits" begin
            kets = States.planar_symmetric_qubit_kets(2)

            @test kets[1] ≈ [1,0]
            @test kets[2] ≈ [0,1]
        end

        @testset "trine states" begin
            kets = States.planar_symmetric_qubit_kets(3)
            trine_kets = States.trine_qubit_kets

            @test kets[1] ≈ trine_kets[1]
            @test kets[2] ≈ trine_kets[2]
            @test kets[3] ≈ -1*trine_kets[3]
        end

        @testset "bb84 states" begin
            kets = States.planar_symmetric_qubit_kets(4)

            @test kets[1] ≈ [1,0]
            @test kets[2] ≈ [1,1]/sqrt(2)
            @test kets[3] ≈ [0,1]
            @test kets[4] ≈ [-1,1]/sqrt(2)
        end
    end

    @testset "scanning over cases verifying right number are computed" begin
        for n in 2:100
            kets = States.planar_symmetric_qubit_kets(n)

            @test length(kets) == n
            @test kets isa Vector{States.QubitKet}
        end
    end

    @test_throws(DomainError, States.planar_symmetric_qubit_kets(1))
end

@testset "planar_symmetric_qubits()" begin
    @testset "simple cases" begin
        @testset "bb84 cases" begin
            states = States.planar_symmetric_qubits(4)

            @test states[1] ≈ States.bb84_qubits[1]
            @test states[2] ≈ States.bb84_qubits[3]
            @test states[3] ≈ States.bb84_qubits[2]
            @test states[4] ≈ States.bb84_qubits[4]
        end
    end

    @testset "scanning over cases verifying right number are computed" begin
        for n in 2:100
            qubits = States.planar_symmetric_qubits(n)

            @test length(qubits) == n
            @test qubits isa Vector{States.Qubit}
        end
    end
end

end
