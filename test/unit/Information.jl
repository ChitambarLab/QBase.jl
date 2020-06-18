using Test, LinearAlgebra

@testset "/src/QBase/Information.jl" begin

using QBase

@testset "shannon_entropy()" begin
    @test Information.shannon_entropy(QMath.Marginals(fill(1/8,8))) == 3
    @test Information.shannon_entropy(QMath.Marginals([1/3,1/3,1/3])) == log2(3)
    @test Information.shannon_entropy(QMath.Marginals([0.5,0.5])) == 1
    @test Information.shannon_entropy(QMath.Marginals([0,0.5,0.5])) == 1
    @test Information.shannon_entropy(QMath.Marginals([0.7,0.3])) ≈ 0.88129089
    @test Information.shannon_entropy(QMath.Marginals([0,0,0,0,0,0,0,0,0,0,0,0,0,1])) == 0
    @test Information.shannon_entropy(QMath.Marginals([1/2,1/4,1/8,1/8])) == 7/4

    @test_throws DomainError Information.shannon_entropy(QMath.Marginals([2,0]))
end

@testset "von_neumann_entropy()" begin
    @test Information.von_neumann_entropy(States.DensityMatrix(diagm(0 => fill(1/8,8)))) == 3
    @test Information.von_neumann_entropy(States.Qubit([1 0;0 0])) == 0

    @test Information.von_neumann_entropy.(States.DensityMatrix.([[1 0;0 0],[0.5 0;0 0.5]])) == [0,1]

    @test_throws DomainError Information.von_neumann_entropy(States.DensityMatrix([1 0;0 1]))
end

@testset "holevo_bound()" begin

    priors = QMath.Marginals([0.5,0.5])
    ρ_set = States.DensityMatrix.([[1 0;0 0],[0 0;0 1]])
    @test Information.holevo_bound(priors,ρ_set) == 1

    priors = QMath.Marginals([1/3,1/3,1/3])
    ρ_set = States.mirror_symmetric_qubits(0)
    @test Information.holevo_bound(priors,ρ_set) ≈ 0

    priors = QMath.Marginals([1/3,1/3,1/3])
    ρ_set = States.mirror_symmetric_qubits(π/2)
    @test isapprox(Information.holevo_bound(priors,ρ_set), 2/3*log2(3/2)+1/3*log2(3),atol=1e-6)

    priors = QMath.Marginals([1/3,1/3,1/3])
    ρ_set = States.mirror_symmetric_qubits(π/3)
    @test Information.holevo_bound(priors,ρ_set) ≈ 1

    priors = QMath.Marginals([1/3,1/3,1/3])
    ρ_set = States.mirror_symmetric_qubits(5π/12)
    @test Information.holevo_bound(priors,ρ_set) ≈ 0.9566111446

    qb1 = States.bloch_qubit_ket(2*π/3,π/2)
    qb2 = States.bloch_qubit_ket(2*π/3,3π/2)
    qb3 = States.bloch_qubit_ket(0,0)

    ρ_set = States.pure_state.([qb1,qb2,qb3])
    @test Information.holevo_bound(priors,ρ_set) ≈ 1

    qb1 = States.bloch_qubit(2*π/3,0)
    qb2 = States.bloch_qubit(5*π/12,π)
    qb3 = States.bloch_qubit(0,0)

    @test Information.holevo_bound(priors,[qb1,qb2,qb3]) ≈ 0.952526321

    states = [States.Qubit([1 0;0 0]), States.Qubit([1 0;0 0]), States.Qubit([0 0;0 1])]
    @test Information.holevo_bound(priors, states) ≈ 0.91829583

    priors = QMath.Marginals(ones(4)/4)
    @test Information.holevo_bound(priors, States.sic_qubits) == 1
end

@testset "holevo_information()" begin
    priors = QMath.Marginals([1/3,1/3,1/3])
    ρ_set = States.trine_qubits
    povm = Observables.sqrt_povm(priors, ρ_set)
    @test Information.holevo_information(priors, ρ_set, povm) ≈ 1/3

    rot_y = Unitaries.qubit_rotation(π,axis="y")

    rot_povm = Observables.QubitPOVM(map(el -> rot_y*el*rot_y', povm))

    @test Information.holevo_information(priors, ρ_set, rot_povm) ≈ 0.5849625

    orth_povm = Observables.QubitPOVM([[0.5 0;0 0],[0.5 0;0 0],[0 0;0 1]])
    @test Information.holevo_information(priors, ρ_set, orth_povm) ≈ 0.459147917
end

@testset "joint_entropy()" begin
    priors = QMath.Marginals([0.5,0.5])
    conditionals = QMath.Conditionals([0.5 0.5;0.5 0.5])

    @test Information.joint_entropy(priors, conditionals) == 2

    priors = QMath.Marginals([1/3,1/3,1/3])
    conditionals = QMath.Conditionals([2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3])
    @test Information.joint_entropy(priors, conditionals) ≈ 2.8365917
end

@testset "mutual_information()" begin
    priors = QMath.Marginals([1/3,1/3,1/3])
    conditionals = QMath.Conditionals([2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3])
    @test Information.mutual_information(priors,conditionals) ≈ 1/3

    priors = QMath.Marginals([1/3,1/3,1/3])
    conditionals = QMath.Conditionals([0 1/2 1/2;1/2 0 1/2;1/2 1/2 0])
    @test Information.mutual_information(priors,conditionals) ≈ 0.5849625

    conditionals = QMath.Conditionals([1 0 1;0 1 0;0 0 0])
    @test Information.mutual_information(priors,conditionals) ≈ 2/3*log2(3/2)+1/3*log2(3)

    conditionals = QMath.Conditionals([1 0 0;0 1 0;0 0 1])
    @test Information.mutual_information(priors,conditionals) ≈ log2(3)

    conditionals = QMath.Conditionals([1/2 1/8 1/8;1/2 1/8 1/8;0 3/4 3/4])
    @test Information.mutual_information(priors,conditionals) ≈ 0.459147917
end

@testset "conditional_entropy()" begin
    priors = QMath.Marginals([1/3, 1/3, 1/3])
    conditionals = QMath.Conditionals([2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3])
    @test Information.conditional_entropy(priors, conditionals) ≈ 1.251629167

    cond_ent = Information.conditional_entropy(priors, conditionals)
    mut_inf = Information.mutual_information(priors,conditionals)
    @test cond_ent + mut_inf ≈ Information.shannon_entropy(priors)

    priors = QMath.Marginals([1/3, 1/3, 1/3])
    conditionals = QMath.Conditionals([3/4 1/8 1/8;1/8 3/4 1/8;1/8 1/8 3/4])
    Information.conditional_entropy(priors, conditionals)

    cond_ent = Information.conditional_entropy(priors, conditionals)
    mut_inf = Information.mutual_information(priors,conditionals)
    @test cond_ent + mut_inf ≈ Information.shannon_entropy(priors)
end

@testset "success_probability()" begin
    @testset "trivial case" begin
        Π = Observables.POVM([[1 0;0 0],[0 0;0 1]])
        ρ_states = States.DensityMatrix.([[1 0;0 0],[0 0;0 1]])
        priors = QMath.Marginals([1/2,1/2])

        @test Information.success_probability(priors, ρ_states, Π) == 1
    end
end

@testset "error_probability()" begin
    @testset "trivial case" begin
        Π = Observables.POVM([[1 0;0 0],[0 0;0 1]])
        ρ_states = States.DensityMatrix.([[1 0;0 0],[0 0;0 1]])
        priors = QMath.Marginals([1/2,1/2])

        @test Information.error_probability(priors, ρ_states, Π) == 0
    end
end

end
