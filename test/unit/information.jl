using Test, LinearAlgebra
using QBase

@testset "/src/information.jl" begin

@testset "shannon_entropy()" begin
    @test shannon_entropy(Probabilities([1,0,0,0])) == 0
    @test shannon_entropy([1,0,0,0]) == 0

    @test shannon_entropy(Probabilities(fill(1/8,8))) == 3
    @test shannon_entropy(fill(1/8,8)) == 3

    @test shannon_entropy(Probabilities([1/3,1/3,1/3])) == log2(3)
    @test shannon_entropy([1/3,1/3,1/3]) == log2(3)

    @test shannon_entropy(Probabilities([0,0.5,0.5])) == 1
    @test shannon_entropy(Probabilities([0.7,0.3])) ≈ 0.88129089
    @test shannon_entropy(Probabilities([1/2,1/4,1/8,1/8])) == 7/4

    @test_throws DomainError shannon_entropy([-0.1,1.1])
end

@testset "von_neumann_entropy()" begin
    @test von_neumann_entropy(State(diagm(0 => fill(1/4,4)))) == 2
    @test von_neumann_entropy(State([1 0;0 0])) == 0

    @test von_neumann_entropy([0.5 0.5;0.5 0.5]) == 0

    @test von_neumann_entropy.(State.([[1 0;0 0],[0.5 0;0 0.5]])) == [0,1]

    @test_throws DomainError von_neumann_entropy([1 0;0 1])
end

@testset "holevo_bound()" begin
    priors = Probabilities([0.5,0.5])
    ρ_set = State.([[1 0;0 0],[0 0;0 1]])
    @test holevo_bound(priors,ρ_set) == 1
    @test holevo_bound([0.5,0.5],[[1 0;0 0],[0 0;0 1]]) == 1
    @test holevo_bound([0.5,0.5], ρ_set) == 1

    priors = Probabilities([1/3,1/3,1/3])
    ρ_set = mirror_symmetric_qubit_states(0)
    @test holevo_bound(priors,ρ_set) ≈ 0

    priors = Probabilities([1/3,1/3,1/3])
    ρ_set = mirror_symmetric_qubit_states(π/2)
    @test isapprox(holevo_bound(priors,ρ_set), 2/3*log2(3/2)+1/3*log2(3),atol=1e-6)

    priors = Probabilities([1/3,1/3,1/3])
    ρ_set = mirror_symmetric_qubit_states(π/3)
    @test holevo_bound(priors,ρ_set) ≈ 1

    priors = Probabilities([1/3,1/3,1/3])
    ρ_set = mirror_symmetric_qubit_states(5π/12)
    @test isapprox(holevo_bound(priors,ρ_set), 0.957,atol=1e-3)

    qb1 = bloch_qubit_ket(2*π/3,π/2)
    qb2 = bloch_qubit_ket(2*π/3,3π/2)
    qb3 = bloch_qubit_ket(0,0)
    ρ_set = pure_state.([qb1,qb2,qb3])
    @test holevo_bound(priors,ρ_set) ≈ 1

    qb1 = bloch_qubit_state(2*π/3,0)
    qb2 = bloch_qubit_state(5*π/12,π)
    qb3 = bloch_qubit_state(0,0)
    @test holevo_bound(priors,[qb1,qb2,qb3]) ≈ 0.952526321

    states = [State([1 0;0 0]), State([1 0;0 0]), State([0 0;0 1])]
    @test holevo_bound(priors, states) ≈ 0.91829583

    priors = Probabilities(ones(4)/4)
    @test holevo_bound(priors, sic_qubit_states()) == 1
end

@testset "holevo_information()" begin
    priors = Probabilities([1/3,1/3,1/3])
    ρ_set = trine_qubit_states()
    povm = sqrt_povm(priors, ρ_set)
    @test holevo_information(priors, ρ_set, povm) ≈ 1/3

    @test holevo_information([0.5,0.5],[[1 0;0 0],[0 0;0 1]],[[1 0;0 0],[0 0;0 1]]) == 1

    rot_y = qubit_rotation(π,axis="y")

    rot_povm = POVM(map(el -> rot_y*el*rot_y', povm))

    @test holevo_information(priors, ρ_set, rot_povm) ≈ 0.5849625

    orth_povm = POVM([[0.5 0;0 0],[0.5 0;0 0],[0 0;0 1]])
    @test holevo_information(priors, ρ_set, orth_povm) ≈ 0.459147917
end

@testset "joint_entropy()" begin
    priors = Probabilities([0.5,0.5])
    conditionals = Conditionals([0.5 0.5;0.5 0.5])

    @test joint_entropy(priors, conditionals) == 2
    @test joint_entropy(priors, [0.5 0.5;0.5 0.5]) == 2

    priors = Probabilities([1/3,1/3,1/3])
    conditionals = Conditionals([2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3])
    @test joint_entropy(priors, conditionals) ≈ 2.8365917
end

@testset "mutual_information()" begin
    priors = Probabilities([1/3,1/3,1/3])
    conditionals = Conditionals([2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3])
    @test mutual_information(priors,conditionals) ≈ 1/3
    @test mutual_information([1/3,1/3,1/3],conditionals) ≈ 1/3

    priors = Probabilities([1/3,1/3,1/3])
    conditionals = Conditionals([0 1/2 1/2;1/2 0 1/2;1/2 1/2 0])
    @test mutual_information(priors,conditionals) ≈ 0.5849625

    conditionals = Conditionals([1 0 1;0 1 0;0 0 0])
    @test mutual_information(priors,conditionals) ≈ 2/3*log2(3/2)+1/3*log2(3)

    conditionals = Conditionals([1 0 0;0 1 0;0 0 1])
    @test mutual_information(priors,conditionals) ≈ log2(3)

    conditionals = Conditionals([1/2 1/8 1/8;1/2 1/8 1/8;0 3/4 3/4])
    @test mutual_information(priors,conditionals) ≈ 0.459147917
end

@testset "conditional_entropy()" begin
    priors = Probabilities([1/3, 1/3, 1/3])
    conditionals = Conditionals([2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3])
    @test conditional_entropy(priors, conditionals) ≈ 1.251629167

    cond_ent = conditional_entropy(priors, conditionals)
    mut_inf = mutual_information(priors,conditionals)
    @test cond_ent + mut_inf ≈ shannon_entropy(priors)

    priors = Probabilities([1/3, 1/3, 1/3])
    conditionals = Conditionals([3/4 1/8 1/8;1/8 3/4 1/8;1/8 1/8 3/4])
    conditional_entropy(priors, conditionals)

    cond_ent = conditional_entropy(priors, conditionals)
    mut_inf = mutual_information(priors,conditionals)
    @test cond_ent + mut_inf ≈ shannon_entropy(priors)
end

@testset "success_probability()" begin
    @testset "trivial case" begin
        Π = POVM([[1 0;0 0],[0 0;0 1]])
        ρ_states = State.([[1 0;0 0],[0 0;0 1]])
        priors = Probabilities([1/2,1/2])

        @test success_probability(priors, ρ_states, Π) == 1
        @test success_probability([1/2,1/2], [[1 0;0 0],[0 0;0 1]], Π) == 1
    end
end

@testset "error_probability()" begin
    @testset "trivial case" begin
        Π = POVM([[1 0;0 0],[0 0;0 1]])
        ρ_states = State.([[1 0;0 0],[0 0;0 1]])
        priors = Probabilities([1/2,1/2])

        @test error_probability(priors, ρ_states, Π) == 0
    end
end

end
