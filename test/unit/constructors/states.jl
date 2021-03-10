using Test, LinearAlgebra
using QBase

@testset "./src/constructors/states.jl" begin

@testset "pure_state()" begin
    @testset "valid inputs" begin
        ρ = pure_state(Ket([im,0]))
        @test ρ isa State{Complex{Int64}}
        @test rank(ρ) == 1
        @test ρ == [1 0;0 0]

        ρ = pure_state(Bra([0.,-1im]))
        @test ρ isa State{Complex{Float64}}
        @test rank(ρ) == 1
        @test ρ == [0 0;0 1]

        ρ = pure_state([1,1]/sqrt(2))
        @test ρ isa State{Float64}
        @test rank(ρ) == 1
        @test ρ ≈ [1 1;1 1]/2

        ρ = pure_state(([1,-im]/sqrt(2))')
        @test ρ isa State{Complex{Float64}}
        @test rank(ρ) == 1
        @test ρ ≈ [1 im;-im 1]/2
    end

    @test_throws DomainError pure_state([1;1])
end

@testset "mixed_state()" begin
    ρ_mix = mixed_state(
        Marginals([0.7,0.2,0.1]),
        State.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 0 0;0 0 1],[0 0 0;0 1 0;0 0 0]])
    )
    @test ρ_mix == [0.7 0 0;0 0.1 0;0 0 0.2]
    @test ρ_mix isa State{Float64}
    @test rank(ρ_mix) == 3

    ρ_mix = mixed_state([0.7,0.2,0.1],[[1 0 0;0 0 0;0 0 0],[0 0 0;0 0 0;0 0 1],[0 0 0;0 1 0;0 0 0]])
    @test ρ_mix == [0.7 0 0;0 0.1 0;0 0 0.2]
    @test ρ_mix isa State{Float64}
    @test rank(ρ_mix) == 3

    @test_throws DomainError mixed_state([0.7,0.2], [[1 0;0 0],[0 0;0 1]])
    @test_throws DomainError mixed_state([0.7,0.3], [[1 0;0 .1],[0 0;0 1]])
end

@testset "trine_qubit_states()" begin
    trine_states = trine_qubit_states()
    @test isapprox(trine_states, mirror_symmetric_qubit_states(π/3), atol=1e-7)
    @test trine_states isa Vector{State{Float64}}
end

@testset "mirror_symmetric_qubit_states()" begin
    @testset "range over θ" begin
        for θ in 0:0.1π:π/2
            @test mirror_symmetric_qubit_states(θ) isa Vector{State{Float64}}
        end
    end

    @test_throws DomainError mirror_symmetric_qubit_states(π)

    @test isapprox(mirror_symmetric_qubit_states(π/2),[[1 0; 0 0], [0 0; 0 1], [0 0; 0 1]], atol=1e-7)
end

@testset "bloch_qubit_state()" begin
    @test [1 0;0 0] == bloch_qubit_state(0,0,1)
    @test 0.5*[1 -im;im 1] == bloch_qubit_state(0,1,0)
    @test 0.5*[1 1;1 1] == bloch_qubit_state(1,0,0)
    @test 0.5*[1 0;0 1] == bloch_qubit_state(0,0,0)

    @test bloch_qubit_state(0,0,1) == bloch_qubit_state(0,0)
    @test isapprox(bloch_qubit_state(0,1,0),bloch_qubit_state(π/2,π/2),atol=1e-7)
    @test isapprox(bloch_qubit_state(1,0,0),bloch_qubit_state(π/2,0),atol=1e-7)

    @test bloch_qubit_state(0,0,1) isa State{Complex{Float64}}
    @test_throws DomainError bloch_qubit_state(1,1,1)
end

@testset "sic_qubit_states()" begin
    sic_states = sic_qubit_states()

    @test sic_states isa Vector{State{Complex{Float64}}}

    θ_sic = 2*acos(1/sqrt(3))

    @test sic_states[1] ≈ bloch_qubit_state(0, 0)
    @test sic_states[2] ≈ bloch_qubit_state(θ_sic, 0)
    @test sic_states[3] ≈ bloch_qubit_state(θ_sic, 2π/3)
    @test sic_states[4] ≈ bloch_qubit_state(θ_sic, 4π/3)
end

@testset "bb84_qubit_states()" begin
    bb84_states = bb84_qubit_states()
    @test bb84_states isa Vector{State{Float64}}

    @test isapprox(bb84_states, planar_symmetric_qubit_states(4), atol=1e-7)
end

@testset "computational_basis_states()" begin
    @test computational_basis_states(3) == [[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]]
    @test computational_basis_states(2) == [[1 0;0 0],[0 0;0 1]]
end

@testset "bell_states()" begin
    states = bell_states()

    @test states isa Vector{State{Float64}}
    @test states[1] ≈ 1/2*[1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1]
    @test states[2] ≈ 1/2*[1 0 0 -1;0 0 0 0;0 0 0 0;-1 0 0 1]
    @test states[3] ≈ 1/2*[0 0 0 0;0 1 1 0;0 1 1 0;0 0 0 0]
    @test states[4] ≈ 1/2*[0 0 0 0;0 1 -1 0;0 -1 1 0;0 0 0 0]
end

@testset "generalized_bell_states()" begin
    @testset "qubit bell states" begin
        states = generalized_bell_states(2)

        @test states isa Vector{State{Complex{Float64}}}
        @test isapprox(states, bell_states(),atol=1e-7)
    end
end

@testset "planar_symmetric_qubit_states()" begin
    @test isapprox(planar_symmetric_qubit_states(2), [[1 0;0 0],[0 0;0 1]], atol=1e-7)
    @test isapprox(planar_symmetric_qubit_states(3), trine_qubit_states(),atol=1e-7)
        
    @testset "scanning over cases verifying right number are computed" begin
        for n in 5:50
            qubits = planar_symmetric_qubit_states(n)

            @test length(qubits) == n
            @test qubits isa Vector{State{Float64}}
        end
    end
end

end
