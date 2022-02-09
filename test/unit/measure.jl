using Test
using QBase

@testset "./src/measure.jl" begin

@testset "measure(::POVM)" begin
    @testset "State" begin
        ρ = State([0.5 0.5im;-0.5im 0.5])
        Π = POVM([[1 0;0 0],[0 0;0 1]])

        marginals = measure(Π, ρ)

        @test marginals isa Probabilities
        @test marginals == [0.5,0.5]
    end

    @testset "Ket" begin
        ψ = Ket(1/sqrt(2)*[1,im])
        Π = POVM([[0.5 0.5im;-0.5im 0.5],[0.5 -0.5im;0.5im 0.5]])

        marginals = measure(Π, ψ)

        @test marginals isa Probabilities
        @test marginals ≈ [0,1]
    end

    @testset "Vector{State}" begin
        @testset "trine_measurements" begin
           conditionals = measure(
               mirror_symmetric_qubit_3povm(π/3), trine_qubit_states()
           )
           @test conditionals ≈ [2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3]
           @test conditionals isa Conditionals
       end

       @testset "classical_measurements" begin
           Π = POVM([[1 0;0 0],[0 0;0 1],[0 0;0 0]])
           ρ_set = State.([[0.5 0.5;0.5 0.5],[1 0;0 0],[0 0;0 1]])

           @test measure(Π, ρ_set) == [0.5 1 0;0.5 0 1;0 0 0]

           Π = POVM([[1 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]])
           ρ_set = State.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]])

           @test measure(Π, ρ_set) == [1 1 0;0 0 1]
       end
    end

    @testset "Vector{Ket}" begin
        ψ_kets = Ket.([[1,0],[0,1]])
        Π = POVM([[0 0;0 1],[1 0;0 0],[0 0;0 0]])

        conditionals = measure(Π, ψ_kets)

        @test conditionals isa Conditionals
        @test conditionals == [0 1;1 0;0 0]
    end
end

@testset "measure(::PVM)" begin
    @testset "States" begin
        ρ_ensemble = State.([[1 0;0 0],[1 -1im;1im 1]/2])
        Π = PVM([[1,0],[0,1]])

        probs = measure(Π, ρ_ensemble[1])
        @test probs isa Probabilities
        @test probs == [1,0]

        conditionals = measure(Π, ρ_ensemble)
        @test conditionals isa Conditionals
        @test conditionals ≈ [1 0.5;0 0.5]
    end

    @testset "Kets" begin
        ψ_ensemble = Ket.([[1,0],[1,1im]/sqrt(2)])
        Π = PVM([[1,0],[0,1]])

        probs = measure(Π, ψ_ensemble[1])
        @test probs isa Probabilities
        @test probs == [1,0]

        conditionals = measure(Π, ψ_ensemble)
        @test conditionals isa Conditionals
        @test conditionals ≈ [1 0.5;0 0.5]
    end
end

end
