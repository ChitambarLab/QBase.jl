using Test

@testset "./src/QBase.jl - Measurement functionality" begin

using QBase

@testset "measurement_probs()" begin
    @testset "DensityMatrix" begin
        ρ = States.DensityMatrix([0.5 0.5im;-0.5im 0.5])
        Π = Observables.POVM([[1 0;0 0],[0 0;0 1]])

        marginals = measurement_probs(Π, ρ)

        @test marginals isa QMath.Marginals
        @test marginals == [0.5,0.5]
    end

    @testset "Ket" begin
        ψ = States.Ket(1/sqrt(2)*[1,im])
        Π = Observables.POVM([[0.5 0.5im;-0.5im 0.5],[0.5 -0.5im;0.5im 0.5]])

        marginals = measurement_probs(Π, ψ)

        @test marginals isa QMath.Marginals
        @test marginals ≈ [0,1]
    end

    @testset "Vector{DensityMatrix}" begin
        @testset "trine_measurements" begin
           conditionals = measurement_probs(
               Observables.mirror_symmetric_qubit_3povm(π/3), States.trine_qubits
           )
           @test conditionals ≈ [2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3]
           @test conditionals isa QMath.Conditionals
       end

       @testset "classical_measurements" begin
           Π = Observables.QubitPOVM([[1 0;0 0],[0 0;0 1],[0 0;0 0]])
           ρ_set = States.Qubit.([[0.5 0.5;0.5 0.5],[1 0;0 0],[0 0;0 1]])

           @test measurement_probs(Π, ρ_set) == [0.5 1 0;0.5 0 1;0 0 0]

           Π = Observables.POVM([[1 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]])
           ρ_set = States.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]])

           @test measurement_probs(Π, ρ_set) == [1 1 0;0 0 1]
       end
    end

    @testset "Vector{Ket}" begin
        ψ_kets = States.QubitKet.([[1,0],[0,1]])
        Π = Observables.POVM([[0 0;0 1],[1 0;0 0],[0 0;0 0]])

        conditionals = measurement_probs(Π, ψ_kets)

        @test conditionals isa QMath.Conditionals
        @test conditionals == [0 1;1 0;0 0]
    end
end

end
