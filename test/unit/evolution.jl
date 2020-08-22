using Test

@testset "./src/QBase.jl, evolution functionality" begin

using QBase

@testset "evolve()" begin
    @testset "qubit evolution" begin
        σx = Unitaries.σx

        ρ = evolve(σx,States.Qubit([1 0;0 0]))

        @test ρ == [0 0;0 1]
        @test ρ isa States.DensityMatrix
    end

    @testset "ket evolution" begin
        σz = Unitaries.σz

        ψ = evolve(σz, States.Ket([1, 1]/sqrt(2)))

        @test ψ isa States.Ket
        @test ψ == [1,-1]/sqrt(2)
    end
end

end
