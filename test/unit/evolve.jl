using Test
using QBase

@testset "./src/evolve.jl" begin

@testset "*(unitary, ket)" begin
    ket = Ket([1.,0], atol=1e-6)
    evo_ket = σx*ket
    @test evo_ket isa Ket{Float64}
    @test evo_ket == [0,1]
    @test evo_ket.atol == 1e-6

    evo_ket2 = σz*σy*Ket([1,0])
    @test evo_ket2 isa Ket{Complex{Int64}}
    @test evo_ket2 == [0,-1im]
    @test evo_ket2.atol == QBase.ATOL
end

@testset "*(bra, unitary)" begin
    bra = Bra([1.,0], atol=1e-6)
    rx = qubit_rotation(π/2)
    evo_bra = bra*rx
    @test evo_bra isa Bra{Complex{Float64}}
    @test evo_bra ≈ [1 -1im]/sqrt(2)
    @test evo_bra.atol == 1e-6

    bra = Bra([1.,0], atol=1e-6)
    evo_bra = bra*rx'
    @test evo_bra isa Bra{Complex{Float64}}
    @test evo_bra ≈[1 1im]/sqrt(2)
end

@testset "*(unitary, state,  unitary)"  begin
    M = σx * State([1. 0;0 0]) * σz
    @test M == [0 0;1 0]
    @test M isa Matrix{Float64}

    M = σy * State([1 0;0 0]) * σy'
    @test M ==  [0 0;0 1]
    @test M isa Matrix{Complex{Int64}}
end

@testset "evolve()" begin
    @testset "qubit evolution" begin
        ρ = evolve(σx,State([1 0;0 0]))

        @test ρ == [0 0;0 1]
        @test ρ isa State{Int64}
    end

    @testset "ket evolution" begin
        ψ = evolve(σz, Ket([1, 1]/sqrt(2)))

        @test ψ isa Ket{Float64}
        @test ψ == [1,-1]/sqrt(2)
    end
end

@testset "evolve(::ChoiOp)" begin
    ρ_mix = [1 0;0 1]/2
    Λ_depol = ChoiOp(x -> ρ_mix, [2,2])
    ρ = [1 0;0 0]

    ρ_out = evolve(Λ_depol, ρ)
    @test ρ_out isa Matrix
    @test ρ_out == ρ_mix

    ρ_out = evolve(Λ_depol, State(ρ))
    @test ρ_out isa State
    @test ρ_out == ρ_mix
end

@testset "evolve(::KrausChannel)" begin
    𝒩_amp_damp = KrausChannel([
        [1 0;0 0.5],
        [0 sqrt(3/4);0 0]
    ])
    ρ = [1 1;1 1] / 2

    ρ_out = evolve(𝒩_amp_damp, ρ)
    @test ρ_out isa Matrix{Float64}
    @test ρ_out == [7/8 1/4;1/4 1/8]

    ρ_out = evolve(𝒩_amp_damp, State(ρ))
    @test ρ_out isa State{Float64}
    @test ρ_out == [7/8 1/4;1/4 1/8]
end

end
