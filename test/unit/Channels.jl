using Test

@testset "./src/Channels.jl" begin

using QBase

@testset "depolarizing()" begin
    @testset "simple Qubit examples" begin
        ρ0 = States.Qubit([1 0;0 0])

        ρ0_out1 = Channels.depolarizing(ρ0, 1)

        @test ρ0_out1 isa States.DensityMatrix
        @test ρ0_out1 == ρ0

        ρ0_out2 = Channels.depolarizing(ρ0, 0)

        @test ρ0_out2 isa States.DensityMatrix
        @test ρ0_out2 == [1/2 0;0 1/2]

        ρ0_out3 = Channels.depolarizing(ρ0, 0.5)

        @test ρ0_out3 isa States.DensityMatrix
        @test ρ0_out3 == [3/4 0;0 1/4]
    end

    @testset "simple DensityMatrix examples" begin
        ρ = States.DensityMatrix([0.5 0 0.5;0 0 0;0.5 0 0.5])

        ρ_out1 = Channels.depolarizing(ρ, 1)

        @test ρ_out1 isa States.DensityMatrix
        @test ρ_out1 == ρ

        ρ_out2 = Channels.depolarizing(ρ, 0.5)

        @test ρ_out2 isa States.DensityMatrix
        @test ρ_out2 == [(1/4+1/6) 0 1/4;0 1/6 0;1/4 0 (1/4+1/6)]
    end

    @testset "errors" begin
        ρ = States.DensityMatrix([1 0;0 0])

        @test_throws DomainError Channels.depolarizing(ρ, 1.1)
        @test_throws DomainError Channels.depolarizing(ρ, -0.1)
    end
end

@testset "erasure()" begin
    @testset "simple Qubit examples" begin
        ρ = States.Qubit([1 0;0 0])

        ρ_out1 = Channels.erasure(ρ, 1)

        @test ρ_out1 isa States.DensityMatrix
        @test ρ_out1 == [1 0 0;0 0 0;0 0 0]

        ρ_out2 = Channels.erasure(ρ, 0)

        @test ρ_out2 isa States.DensityMatrix
        @test ρ_out2 == [0 0 0;0 0 0;0 0 1]

        ρ_out3 = Channels.erasure(ρ, 0.5)

        @test ρ_out3 isa States.DensityMatrix
        @test ρ_out3 == [0.5 0 0;0 0 0;0 0 0.5]
    end

    @testset "simple DensityMatrix examples" begin
        ρ = States.DensityMatrix([0.5 0 0.5;0 0 0;0.5 0 0.5])

        ρ_out1 = Channels.erasure(ρ, 1)

        @test ρ_out1 isa States.DensityMatrix
        @test ρ_out1 == [0.5 0 0.5 0;0 0 0 0;0.5 0 0.5 0;0 0 0 0]

        ρ_out2 = Channels.erasure(ρ, 0.5)

        @test ρ_out2 isa States.DensityMatrix
        @test ρ_out2 == [1/4 0 1/4 0;0 0 0 0;1/4 0 1/4 0;0 0 0 1/2]
    end

    @testset "errors" begin
        ρ = States.DensityMatrix([1 0;0 0])

        @test_throws DomainError Channels.erasure(ρ, 1.1)
        @test_throws DomainError Channels.erasure(ρ, -0.1)
    end
end

end
