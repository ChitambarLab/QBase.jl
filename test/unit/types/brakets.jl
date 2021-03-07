using Test, QBase

using LinearAlgebra

@testset "./src/types/brakets.jl" begin

    @test true

    @testset "is_wave_vector()" begin
        @testset "simple valid wavefunction" begin
            @test is_wave_vector([1;0])
            @test is_wave_vector([1])
            @test is_wave_vector([1;1;-1]/sqrt(3))
            @test is_wave_vector(1;im]/sqrt(2))
        end

        @testset "invalid wavefunctions" begin
            @test !is_wave_vector([0;0])
            @test !is_wave_vector([1;1])
            @test !is_wave_vector([1 0;0 0])
        end

        @testset "valid input vector types" begin
            @test is_wave_vector([1,0])
            @test is_wave_vector([1,0]')
            @test is_wave_vector([1;0])
            @test is_wave_vector(ones(Int64, 1,2)/sqrt(2))
            @test is_wave_vector(ones(Int64, 1,2)'/sqrt(2))
        end

        @testset "valid inputs eltypes" begin
            @test is_wave_vector([1,0]::Vector{Int64})
            @test is_wave_vector([1.0,0.0]::Vector{Float64})
            @test is_wave_vector([1//1,0//1]::Vector{Rational{Int64}})
            @test is_wave_vector([1.0+im*0.0,0.0+im*0.0 ]::Vector{Complex{Float64}})
        end
    end
end
