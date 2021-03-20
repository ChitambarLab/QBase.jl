using Test, LinearAlgebra
using QBase

@testset "./src/constructors/unitaries.jl" begin

@testset "pauli constants" begin
    @test σI == [1 0;0 1]
    @test σI isa Unitary{Int64}

    @test σx == [0 1; 1 0]
    @test σx isa Unitary{Int64}

    @test σy == [0 -im; im 0]
    @test σy isa Unitary{Complex{Int64}}

    @test σz == [1 0; 0 -1]
    @test σz isa Unitary{Int64}
end

@testset "qubit_rotation()" begin
    ρx = 0.5*[1 1; 1 1]
    ρz = [1 0; 0 0]

    @testset "rotation about x-axis" begin
        for θ in -π:0.1π:π
            Rx = qubit_rotation(θ, axis="x")
            @test isapprox(Rx*ρz*Rx',  [cos(θ/2)^2 im*sin(θ/2)cos(θ/2); -im*sin(θ/2)cos(θ/2) sin(θ/2)^2], atol=1e-7)
        end
    end

    @testset "rotation about y-axis" begin
        for θ in -π:0.1π:π
            Ry = qubit_rotation(θ, axis="y")
            @test isapprox(Ry*ρz*Ry',  [cos(θ/2)^2 sin(θ/2)cos(θ/2); sin(θ/2)cos(θ/2) sin(θ/2)^2], atol=1e-7)
        end
    end

    @testset "rotation about z-axis" begin
        for θ in -π:0.1π:π
            Rz = qubit_rotation(θ, axis="z")
            @test isapprox(Rz*ρz*Rz', ρz)
            @test isapprox(Rz*ρx*Rz', [0.5 0.5*(cos(θ)-im*sin(θ)); 0.5*(cos(θ)+im*sin(θ)) 0.5], atol=1e-7)
        end
    end
end

@testset "random_unitary()" begin
    u2 = random_unitary(2)
    @test u2 isa Unitary{Complex{Float64}}
    @test size(u2) == (2,2)

    u3 = random_unitary(3)
    @test u3 isa Unitary{Complex{Float64}}
    @test size(u3) == (3,3)

    u8 = random_unitary(8)
    @test u8 isa Unitary{Complex{Float64}}
    @test size(u8) == (8,8)

    u3_samples = random_unitary.(fill(3,10))
    @test u3_samples isa Vector{Unitary{Complex{Float64}}}
    @test length(u3_samples) == length(unique(u3_samples))
end

end
