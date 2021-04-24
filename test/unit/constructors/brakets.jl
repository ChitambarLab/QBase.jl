using Test, QBase

@testset "./src/constructors/brakets.jl" begin

@testset "bloch_qubit_ket()" begin
    @testset "complex bloch qubit" begin
        ψ1 = bloch_qubit_ket(π/2,π/2)

        @test ψ1 isa Ket{Complex{Float64}}
        @test ψ1 ≈ [1,im]/sqrt(2)

        ψ2 = bloch_qubit_ket(π/2,3π/2)

        @test ψ2 isa Ket{Complex{Float64}}
        @test ψ2 ≈ [1,-im]/sqrt(2)

        @test_throws DomainError bloch_qubit_ket(-0.1,π/2)
        @test_throws DomainError bloch_qubit_ket(15π/7,π/2)
        @test_throws DomainError bloch_qubit_ket(π/2,-0.1)
        @test_throws DomainError bloch_qubit_ket(π/2,15π/7)
    end

    @testset "real bloch qubit" begin
        ψ1 = bloch_qubit_ket(π/2)

        @test ψ1 isa Ket{Float64}
        @test ψ1 ≈ [1,1]/sqrt(2)

        ψ2 = bloch_qubit_ket(-π/2)

        @test ψ2 isa Ket{Float64}
        @test ψ2 ≈ [1,-1]/sqrt(2)

        @test_throws DomainError bloch_qubit_ket(-2π - 0.1)
        @test_throws DomainError bloch_qubit_ket(2π + 0.1)
    end

    @testset "equivalence between real and complex versions" for θ in 0:0.2:2π
        @test bloch_qubit_ket(θ) ≈ bloch_qubit_ket(θ,0)
        @test isapprox(bloch_qubit_ket(-θ), bloch_qubit_ket(θ,π), atol=1e-7)
    end
end

@testset "computational_basis_kets()" begin
    comp_basis_kets = computational_basis_kets(3)
    @test comp_basis_kets isa Vector{Ket{Int64}}
    @test comp_basis_kets == [[1,0,0],[0,1,0],[0,0,1]]
    @test computational_basis_kets(8) == [
        [1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]
    ]
end

@testset "mirror_symmetric_qubit_kets()" begin
    @test mirror_symmetric_qubit_kets(0) isa Vector{Ket{Float64}}

    @test isapprox(mirror_symmetric_qubit_kets(0), [[1,0],[1,0],[1,0]])
    @test isapprox(mirror_symmetric_qubit_kets(π/4), [[1,0],[1,1]/sqrt(2),[1,-1]/sqrt(2)])
    @test isapprox(mirror_symmetric_qubit_kets(π/3),[[1,0],[0.5,(sqrt(3)/2)],[0.5,-sqrt(3)/2]])
    @test isapprox(mirror_symmetric_qubit_kets(π/2), [[1,0],[0,1],[0,-1]], atol=1e-16)

    @test_throws DomainError mirror_symmetric_qubit_kets(-0.1)
    @test_throws DomainError mirror_symmetric_qubit_kets(3π/5)
end

@testset "planar_symmetric_qubit_kets()" begin
    @testset "orthogonal qubits" begin
        kets = planar_symmetric_qubit_kets(2)
        @test kets isa Vector{Ket{Float64}}
        @test isapprox(kets, [[1,0],[0,1]], atol=1e-7)
    end

    @testset "scanning over cases verifying right number are computed" begin
        for n in 2:100
            kets = planar_symmetric_qubit_kets(n)

            @test length(kets) == n
            @test kets isa Vector{Ket{Float64}}
        end
    end

    @test_throws(DomainError, planar_symmetric_qubit_kets(1))
end

@testset "trine_qubit_kets()" begin
    trine_kets = trine_qubit_kets()
    @test trine_kets isa Vector{Ket{Float64}}
    @test isapprox(trine_kets, mirror_symmetric_qubit_kets(π/3), atol=1e-7)
    @test isapprox(trine_kets, planar_symmetric_qubit_kets(3), atol=1e-7)
end

@testset "sic_qubit_kets()" begin
    sic_kets = sic_qubit_kets()
    @test sic_kets isa Vector{Ket{Complex{Float64}}}

    θ_sic = 2*acos(1/sqrt(3))

    @test sic_kets[1] ≈ bloch_qubit_ket(0, 0)
    @test sic_kets[2] ≈ bloch_qubit_ket(θ_sic, 0)
    @test sic_kets[3] ≈ bloch_qubit_ket(θ_sic, 2π/3)
    @test sic_kets[4] ≈ bloch_qubit_ket(θ_sic, 4π/3)
end

@testset "bb84_qubit_kets" begin
    bb84_kets = bb84_qubit_kets()
    @test bb84_kets isa Vector{Ket{Float64}}
    @test isapprox(bb84_kets, planar_symmetric_qubit_kets(4), atol=1e-7)
end

@testset "bell_kets()" begin
    kets = bell_kets()
    k0 = [1,0]
    k1 = [0,1]

    @test kets isa Vector{Ket{Float64}}
    @test length(kets) == 4
    @test kets[1] == 1/sqrt(2)*(kron(k0,k0)+kron(k1,k1))
    @test kets[2] == 1/sqrt(2)*(kron(k0,k0)-kron(k1,k1))
    @test kets[3] == 1/sqrt(2)*(kron(k0,k1)+kron(k1,k0))
    @test kets[4] == 1/sqrt(2)*(kron(k0,k1)-kron(k1,k0))
end

@testset "generalized_bell_kets()" begin
    @testset "qubit bell states" begin
        kets = generalized_bell_kets(2)

        @test kets isa Vector{Ket{Complex{Float64}}}
        @test isapprox(kets, bell_kets(), atol=1e-7)
    end

    @testset "qutrit bell states" begin
        kets = generalized_bell_kets(3)

        basis = computational_basis_kets(3)
        ψ1 = (kron(basis[1],basis[1]) + kron(basis[2],basis[2]) + kron(basis[3],basis[3]))/sqrt(3)
        ψ2 = (kron(basis[1],basis[1]) + exp(im*2*π/3)*kron(basis[2],basis[2]) + exp(im*4*π/3)*kron(basis[3],basis[3]))/sqrt(3)
        ψ3 = (kron(basis[1],basis[1]) + exp(im*4*π/3)*kron(basis[2],basis[2]) + exp(im*2*π/3)*kron(basis[3],basis[3]))/sqrt(3)

        ψ4 = (kron(basis[1],basis[2]) + kron(basis[2],basis[3]) + kron(basis[3],basis[1]))/sqrt(3)
        ψ5 = (kron(basis[1],basis[2]) + exp(im*2*π/3)*kron(basis[2],basis[3]) + exp(im*4*π/3)*kron(basis[3],basis[1]))/sqrt(3)
        ψ6 = (kron(basis[1],basis[2]) + exp(im*4*π/3)*kron(basis[2],basis[3]) + exp(im*2*π/3)*kron(basis[3],basis[1]))/sqrt(3)

        ψ7 = (kron(basis[1],basis[3]) + kron(basis[2],basis[1]) + kron(basis[3],basis[2]))/sqrt(3)
        ψ8 = (kron(basis[1],basis[3]) + exp(im*2*π/3)*kron(basis[2],basis[1]) + exp(im*4*π/3)*kron(basis[3],basis[2]))/sqrt(3)
        ψ9 = (kron(basis[1],basis[3]) + exp(im*4*π/3)*kron(basis[2],basis[1]) + exp(im*2*π/3)*kron(basis[3],basis[2]))/sqrt(3)

        ψ_match = [ψ1,ψ2,ψ3,ψ4,ψ5,ψ6,ψ7,ψ8,ψ9]

        @test kets isa Vector{Ket{Complex{Float64}}}
        @test kets ≈ ψ_match
    end
end

end
