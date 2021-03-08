using Test, QBase

using LinearAlgebra

@testset "./src/types/brakets.jl" begin

@testset "is_wave_vector()" begin
    @testset "simple valid wavefunction" begin
        @test is_wave_vector([1;0])
        @test is_wave_vector([1])
        @test is_wave_vector([1;1;-1]/sqrt(3))
        @test is_wave_vector([1;im]/sqrt(2))
    end

    @testset "invalid wavefunctions" begin
        @test !is_wave_vector([0;0])
        @test !is_wave_vector([1;1])
        @test !is_wave_vector([1 0;0 0])
        @test !is_wave_vector([1 0;0 0]')
    end

    @testset "valid input vector types" begin
        @test is_wave_vector([1,0])
        @test is_wave_vector([1,0]')
        @test is_wave_vector([1;0])
        @test is_wave_vector(ones(Int64, 1,2)/sqrt(2))
        @test is_wave_vector(ones(Int64, 1,2)'/sqrt(2))
        @test is_wave_vector(Ket([1,0]))
        @test is_wave_vector(Bra([1,0]))
    end

    @testset "valid inputs eltypes" begin
        @test is_wave_vector([1,0]::Vector{Int64})
        @test is_wave_vector([1.0,0.0]::Vector{Float64})
        @test is_wave_vector([1//1,0//1]::Vector{Rational{Int64}})
        @test is_wave_vector([1.0+im*0.0,0.0+im*0.0 ]::Vector{Complex{Float64}})
    end

    @testset "tuning atol" begin
        ϵ = 1e-7

        @test is_wave_vector([1,1], atol=1)
        @test is_wave_vector([1,1], atol=sqrt(2)-1+ϵ)
        @test !is_wave_vector([1,1], atol=sqrt(2)-1-ϵ)
        @test !is_wave_vector([1,1], atol=0)
    end
end

@testset "_wave_vector_error()" begin
    @test_throws DomainError QBase._wave_vector_error([1,0])

    try
        QBase._wave_vector_error([1,0])
        @test false
    catch err
        @test err isa DomainError
        @test err.val == [1, 0]
        @test err.msg ==  "ψ must be a vector that satisfies `norm(ψ,2) != 1`."
    end
end

@testset "Ket()" begin
    @testset "valid input types" begin
        ψ = Ket([1,0])
        @test ψ isa Ket{Int64}
        @test ψ isa AbstractKet{Int64}
        @test ψ == [1,0]

        ψ = Ket([1.,0.])
        @test ψ isa Ket{Float64}
        @test ψ isa AbstractKet{Float64}
        @test ψ == [1,0]

        ψ = Ket([1//1,0//1])
        @test ψ isa Ket{Rational{Int64}}
        @test ψ isa AbstractKet{Rational{Int64}}
        @test ψ == [1,0]

        ψ = Ket([im,0])
        @test ψ isa Ket{Complex{Int64}}
        @test ψ isa AbstractKet{Complex{Int64}}
        @test ψ == [im,0]

        ψ = Ket([im,0.])
        @test ψ isa Ket{Complex{Float64}}
        @test ψ isa AbstractKet{Complex{Float64}}
        @test ψ == [im,0]
    end

    @testset "simple test cases" begin
        valid_cases = [
            [1,0,0],
            [1,0,0,0],
            [1,1,im]/sqrt(3),
            [0.5im,0//1,-sqrt(3)/2]
        ]

        for case in valid_cases
            @test Ket(case) == case
        end
    end

    @testset "input vector formats" begin
        ψ_match = [1,0]

        ψ = Ket([1;0] :: Vector{Int64})
        @test ψ == ψ_match

        ψ = Ket([1,0]' :: Adjoint{Int64,Vector{Int64}})
        @test ψ == ψ_match

        ψ = Ket([1,0] :: Vector{Int64})
        @test ψ == ψ_match

        ψ = Ket(Bra([im,0]))
        @test ψ == [-im,0]

        ψ = Ket([1 0] :: Matrix{Int64})
        @test ψ == ψ_match

        ψ = Ket([1 0]' :: Adjoint{Int64,Matrix{Int64}})
        @test ψ == ψ_match
    end

    @testset "invalid inputs" begin
        @test_throws DomainError Ket([1,1,1])
        @test_throws MethodError Ket(["a","b","c"])
    end

    @testset "tuning atol" begin
        ϵ = 1e-7
        ψ = Ket([1,1], atol=sqrt(2)-1+ϵ)

        @test ψ.atol ≈ sqrt(2)-1+ϵ
        @test ψ isa Ket{Int64}
        @test_throws DomainError Ket([1,1], atol=sqrt(2)-1-ϵ)
    end
end

@testset "Bra()" begin
    @testset "valid input types" begin
        ψ = Bra([1,0])
        @test ψ isa Bra{Int64}
        @test ψ isa AbstractBra{Int64}
        @test ψ == [1 0]

        ψ = Bra([1.,0.])
        @test ψ isa Bra{Float64}
        @test ψ isa AbstractBra{Float64}
        @test ψ == [1 0]

        ψ = Bra([1//1,0//1])
        @test ψ isa Bra{Rational{Int64}}
        @test ψ isa AbstractBra{Rational{Int64}}
        @test ψ == [1 0]

        ψ = Bra([im,0])
        @test ψ isa Bra{Complex{Int64}}
        @test ψ isa AbstractBra{Complex{Int64}}
        @test ψ == [im 0]

        ψ = Bra([im,0.])
        @test ψ isa Bra{Complex{Float64}}
        @test ψ isa AbstractBra{Complex{Float64}}
        @test ψ == [im 0]
    end

    @testset "simple test cases" begin
        valid_cases = [
            [1 0 0],
            [1 0 0 0],
            [1 1 im]/sqrt(3),
            [0.5im 0//1 -sqrt(3)/2]
        ]

        for case in valid_cases
            @test Bra(case) == case
        end
    end

    @testset "input vector formats" begin
        ψ_match = [1 0]

        ψ = Bra([1;0] :: Vector{Int64})
        @test ψ == ψ_match

        ψ = Bra([1,0]' :: Adjoint{Int64,Vector{Int64}})
        @test ψ == ψ_match

        ψ = Bra([1,0] :: Vector{Int64})
        @test ψ == ψ_match

        ψ = Bra(Ket([im,0]))
        @test ψ == [-im 0]

        ψ = Bra([1 0] :: Matrix{Int64})
        @test ψ == ψ_match

        ψ = Bra((ones(Int64, 2,1)/sqrt(2))' :: Adjoint{Float64,Matrix{Float64}})
        @test ψ == ones(Int64, 1,2)/sqrt(2)
    end

    @testset "invalid inputs" begin
        @test_throws DomainError Bra([1,1,1])
        @test_throws MethodError Bra(["a","b","c"])
    end

    @testset "tuning atol" begin
        ϵ = 1e-7
        ψ = Bra([1,1], atol=sqrt(2)-1+ϵ)

        @test ψ.atol ≈ sqrt(2)-1+ϵ
        @test ψ isa Bra{Int64}
        @test_throws DomainError Bra([1,1], atol=sqrt(2)-1-ϵ)
    end
end

@testset "bra and ket multiplication" begin
    @testset "product of Complex{Int64} types" begin
        ket = Ket([im,0])
        bra = Bra([-im,0])

        inner_product = bra*ket
        @test inner_product isa Complex{Int64}
        @test inner_product == 1

        outer_product = ket*bra
        @test outer_product isa Matrix{Complex{Int64}}
        @test outer_product == [1 0;0 0]
    end

    @testset "product of Complex{Int64} types" begin
        ket = Ket([1,0])
        bra = Bra([1,-im]/sqrt(2))

        inner_product = bra*ket
        @test inner_product isa Complex{Float64}
        @test inner_product ≈ 1/sqrt(2)

        outer_product = ket*bra
        @test outer_product isa Matrix{Complex{Float64}}
        @test outer_product ≈ [1 -im;0 0]/sqrt(2)
    end
end

@testset "bra and ket conversion via adjoint" begin
    ket = Ket([im,0])
    bra = Bra([-im,0])

    ψ = ket'
    @test ψ isa Bra{Complex{Int64}}
    @test ψ == [-im 0]

    ψ = adjoint(ket)
    @test ψ isa Bra{Complex{Int64}}
    @test ψ == [-im 0]

    ψ = bra'
    @test ψ isa Ket{Complex{Int64}}
    @test ψ == [im,0]

    ψ = adjoint(bra)
    @test ψ isa Ket{Complex{Int64}}
    @test ψ == [im,0]
end

@testset "kron() for bras and kets" begin
    @testset "ket" begin
        ket0 = Ket([1,0])
        ket1 = Ket([0,1])

        ket01 = kron(ket0,ket1)
        @test ket01 isa Ket{Int64}
        @test ket01 == [0,1,0,0]

        ket010 = kron(ket0,ket1,ket0)
        @test ket010 isa Ket{Int64}
        @test ket010 == [0,0,1,0,0,0,0,0]
    end

    @testset "ket using atol" begin
        ϵ = 1e-7
        ket0 = Ket([1,1], atol=sqrt(2)-1+ϵ)
        ket1 = Ket([1,-1], atol=sqrt(2)-1+ϵ)

        ket01 = kron(ket0,ket1, atol=2-1+ϵ)
        @test ket01 isa Ket{Int64}
        @test ket01 == [1,-1,1,-1]

        @test_throws DomainError kron(ket0,ket1, atol=2-1-ϵ)
        @test_throws DomainError kron(ket0,ket1)
    end

    @testset "bra" begin
        bra0 = Bra([1,0])
        bra1 = Bra([0,1])

        bra01 = kron(bra0,bra1)
        @test bra01 isa Bra{Int64}
        @test bra01 == [0 1 0 0]

        bra010 = kron(bra0,bra1,bra0)
        @test bra010 isa Bra{Int64}
        @test bra010 == [0 0 1 0 0 0 0 0]
    end

    @testset "bra using atol" begin
        ϵ = 1e-7
        bra0 = Bra([1,1], atol=sqrt(2)-1+ϵ)
        bra1 = Bra([1,-1], atol=sqrt(2)-1+ϵ)

        bra01 = kron(bra0,bra1, atol=2-1+ϵ)
        @test bra01 isa Bra{Int64}
        @test bra01 == [1 -1 1 -1]

        @test_throws DomainError kron(bra0,bra1, atol=2-1-ϵ)
        @test_throws DomainError kron(bra0,bra1)
    end
end

end
