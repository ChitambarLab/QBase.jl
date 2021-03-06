using Test, LinearAlgebra
using QBase

@testset "./src/types/states.jl" begin

@testset "is_density_matrix()" begin
    @testset "invalid matrices" begin
        @test !is_density_matrix([1 0;0 1])
        @test !is_density_matrix([0.5 1;-1 0.5])
        @test !is_density_matrix([0.5 1;1 0.5])
        @test !is_density_matrix([1.5 0;0 -0.5])
        @test !is_density_matrix([0.5 im;-im 0.5])
        @test !is_density_matrix([1 4;4 1])
        @test !is_density_matrix([0.5 0;0 0.5; 0 0])
    end
    @testset "valid matrices" begin
        @test is_density_matrix([0.5 0;0 0.5]::Matrix{Float64})
        @test is_density_matrix([1 0;0 0]::Matrix{Int64})
        @test is_density_matrix([0 0;0 1])
        @test is_density_matrix([0.5 0.5;0.5 0.5])
        @test is_density_matrix([0.5 0.5im;-0.5im 0.5]::Matrix{Complex{Float64}})
        @test is_density_matrix(State([1 0;0 0]))
    end
    @testset "using atol" begin
        ϵ = 1e-6

        ρ_almost_tr1 = [0.5 0;0 0.5+1e-6]
        @test !is_density_matrix(ρ_almost_tr1)
        @test is_density_matrix(ρ_almost_tr1, atol=1e-5)

        ρ_almost_herm = [1+1e-6 0;0 -1e-6]
        @test !is_density_matrix(ρ_almost_herm)
        @test is_density_matrix(ρ_almost_herm, atol=1e-5)
    end
end

@testset "_density_matrix_error()" begin
    @testset "non-hermitian error" begin
        try
            QBase._density_matrix_error([0.5 0;0.5 0.5])
            @test false
        catch err
            @test err isa DomainError
            @test err.val == [0.5 0;0.5 0.5]
            @test err.msg == "Density matrix `ρ` is not hermitian."
        end
    end

    @testset "non-trace-one error" begin
        try
            QBase._density_matrix_error([0.5 0;0 0.6])
            @test false
        catch err
            @test err isa DomainError
            @test err.val == [0.5 0;0 0.6]
            @test err.msg == "Density matrix `ρ` is not trace-one."
        end
    end

    @testset "non-positive semi-definite error" begin
        try
            QBase._density_matrix_error([1.5 0;0 -0.5])
            @test false
        catch err
            @test err isa DomainError
            @test err.val == [1.5 0;0 -0.5]
            @test err.msg == "Density matrix `ρ` is not positive semi-definite."
        end
    end
end

@testset "State()" begin
    @testset "type inheritance" begin
        ρ_int = State([1 0;0 0])
        @test ρ_int isa State{Int64}
        @test ρ_int isa Operator
        @test ρ_int == [1 0;0 0]

        ρ_float = State([1. 0.;0. 0.])
        @test ρ_float isa State{Float64}
        @test ρ_float isa Operator
        @test ρ_float == [1 0;0 0]

        ρ_complex = State([0.5 0.5im;-0.5im 0.5])
        @test ρ_complex isa State{Complex{Float64}}
        @test ρ_complex isa Operator
        @test ρ_complex == [0.5 0.5im;-0.5im 0.5]

        ρ_rational = State([1//2 1//2; 1//2 1//2])
        @test ρ_rational isa State{Rational{Int64}}
        @test ρ_rational isa Operator
        @test ρ_rational == [0.5 0.5;0.5 0.5]
    end

    @testset "adjoint matrix input" begin
        ρ = State([0.5 -0.5im;0.5im 0.5]')
        @test ρ isa State{Complex{Float64}}
        @test ρ == [0.5 -0.5im;0.5im 0.5]
    end

    @testset "using atol" begin
        ϵ = 1e-5

        @test_throws DomainError State([1 0;0 ϵ])
        @test State([1 0;0 ϵ], atol=1e-4) isa State{Float64}
    end
end

@testset "kron() for states" begin
    ρ0_int = State([1 0;0 0])
    ρ0_float = State([1. 0;0 0])

    ρ1_int = State([0 0;0 1])
    ρ1_float = State([0 0;0 1.])

    ρ01 = kron(ρ0_int,ρ1_int)
    @test ρ01 isa State{Int64}
    @test ρ01 == [0 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0]

    ρ11 = kron(ρ1_int, ρ1_float)
    @test ρ11 isa State{Float64}
    @test ρ11 == [0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 1]

    ρ100 = kron(ρ1_float, ρ0_float, ρ0_int)
    @test ρ100 isa State{Float64}
    @test ρ100 ==  [
        0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0;
    ]

    @testset "using atol" begin
        ϵ = 1e-7

        ρ = State([0.5 0.5;0.5 0.5 + ϵ])

        @test_throws DomainError kron(ρ,ρ)

        ρρ = kron(ρ, ρ, atol=1e-6)
        @test ρρ isa State{Float64}
        @test isapprox(ρρ, ones(Float64,4,4)/4, atol=1e-6)
        @test ρρ.atol == 1e-6
    end
end

@testset "partial_trace() for states" begin
    @testset "maximally entangled bell state" begin
        ρ_AB = State([1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]/2)

        ρ_A = partial_trace(ρ_AB, [2,2], 2)
        @test ρ_A isa State{Float64}
        @test ρ_A == [1 0;0 1]/2

        ρ_B = partial_trace(ρ_AB, [2,2], 1)
        @test ρ_B isa State{Float64}
        @test ρ_B == [1 0;0 1]/2
    end

    @testset "separable qubit states" begin
        ρ_x0 = State([1 1;1 1]/2)
        ρ_0 = State([1 0;0 0])

        ρ_AB = kron(ρ_x0, ρ_0)
        @test ρ_AB isa State{Float64}

        ρ_A = partial_trace(ρ_AB, [2,2], 2)
        @test ρ_A isa State{Float64}
        @test ρ_A == ρ_x0

        ρ_B = partial_trace(ρ_AB, [2,2], 1)
        @test ρ_B isa State{Float64}
        @test ρ_B == ρ_0
    end

    @testset "using atol" begin
        ρ_AB = State([1 1 1 1;1 1 1 1;1 1 1 1;1 1 1 1], atol= 3)

        @test_throws DomainError partial_trace(ρ_AB, [2,2], 2)

        ρ_A = partial_trace(ρ_AB, [2,2], 2, atol=3)
        @test ρ_A isa State{Int64}
        @test ρ_A == [2 2;2 2]
        @test ρ_A.atol == 3
    end
end

@testset "rank() for State types" begin
    @test 1 == rank(pure_state([1,-im]/sqrt(2)))
    @test 2 == rank(State([1 0;0 1]/2))
    @test 1 == rank(State([0.9 0;0 0.1], atol=0.11))
    @test 2 == rank(State([0.9 0;0 0.1]))
end

@testset "eigvals() for State types" begin
    priors = Probabilities([1/3,1/3,1/3])

    qb1 = bloch_qubit_state(2*π/3,0)
    qb2 = bloch_qubit_state(5*π/12,π)
    qb3 = bloch_qubit_state(0,0)
    ρ_states = [qb1,qb2,qb3]

    ρ_mix = mixed_state(priors, ρ_states)
    @test isapprox(eigvals(ρ_mix), [0.372, 0.628], atol=1e-3)
end

@testset "is_pure()" begin
    @test is_pure(State([1 0;0 0]))
    @test is_pure([0 0;0 1])

    @test !is_pure(State([0.5 0;0 0.5]))
    @test_throws DomainError is_pure(State([1 0;0 1]))

    @test !is_pure(State([1-1e-5 0;0 1e-5]))
    @test is_pure(State([1-1e-5 0;0 1e-5], atol=1e-4))

end

@testset "is_mixed()" begin
    @test !is_mixed(State([1 0;0 0]))
    @test !is_mixed([0 0;0 1])

    @test is_mixed(State([0.5 0;0 0.5]))
    @test_throws DomainError is_mixed(State([1 0;0 1]))

    @test is_mixed(State([1-1e-5 0;0 1e-5]))
    @test !is_mixed(State([1-1e-5 0;0 1e-5], atol=1e-4))
end

end
