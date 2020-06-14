using Test

@testset "QMath/probability.jl" begin

using QBase: QMath

@testset "QMath.is_probability_distribution()" begin
    @testset "valid cases" begin
        @test QMath.is_probability_distribution([1])
        @test QMath.is_probability_distribution([0,0,0,1,0])
        @test QMath.is_probability_distribution(1/100*fill(1,100))
        @test QMath.is_probability_distribution([0.2,0,0,0.1,0.1,0.1,0.3,0.2,0])
    end

    @testset "invalid cases" begin
        @test !QMath.is_probability_distribution([0])
        @test !QMath.is_probability_distribution([0.2,0.1,0.11,0.4])
        @test !QMath.is_probability_distribution([-0.1,1,0.1])
        @test !QMath.is_probability_distribution([-0.1,-0.2])
    end
end

@testset "QMath.Marginals()" begin
    @testset "valid cases" begin
        pdf = QMath.Marginals([0.5,0.5])
        @test pdf isa QMath.Marginals
        @test pdf == [0.5,0.5]

        @test QMath.Marginals([1,0]) == [1.0,0]
    end

    @testset "invalid cases" begin
        @test_throws DomainError QMath.Marginals([0.5,0.4,0.2])
        @test_throws DomainError QMath.Marginals([0.5,-1,1.5])
    end
end

@testset "is_conditional_distribution()" begin
    @testset "valid cases" begin
        @test QMath.is_conditional_distribution([1 0 0;0 1 0;0 0 1])
        @test QMath.is_conditional_distribution([1/3 1/2 1/8;1/3 1/4 0;1/3 1/4 7/8])
    end

    @testset "invalid cases" begin
        @test !QMath.is_conditional_distribution([1 1 1;1 0 0;1 0 0])
        @test !QMath.is_conditional_distribution([1 0 0;0 1 0;0 0 0])
    end
end

@testset "Conditionals()" begin
    @testset "valid cases" begin
        conditionals = QMath.Conditionals([1/3 1/2 1/8;1/3 1/4 0;1/3 1/4 7/8])
        @test conditionals == [1/3 1/2 1/8;1/3 1/4 0;1/3 1/4 7/8]
        @test conditionals isa QMath.Conditionals
    end

    @testset "invalid cases" begin
        @test_throws DomainError QMath.Conditionals([0.5 -0.5;0.5 1.5])
    end
end

@testset "cvx_combo()" begin
    @test_throws ArgumentError QMath.cvx_combo([0.1, 0.8, 0.2], [[1 0 0],[0 1 0],[0 0 1]])
    @test_throws ArgumentError QMath.cvx_combo([0.8, 0.2], [[1 0 0],[0 1 0],[0 0 1]])

    @test QMath.cvx_combo([0.5, 0.5], [[1 1 1;0 0 0;0 0 0],[0 0 0;1 1 1;0 0 0]]) == [0.5 0.5 0.5;0.5 0.5 0.5;0 0 0]
    @test QMath.cvx_combo([0.2, 0.3, 0.5], [[1 0 0 0 1 0 0 0 1],[1 0 0 0 0 1 0 1 0],[0 0 1 1 0 0 0 1 0]]) == [0.5 0 0.5 0.5 0.2 0.3 0 0.8 0.2]
end

@testset "joint_probabilities()" begin
    priors = QMath.Marginals([0.3,0.6,0.1])
    conditionals = QMath.Conditionals([0.3 0.5 0.4;0.7 0.5 0.6])
    @test QMath.joint_probabilities(priors,conditionals) ≈ [0.09 0.30 0.04;0.21 0.30 0.06]
    @test QMath.is_probability_distribution(QMath.joint_probabilities(priors,conditionals)[:])
end

@testset "outcome_probabilities()" begin
    marginals = QMath.Marginals([0.3,0.6,0.1])
    conditionals = QMath.Conditionals([0.3 0.5 0.4;0.7 0.5 0.6])
    @test QMath.outcome_probabilities(marginals,conditionals) ≈ [0.43,0.57]
    @test QMath.is_probability_distribution(QMath.outcome_probabilities(marginals,conditionals))
end

end
