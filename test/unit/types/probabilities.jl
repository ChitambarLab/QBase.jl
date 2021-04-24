using Test, QBase

@testset "./src/types/probabilities.jl" begin

@testset "is_probability_distribution()" begin
    @testset "valid cases" begin
        @test is_probability_distribution([1])
        @test is_probability_distribution([0,0,0,1,0])
        @test is_probability_distribution(1/100*fill(1,100))
        @test is_probability_distribution([0.2,0,0,0.1,0.1,0.1,0.3,0.2,0])
    end

    @testset "invalid cases" begin
        @test !is_probability_distribution([0])
        @test !is_probability_distribution([0.2,0.1,0.11,0.4])
        @test !is_probability_distribution([-0.1,1,0.1])
        @test !is_probability_distribution([-0.1,-0.2])
    end

    @testset "atol" begin
        @test !is_probability_distribution([0.3,0.3,0.39])
        @test is_probability_distribution([0.3,0.3,0.39], atol=2e-2)
        @test !is_probability_distribution([0.3,0.3,0.41,-0.01])
        @test is_probability_distribution([0.3,0.3,0.41,-0.01], atol=2e-2)
    end

    @testset "distribution types" begin
        @test is_probability_distribution(Probabilities([1,0,0]))
        @test is_probability_distribution(JointProbabilities([1 0 0;0 0 0;0 0 0]))
    end
end

@testset "is_conditional_distribution()" begin
    @testset "valid cases" begin
        @test is_conditional_distribution([1 0 0;0 1 0;0 0 1])
        @test is_conditional_distribution([1/3 1/2 1/8;1/3 1/4 0;1/3 1/4 7/8])
    end

    @testset "invalid cases" begin
        @test !is_conditional_distribution([1 1 1;1 0 0;1 0 0])
        @test !is_conditional_distribution([1 0 0;0 1 0;0 0 0])
    end

    @testset "atol" begin
        @test !is_conditional_distribution([1 0 0;0 1 0;0 0 1.1])
        @test is_conditional_distribution([1 0 0;0 1 0;0 0 1.1], atol=0.11)
    end

    @test is_conditional_distribution(Conditionals([1 0 0;0 1 0;0 0 1]))
end

@testset "Probabilities()" begin
    @testset "valid cases" begin
        probs = Probabilities([0.5,0.5])
        @test probs isa Probabilities{Float64}
        @test probs isa ProbabilityDistribution
        @test probs == [0.5,0.5]

        probs_int = Probabilities([1,0])
        @test probs_int isa Probabilities{Int64}
        @test probs_int isa ProbabilityDistribution
        @test probs_int == [1,0]

        probs_atol = Probabilities([1,0.1],atol=.11)
        @test probs_atol isa Probabilities{Float64}
        @test probs_atol isa ProbabilityDistribution
        @test probs_atol == [1,0.1]
    end

    @testset "invalid cases" begin
        @test_throws DomainError Probabilities([0.5,0.4,0.2])
        @test_throws DomainError Probabilities([0.5,-1,1.5])
    end
end


@testset "Conditionals()" begin
    @testset "valid cases" begin
        conditionals = Conditionals([1/3 1/2 1/8;1/3 1/4 0;1/3 1/4 7/8])
        @test conditionals == [1/3 1/2 1/8;1/3 1/4 0;1/3 1/4 7/8]
        @test conditionals isa Conditionals{Float64}
        @test conditionals isa ConditionalDistribution
    end

    @testset "inexact edge case" begin
        @test_throws DomainError Conditionals([
            0.666667 0.166667 0.166667; 0.166667 0.666667 0.166667; 0.166667 0.166667 0.666667
        ])

        conditionals = Conditionals([
            0.666667 0.166667 0.166667; 0.166667 0.666667 0.166667; 0.166667 0.166667 0.666667
        ], atol=1e-5)
        @test conditionals == [0.666667 0.166667 0.166667; 0.166667 0.666667 0.166667; 0.166667 0.166667 0.666667]
    end

    @testset "invalid cases" begin
        @test_throws DomainError Conditionals([0.5 -0.5;0.5 1.5])
    end
end

@testset "JointProbabilities()" begin
    @testset "matrix constructor" begin
        joints = JointProbabilities([1 0 0;0 0 0;0 0 0])
        @test joints isa JointProbabilityDistribution
        @test joints isa JointProbabilities{Int64}
        @test joints == [1 0 0;0 0 0;0 0 0]
    end

    @testset "errors" begin
        @test_throws DomainError JointProbabilities([1 1 0;0 0 0;0 0 0])
    end

    @testset "atol" begin
        @test_throws DomainError JointProbabilities([1.1 0 0;0 0 0;0 0 -0.1])
        @test JointProbabilities([1.1 0 0;0 0 0;0 0 -0.1], atol=0.11) isa JointProbabilityDistribution

        @test_throws DomainError JointProbabilities([1,0,0],[1,0,0.1])
        @test JointProbabilities([1,0,0],[1,0,0.1], atol=0.11) isa JointProbabilityDistribution
    end

    @testset "priors and conditionals" begin
        priors_vec = [0.3,0.6,0.1]
        priors = Probabilities(priors_vec)

        conditionals_mat = [0.3 0.5 0.4;0.7 0.5 0.6]
        conditionals = Conditionals(conditionals_mat)

        joint_probs = JointProbabilities(priors,conditionals)
        @test joint_probs ≈ [0.09 0.30 0.04;0.21 0.30 0.06]
        @test joint_probs isa JointProbabilityDistribution
        @test joint_probs isa JointProbabilities{Float64}

        joint_probs = JointProbabilities(priors_vec,conditionals_mat)
        @test joint_probs ≈ [0.09 0.30 0.04;0.21 0.30 0.06]
        @test joint_probs isa JointProbabilityDistribution
        @test joint_probs isa JointProbabilities{Float64}

        joint_probs = JointProbabilities(conditionals_mat, priors_vec)
        @test joint_probs ≈ [0.09 0.30 0.04;0.21 0.30 0.06]
        @test joint_probs isa JointProbabilityDistribution
        @test joint_probs isa JointProbabilities{Float64}

        joint_probs = JointProbabilities(priors_vec,conditionals)
        @test joint_probs ≈ [0.09 0.30 0.04;0.21 0.30 0.06]
        @test joint_probs isa JointProbabilityDistribution
        @test joint_probs isa JointProbabilities{Float64}
    end

    @testset "priors and priors" begin
        priors1_vec = [0.9,0,0.1]
        priors1 = Probabilities(priors1_vec)

        priors2_vec = [0,0,0.5,0.5]
        priors2 = Probabilities(priors2_vec)

        joint_probs = JointProbabilities(priors1, priors2)
        @test joint_probs == [0 0 0;0 0 0;0.45 0 0.05;0.45 0 0.05]
        @test joint_probs isa JointProbabilities{Float64}

        joint_probs = JointProbabilities(priors2, priors1)
        @test joint_probs == [0 0 0.45 0.45;0 0 0 0;0 0 0.05 0.05]

        joint_probs = JointProbabilities(priors1_vec, priors2_vec)
        @test joint_probs == [0 0 0;0 0 0;0.45 0 0.05;0.45 0 0.05]
    end
end

@testset "outcome_probabilities()" begin
    @testset "multiplication" begin
        priors = Probabilities([0.3,0.3,0.4])
        conditionals = Conditionals([1 0 0;0 1 0;0 0 1])

        probs = conditionals*priors
        @test probs isa Probabilities{Float64}
        @test probs isa ProbabilityDistribution
        @test probs == [0.3,0.3,0.4]
    end

    @testset "ordering of inputs" begin
        priors_vec = [0.3,0.6,0.1]
        priors = Probabilities(priors_vec)

        conditionals_mat = [0.3 0.5 0.4;0.7 0.5 0.6]
        conditionals = Conditionals(conditionals_mat)

        out_probs = outcome_probabilities(priors,conditionals)
        @test out_probs ≈ [0.43,0.57]
        @test out_probs isa Probabilities{Float64}

        out_probs = outcome_probabilities(conditionals,priors)
        @test out_probs ≈ [0.43,0.57]
        @test out_probs isa Probabilities{Float64}

        out_probs = outcome_probabilities(conditionals_mat,priors_vec)
        @test out_probs ≈ [0.43,0.57]
        @test out_probs isa Probabilities{Float64}

        out_probs = outcome_probabilities(priors_vec,conditionals_mat)
        @test out_probs ≈ [0.43,0.57]
        @test out_probs isa Probabilities{Float64}
    end
end

end
