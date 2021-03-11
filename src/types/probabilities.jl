export ProbabilityDistribution, ConditionalDistribution, JointProbabilityDistribution
export is_probability_distribution, is_conditional_distribution
export Probabilities, Conditionals, JointProbabilities
export outcome_probabilities

"""
    abstract type ProbabilityDistribution <: AbstractVector{Real} end

The abstract type representing a discrete probability distribution.
"""
abstract type ProbabilityDistribution <: AbstractVector{Real} end
Base.size(probs::ProbabilityDistribution) = size(probs.distribution)
Base.getindex(probs::ProbabilityDistribution, id::Int) = getindex(probs.distribution, id...)
Base.setindex!(probs::ProbabilityDistribution, v, id::Int) = (probs.distribution[id...] = v)

"""
    abstract type ConditionalDistribution <: AbstractMatrix{Real} end

The abstract type representing a discrete conditional probability distribution.
"""
abstract type ConditionalDistribution <: AbstractMatrix{Real} end
Base.size(probs::ConditionalDistribution) = size(probs.distribution)
Base.getindex(probs::ConditionalDistribution, id::Vararg{Int,2}) = getindex(probs.distribution, id...)
Base.setindex!(probs::ConditionalDistribution, v, id::Vararg{Int,2}) = (probs.distribution[id...] = v)

"""
    abstract type JointProabilityDistribution <: AbstractMatrix{Real} end

The abstract type representing a discrete joint probability distribution.
"""
abstract type JointProbabilityDistribution <: AbstractMatrix{Real} end
Base.size(probs::JointProbabilityDistribution) = size(probs.distribution)
Base.getindex(probs::JointProbabilityDistribution, id::Vararg{Int,2}) = getindex(probs.distribution, id...)
Base.setindex!(probs::JointProbabilityDistribution, v, id::Vararg{Int,2}) = (probs.distribution[id...] = v)

"""
    is_probability_distribution(
        probabilities :: Vector{<:Real};
        atol=ATOL :: Float64
    ) :: Bool

Returns `true` if the provided vector is a valid probability distribution:
* `sum(probabilities) ≈ 1`
* `probabilities[i] ≥ 0 ∀ i`
"""
function is_probability_distribution(probabilities :: AbstractVector{<:Real}; atol=ATOL :: Float64) :: Bool
    _is_prob(prob :: Real) :: Bool = begin
        (prob ≥ 0) || isapprox(prob, 0, atol=atol)
    end

    if findfirst(!_is_prob, probabilities) !== nothing
        return false
    end

    if !isapprox(sum(probabilities), 1, atol=atol)
        return false
    end

    return true
end
is_probability_distribution(::ProbabilityDistribution) = true
is_probability_distribution(::JointProbabilityDistribution) = true


"""
    is_conditional_distribution( probabilities :: AbstractMatrix{<:Real}; atol=ATOL :: Flaot64 ) :: Bool

Returns `true` if the provided matrix is column stochastic. That is, each column
is a valid probability distribution.
"""
function is_conditional_distribution(conditionals :: Matrix{<:Real}; atol=ATOL :: Float64) :: Bool
    all(col -> is_probability_distribution(col, atol=atol), eachcol(conditionals))
end
is_conditional_distribution(::ConditionalDistribution) = true

"""
    Probabilities( distribution :: Vector{<:Real}; atol=ATOL :: Float64 ) <: ProbabilityDistribution

A struct representing a discrete probability distribution. All elements
in the marginal distribution must be positive and their sum must be one.
"""
struct Probabilities{T} <: ProbabilityDistribution
    distribution :: Vector{T}
    atol :: Float64
    Probabilities(
        distribution :: Vector{<:Real};
        atol=ATOL :: Float64
    ) = is_probability_distribution(distribution, atol=atol) ? new{eltype(distribution)}(distribution,atol) : throw(
        DomainError(distribution, "is not a valid probability distribution.")
    )
end

"""
    Conditionals( distribution :: Matrix{<:Real}; atol=ATOL :: Float64 ) <: ConditionalDistribution

A struct representing a conditional probability distribution function. `Conditionals`
are organized in a matrix with rows correpsonding to outputs and columns corresponding
to inputs. For example if there are ``M`` inputs and ``N`` outputs, the corresponding
conditionals matrix takes the form:

```math
\\begin{bmatrix}
p(1|1) & \\dots & p(1|M) \\\\
\\vdots & \\ddots & \\vdots \\\\
p(N|1) & \\dots & p(N|M) \\\\
\\end{bmatrix}
```
"""
struct Conditionals{T} <: ConditionalDistribution
    distribution :: Matrix{T}
    atol :: Float64
    Conditionals(
        distribution :: Matrix{<:Real};
        atol=ATOL :: Float64
    ) = is_conditional_distribution(
            distribution, atol=atol
        ) ? new{eltype(distribution)}(distribution, atol) : throw(
            DomainError(distribution, "is not a valid conditional probability distribution.")
        )
end

"""
    JointProbabilities(
        distribution :: Matrix{<:Real};
        atol=ATOL :: Float64
    ) <: JointProbabilityDistribution

A struct representing a discrete probability distribution. All elements
in the marginal distribution must be positive and their sum must be one.
For convenience, the `JointProbabilities` can also be constructed by passing in
a `ConditionalDistribution` and a `ProbabilityDistribution`,

    JointProbabilities(
        priors :: AbstractVector{<:Real},
        conditionals :: AbstractMatrix{<:Real};
        atol=ATOL :: Float64
    )

Or, two `ProbabilityDistribution`s can be provided

    JointProbabilities(
        priors1 :: AbstractVector{<:Real},
        priors2 :: AbstractVector{<:Real};
        atol=ATOL :: Float64
    )

A `DomainError` is thrown if the joint probability distribution is invalid.
"""
struct JointProbabilities{T} <: JointProbabilityDistribution
    distribution :: Matrix{T}
    atol :: Float64
    JointProbabilities(
        distribution :: Matrix{<:Real};
        atol=ATOL :: Float64
    ) = is_probability_distribution(
            distribution[:], atol=atol
        ) ? new{eltype(distribution)}(distribution, atol) : throw(
            DomainError(distribution, "is not a valid joint probability distribution.")
        )
    JointProbabilities(
        priors :: AbstractVector{<:Real},
        conditionals :: AbstractMatrix{<:Real};
        atol=ATOL :: Float64
    ) = JointProbabilities(priors' .* conditionals, atol=atol)
    JointProbabilities(
        conditionals :: AbstractMatrix{<:Real},
        priors :: AbstractVector{<:Real};
        atol=ATOL :: Float64
    ) = JointProbabilities(priors' .* conditionals, atol=atol)
    JointProbabilities(
        priors1 :: AbstractVector{<:Real},
        priors2 :: AbstractVector{<:Real};
        atol=ATOL :: Float64
    ) = JointProbabilities(priors2 * priors1', atol=atol)
end

"""
    *(
        C :: ConditionalDistribution,
        P :: ProbabilityDistribution;
        atol=ATOL :: Float64
    ) :: ProbabilityDistribution

Multiplication of a conditional with a prior distribution yields the outcome probabilites.
"""
*(
    conditional :: ConditionalDistribution,
    prior :: ProbabilityDistribution;
    atol=ATOL :: Float64
) :: ProbabilityDistribution = Probabilities(conditional.distribution * prior.distribution, atol=atol)

"""
    outcome_probabilities(
        conditionals::ConditionalDistribution,
        priors :: ProbabilityDistribution;
        atol=ATOL :: Float64
    ) :: ProbabilityDistribution

Returns the probability of each outcome given priors and conditional probabilities.
For convenience, this method can be called with an abstract vector/matrix

    outcome_probabilities(
        conditionals :: AbstractMatrix{<:Real},
        priors :: AbstractVector{<:Real};
        atol=ATOL :: Float64
    )

Furthermore, the ordering of the arguments can be reversed.
A `DomainError` is thrown if the constructed vector is not a valid probability
distribution.
"""
function outcome_probabilities(
    conditionals :: ConditionalDistribution,
    priors :: ProbabilityDistribution;
    atol=ATOL :: Float64
) :: ProbabilityDistribution
    *(conditionals,priors,atol=atol)
end
function outcome_probabilities(
    priors :: ProbabilityDistribution,
    conditionals :: ConditionalDistribution;
    atol=ATOL :: Float64
) :: ProbabilityDistribution
    *(conditionals,priors, atol=atol)
end
function outcome_probabilities(
    conditionals :: AbstractMatrix{<:Real},
    priors :: AbstractVector{<:Real};
    atol=ATOL :: Float64
) :: ProbabilityDistribution
    Probabilities(conditionals*priors, atol=atol)
end
function outcome_probabilities(
    priors :: AbstractVector{<:Real},
    conditionals :: AbstractMatrix{<:Real};
    atol=ATOL :: Float64
) :: ProbabilityDistribution
    Probabilities(conditionals*priors, atol=atol)
end
