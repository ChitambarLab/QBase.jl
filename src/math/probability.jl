export is_probability_distribution, Marginals
export is_conditional_distribution, Conditionals
export joint_probabilities, outcome_probabilities, cvx_combo

"""
    is_probability_distribution( probabilities :: Vector ) :: Bool

Returns `true` if the provided vector is a valid probability distribution:
* `sum(probabilities) ≈ 1`
* `p[i] ≥ 0 ∀ i`
"""
function is_probability_distribution(probabilities::Vector)::Bool
    is_positive = all(p -> (p > 0) || isapprox(p, 0, atol=10e-7), probabilities)
    is_normalized = isapprox(sum(probabilities), 1, atol=10e-6)

    (is_positive & is_normalized)
end

"""
    Marginals( probabilities :: Vector ) <: AbstractVector{Float64}

A struct representing a marginal probability distribution function. All elements
in the marginal distribution must be positive and their sum must be one.
"""
struct Marginals <: AbstractVector{Float64}
    probabilities :: Vector{Float64}
    Base.size(p::Marginals) = size(p.probabilities)
    Base.getindex(p::Marginals, id::Int) = getindex(p.probabilities, id...)
    Base.setindex!(p::Marginals, v, id::Int) = (p.probabilities[id...] = v)
    Marginals(probabilities) = is_probability_distribution(probabilities) ? new(probabilities) : throw(DomainError(probabilities, "is not a valid probability distribution."))
end

"""
    is_conditional_distribution( probabilities :: Matrix ) :: Bool

Returns `true` if the provided matrix is a valid conditional probability distribution.
Each column corresponds to probabilisitic outcomes conditioned for a distinct input.
"""
function is_conditional_distribution(conditionals :: Matrix) :: Bool
    all(id -> is_probability_distribution(conditionals[:,id]), 1:size(conditionals)[2])
end

@doc raw"""
    Conditionals( probabilities :: Matrix ) <: AbstractMatrix{Float64}

A struct representing a conditional probability distribution function. `Conditionals`
are organized in a matrix with rows correpsonding to outputs and columns corresponding
to inputs. For example if there are ``M`` inputs and ``N`` outputs, the corresponding
conditionals matrix takes the form:

```math
\begin{bmatrix}
p(1|1) & \dots & p(1|M) \\
\vdots & \ddots & \vdots \\
p(N|1) & \dots & p(N|M) \\
\end{bmatrix}
```
"""
struct Conditionals <: AbstractMatrix{Float64}
    probabilities :: Matrix{Float64}
    Base.size(C::Conditionals) = size(C.probabilities)
    Base.getindex(C::Conditionals, id::Vararg{Int,2}) = getindex(C.probabilities, id...)
    Base.setindex!(C::Conditionals, v, id::Vararg{Int,2}) = (C.probabilities[id...] = v)
    Conditionals(probabilities) = is_conditional_distribution(probabilities) ? new(probabilities) : throw(DomainError(probabilities, "is not a valid conditional probability distribution"))
end

# """
# cvx_combo(probabilities, elements):
#
#     Constructs the convex combination of elements as weighted by the probabilities.
#
# Inputs:
#     probabilities: Array, contains a normalized probability distribution of
#                    length equal to the number of elements.
#     elements: Array, a set of data structures which can be scaled and summed.
# Output:
#     cvx_combo: A single data structure representing the convex combination.
# """
function cvx_combo(probabilities, elements)
    if !(QMath.is_probability_distribution(probabilities))
        throw(ArgumentError("probabilities are not normalized/positive"))
    elseif length(probabilities) != length(elements)
        throw(ArgumentError("The number of elements does not match the number of weights"))
    end

    sum(probabilities .* elements)
end

# """
# joint_probability(marginal, conditional)
#
#     Returns the probability of prior and conditioned events both occurring .
#
# Input:
#     marginals: Array, contains floats
#     conditionals: Matrix, rows are valid pdf for a given prior, cols correspond
#                   to distinct outputs.
# Output:
#     joint_prob: float, joint probability
# """
function joint_probabilities(marginals::Marginals, conditionals::Conditionals) :: Matrix{Float64}
    joint_probs = marginals' .* conditionals
    joint_probs
end

# """
# outcome_probabilities(priors, conditionals)
#
#     Returns the probability of each outcome given priors and conditional probabilities.
#
# Inputs:
#     priors: Array, float, is a valid pdf
#     conditionals: Matrix, rows are valid pdf for a given prior, cols correspond
#                   to distinct outputs.
# Outputs:
#     outcome_probs: Array, floats, is a valid pdf
# """
function outcome_probabilities(marginals::Marginals, conditionals::Conditionals)
    outcome_probs = conditionals * marginals
    outcome_probs
end
