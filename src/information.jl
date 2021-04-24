export shannon_entropy, von_neumann_entropy, joint_entropy, conditional_entropy
export holevo_bound, mutual_information
export success_probability, error_probability

"""
    shannon_entropy( probabilities :: AbstractVector ) :: Float64

The classical entropy of a probability distribution:

```math
S = -\\sum_{i=1}^n p_i \\log_2(p_i)
```

A `DomainError` is thrown if input `probabilities` does not satisfy [`is_probability_distribution`](@ref).
"""
function shannon_entropy( probabilities :: AbstractVector ) :: Float64
    is_probability_distribution(probabilities) || Probabilities(probabilities)

    entropy = -1*sum(map(
        (p) -> isapprox(p, 0, atol=1e-7) ? 0 : p*log2(p),
        probabilities
    ))

    entropy
end

"""
    von_neumann_entropy( ρ :: AbstractMatrix ) :: Float64

The von neumann entropy of a density matrix:


```math
S(\\rho) = - \\sum_j  \\lambda_j \\log_2(\\lambda_j)
```

where ``\\lambda_j`` are the eigenvalues of quantum  state ``\\rho``.

A `DomainError` is thrown if `ρ` does not satisfy [`is_density_matrix`](@ref).
"""
function von_neumann_entropy(ρ :: AbstractMatrix) :: Float64
    is_density_matrix(ρ) || State(ρ)

    λs = Probabilities(eigvals(ρ))
    entropy = shannon_entropy(λs)

    entropy
end

"""
    holevo_bound(
        priors :: AbstractVector,
        ρ_states :: Vector{<:AbstractMatrix}
    ) :: Float64

The Holevo theorem places a bound on the [`mutual_information`](@ref) ``I(X : Y) \\leq \\mathcal{X}``.
It places a limit on the amount information that can be decoded from a set of quantum states.
For a mixed state ``\\rho = \\sum_i p_i \\rho_i`` the Holevo bound ``\\mathcal{X}`` is

```math
\\mathcal{X}  := S(\\rho) - \\sum_i p_i S(\\rho_i)
```

where ``S(\\rho)`` is the [`von_neumann_entropy`](@ref).


A `DomainError` is thrown  if:
* `priors` does not satisfy [`is_probability_distribution`](@ref).
* Any `ρ ∈ ρ_states` does not satisfy [`is_density_matrix`](@ref).
"""
function holevo_bound(priors::AbstractVector, ρ_states::Vector{<:AbstractMatrix}) :: Float64
    is_probability_distribution(priors) || Probabilities(priors)
    all(is_density_matrix.(ρ_states)) || State.(ρ_states)

    ρ = mixed_state(priors, ρ_states)
    ρ_ent = von_neumann_entropy(ρ)

    ρ_states_ent = von_neumann_entropy.(ρ_states)
    ρ_sum_ent = sum( priors .* ρ_states_ent )

    bound = ρ_ent - ρ_sum_ent

    bound
end

"""
    joint_entropy(priors :: AbstractVector, conditionals :: AbstractMatrix) :: Float64

Returns the entropy for the union of pdf ``P(x,y)``. The joint entropy is the [`shannon_entropy`](@ref)
of the joint probability distribution:

```math
S(X,Y) = - \\sum_{x,y} p(x,y) \\log_2(p(x,y))
```
"""
function joint_entropy(priors::AbstractVector, conditionals::AbstractMatrix) :: Float64
    is_probability_distribution(priors) || Probabilities(priors)
    is_conditional_distribution(conditionals) || Conditionals(conditionals)

    joint_probabilities = JointProbabilities(priors, conditionals)
    shannon_entropy(joint_probabilities[:])
end

"""
    conditional_entropy(priors::AbstractVector, conditionals::AbstractMatrix) :: Float64

Returns the conditional entropy for the system with specified `priors` ``p(x)`` and `conditionals` ``p(y|x)``:

```math
S(Y|X) = - \\sum_{x,y} p(x,y)\\log_2\\left(\\frac{p(y|x)}{p(x)}\\right)
```
"""
function conditional_entropy(priors::AbstractVector, conditionals::AbstractMatrix) :: Float64
    is_probability_distribution(priors) || Probabilities(priors)
    is_conditional_distribution(conditionals) || Conditionals(conditionals)

    joint_probabilities = JointProbabilities(priors, conditionals)
    (num_row, num_col) = size(joint_probabilities)

    conditional_entropy = shannon_entropy(Probabilities(joint_probabilities[:])) + sum(sum(map( row -> map( col ->
        joint_probabilities[row,col]*( isapprox(priors[col], 0, atol=1e-7) ? 0 : log2(priors[col]) ),
    1:num_col), 1:num_row)))

    conditional_entropy
end

"""
    mutual_information(
        priors :: AbstractVector,
        conditionals :: AbstractMatrix
    ) :: Float64

The entropy of the overlap between ``p(x)`` and ``p(y)``. The mutual information
is directly computed from the [`shannon_entropy`](@ref) and [`joint_entropy`](@ref):


```math
I(X : Y) = S(X) + S(Y) - S(X,Y)
```
"""
function mutual_information(priors::AbstractVector, conditionals::AbstractMatrix) :: Float64
    is_probability_distribution(priors) || Probabilities(priors)
    is_conditional_distribution(conditionals) || Conditionals(conditional)

    prior_ent = shannon_entropy(priors)

    outcome_probs = outcome_probabilities(priors,conditionals)
    outcome_ent = shannon_entropy( outcome_probs )
    joint_ent = joint_entropy(priors, conditionals)

    mutual_information = prior_ent + outcome_ent - joint_ent

    mutual_information
end

"""
    mutual_information(
        priors :: AbstractVector,
        ρ_states :: Vector{<:AbstractMatrix},
        Π :: AbstractVector
    ) :: Float64

Computes the classical mutual information for a quantum state and measurement
encoding and decoding.
The conditional probabilities are obtained from quantum states and measurements
using [`measure`](@ref).

A `DomainError` is thrown if:
* `priors` does not satisfy [`is_probability_distribution`](@ref).
* Any `ρ ∈ ρ_states` does not satisfy [`is_density_matrix`](@ref).
* `Π` does not satisfy [`is_povm`](@ref).
"""
function mutual_information(
    priors::AbstractVector,
    states::Vector{<:AbstractMatrix},
    Π::AbstractVector{<:AbstractMatrix}
) :: Float64
    is_probability_distribution(priors) || Probabilities(priors)

    if !all(ρ -> ρ isa State, states)
        states = State.(states)
    end
    if !(Π isa POVM)
        Π = POVM(Π)
    end

    conditionals = measure(Π, states)
    mutual_information(priors, conditionals)
end

"""
    success_probability(
        priors::AbstractVector,
        states::Vector{<:AbstractMatrix},
        Π::Vector{<:AbstractMatrix}
    ) :: Float64

The probability of correctly distinguishing quantum states with the specifed POVM:

```math
P_{\\text{Success}} = \\sum_{i=1}^n p_i \\text{Tr}[\\Pi_i \\rho_i].
```

The number of states must match the number POVMs.
"""
function success_probability(priors::AbstractVector,states::Vector{<:AbstractMatrix}, Π::AbstractVector{<:AbstractMatrix}) :: Float64
    is_probability_distribution(priors) || Probabilities(priors)
    all(is_density_matrix.(states)) || State.(states)
    is_povm(Π) || POVM(Π)

    sum(priors .* tr.(states .* Π))
end

"""
    error_probability(
        priors::AbstractVector,
        ρ_states::Vector{<:AbstractMatrix},
        Π::AbstractVector{<:AbstractMatrix}
    ) :: Float64

The probability of incorrectly distinguishing quantum states with the specifed POVM.
This quantity is simply obtained as ``P_{\\text{Error}} = 1 - P_{\\text{Success}}``.
"""
function error_probability(priors::AbstractVector, ρ_states::Vector{<:AbstractMatrix}, Π::AbstractVector{<:AbstractMatrix}) :: Float64
    1 .- success_probability(priors, ρ_states, Π)
end
