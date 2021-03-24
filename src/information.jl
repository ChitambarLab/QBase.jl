export shannon_entropy, von_neumann_entropy, joint_entropy, conditional_entropy
export holevo_bound, holevo_information
export mutual_information
export success_probability, error_probability

"""
    shannon_entropy( probabilities :: AbstractVector ) :: Float64

The classical entropy of a probability distribution.

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

The von neumann entropy of a density matrix.

A `DomainError` is thrown if `ρ` does not satisyd [`is_density_matrix`](@ref).
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

Computes the upper bound of a quantum channel's information capacity. The information
shared through a quantum channel cannot exceed a classical channel of the same dimension.

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
    holevo_information(
        priors :: AbstractVector,
        ρ_states :: Vector{<:AbstractMatrix,
        Π :: AbstractVector
    ) :: Float64

Computes the holevo (mutual) information shared through a quantum channel.

A `DomainError` is thrown if:
* `priors` does not satisfy [`is_probability_distribution`](@ref).
* Any `ρ ∈ ρ_states` does not satisfy [`is_density_matrix`](@ref).
* `Π` does not satisfy [`is_povm`](@ref).
"""
function holevo_information(
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
    joint_entropy(priors :: AbstractVector, conditionals :: AbstractMatrix) :: Float64

Returns the entropy for the union of pdf ``P(x,y)``.
"""
function joint_entropy(priors::AbstractVector, conditionals::AbstractMatrix) :: Float64
    is_probability_distribution(priors) || Probabilities(priors)
    is_conditional_distribution(conditionals) || Conditionals(conditionals)

    joint_probabilities = JointProbabilities(priors, conditionals)
    shannon_entropy(joint_probabilities[:])
end

"""
    conditional_entropy(priors::AbstractVector, conditionals::AbstractMatrix) :: Float64

Returns the conditional entropy for the system with specified `priors` and `conditionals`.
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

The entropy of the overlap between p(x) and p(y). The information shared from y to x.
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
    success_probability(
        priors::AbstractVector,
        states::Vector{<:AbstractMatrix},
        Π::Vector{<:AbstractMatrix}
    ) :: Float64

The probability of correctly distinguishing quantum states with the specifed POVM.
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
"""
function error_probability(priors::AbstractVector, ρ_states::Vector{<:AbstractMatrix}, Π::AbstractVector{<:AbstractMatrix}) :: Float64
    1 .- success_probability(priors, ρ_states, Π)
end
