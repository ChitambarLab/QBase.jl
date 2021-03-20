"""
Methods for quantifying information and randomness.

Entropy quantifiers are taken with respect to base-2 logarithms. The entropy is
understood as the number of bits {0,1} required to communicate a random result with
certainty.


"""
module Information

using ..QMath
using ..States
using ..Observables

using LinearAlgebra

"""
    shannon_entropy( probabilities :: QMath.Marginals ) :: Float64

    shannon_entropy( probabilities :: Vector{<:Real}) :: Float64

The classical entropy of a probability distribution.
"""
function shannon_entropy( probabilities :: QMath.Marginals ) :: Float64
    entropy = -1*sum(map(
        (p) -> isapprox(p, 0, atol=1e-7) ? 0 : p*log2(p),
        probabilities
    ))

    entropy
end
shannon_entropy(probs::Vector{<:Real}) :: Float64 = shannon_entropy(QMath.Marginals(probs))

"""
    von_neumann_entropy( ρ :: States.AbstractDensityMatrix ) :: Float64

    von_neumann_entropy( ρ :: Matrix{<:Number} ) :: Float64

The von neumann entropy of a density matrix.
"""
function von_neumann_entropy(ρ::States.AbstractDensityMatrix) :: Float64
    λs = QMath.Marginals(eigvals(ρ))
    entropy = Information.shannon_entropy(λs)

    entropy
end
von_neumann_entropy(ρ::Matrix{<:Number}) :: Float64 = von_neumann_entropy(States.DensityMatrix(ρ))

"""
    holevo_bound(
        priors :: QMath.Marginals,
        ρ_states :: Vector{<:AbstractDensityMatrix}
    ) :: Float64

    holevo_bound(
        priors :: Vector{<:Real},
        ρ_states :: Vector{Matrix{<:Number}}
    ) :: Float64

Computes the upper bound of a quantum channel's information capacity. The information
shared through a quantum channel cannot exceed a classical channel of the same dimension.
"""
function holevo_bound(priors::QMath.Marginals, ρ_states::Vector{<:States.AbstractDensityMatrix}) :: Float64
    ρ = States.mixed_state(priors, ρ_states)
    ρ_ent = Information.von_neumann_entropy(ρ)

    ρ_states_ent = Information.von_neumann_entropy.(ρ_states)
    ρ_sum_ent = sum( priors .* ρ_states_ent )

    bound = ρ_ent - ρ_sum_ent

    bound
end
holevo_bound(priors::AbstractVector{<:Real}, ρ_states::Vector{<:AbstractMatrix{T}} where T <: Number ) :: Float64 = begin
    holevo_bound(
        (priors isa QMath.Marginals) ? priors : QMath.Marginals(priors),
        (ρ_states isa Vector{<:States.AbstractDensityMatrix}) ? ρ_states : States.DensityMatrix.(ρ_states)
    )
end

"""
    holevo_information(
        priors :: QMath.Marginals,
        ρ_states :: Vector{<:AbstractDensityMatrix},
        Π :: Observables.AbstractPOVM
    ) :: Float64

Computes the holevo (mutual) information shared through a quantum channel.
"""
function holevo_information(
    priors::QMath.Marginals,
    ρ_states::Vector{<:States.AbstractDensityMatrix},
    Π::Observables.AbstractPOVM
) :: Float64
    x_num = length(priors)
    y_num = length(Π)

    joint_probs = zeros(y_num,x_num)
    for x_id in 1:x_num
        p = priors[x_id]
        ρ = ρ_states[x_id]
        for y_id in 1:y_num
            joint_probs[y_id,x_id] += p*tr(ρ*Π[y_id])
        end
    end

    conditionals = QMath.Conditionals(joint_probs .* (1 ./ priors)')

    h_information = Information.mutual_information(priors, conditionals)

    h_information
end
holevo_information(
    priors::AbstractVector{<:Real},
    ρ_states::Vector{<:AbstractMatrix{T}} where T <: Number,
    Π :: AbstractVector{Matrix{S}} where S <: Number
) = begin
    holevo_information(
        (priors isa QMath.Marginals) ? priors : QMath.Marginals(priors),
        (ρ_states isa Vector{<:States.AbstractDensityMatrix}) ? ρ_states : States.DensityMatrix.(ρ_states),
        (Π isa Observables.AbstractPOVM) ? Π : Observables.POVM(Π)
    )
end

"""
    joint_entropy(priors :: QMath.Marginals, conditionals :: QMath.Conditionals) :: Float64

Returns the entropy for the union of pdf ``P(x,y)``.
"""
function joint_entropy(priors::QMath.Marginals, conditionals::QMath.Conditionals) :: Float64
    joint_probabilities = QMath.joint_probabilities(priors, conditionals)

    shannon_entropy(QMath.Marginals(joint_probabilities[:]))
end
joint_entropy(
    priors::AbstractVector{<:Number},
    conditionals::AbstractMatrix{<:Number}
) :: Float64 = begin
    joint_entropy(
        (priors isa QMath.Marginals) ? priors : QMath.Marginals(priors),
        (conditionals isa QMath.Conditionals) ? conditionals : QMath.Conditionals(conditionals)
    )
end

"""
    conditional_entropy(priors::QMath.Marginals, conditionals::QMath.Conditionals) :: Float64

Returns the conditional entropy for the system with specified `priors` and `conditionals`.
"""
function conditional_entropy(priors::QMath.Marginals, conditionals::QMath.Conditionals) :: Float64
    joint_probabilities = QMath.joint_probabilities(priors, conditionals)

    (num_row, num_col) = size(joint_probabilities)

    conditional_entropy = shannon_entropy(QMath.Marginals(joint_probabilities[:])) + sum(sum(map( row -> map( col ->
        joint_probabilities[row,col]*( isapprox(priors[col], 0, atol=1e-7) ? 0 : log2(priors[col]) ),
    1:num_col), 1:num_row)))

    conditional_entropy
end
conditional_entropy(
    priors::AbstractVector{<:Number},
    conditionals::AbstractMatrix{<:Number}
) :: Float64 = begin
    conditional_entropy(
        (priors isa QMath.Marginals) ? priors : QMath.Marginals(priors),
        (conditionals isa QMath.Conditionals) ? conditionals : QMath.Conditionals(conditionals)
    )
end

"""
    mutual_information(
        priors :: QMath.Marginals,
        conditionals :: QMath.Conditionals
    ) :: Float64

The entropy of the overlap between p(x) and p(y). The information shared from y to x.
"""
function mutual_information(priors::QMath.Marginals, conditionals::QMath.Conditionals) :: Float64

    prior_ent = Information.shannon_entropy(priors)

    outcome_probs = QMath.outcome_probabilities(priors,conditionals)
    outcome_ent = Information.shannon_entropy( QMath.Marginals(outcome_probs) )
    joint_ent = Information.joint_entropy(priors, conditionals)

    mutual_information = prior_ent + outcome_ent - joint_ent

    mutual_information
end
mutual_information(
    priors::AbstractVector{<:Number},
    conditionals::AbstractMatrix{<:Number}
) :: Float64 = begin
    mutual_information(
        (priors isa QMath.Marginals) ? priors : QMath.Marginals(priors),
        (conditionals isa QMath.Conditionals) ? conditionals : QMath.Conditionals(conditionals)
    )
end

"""
    success_probability(
        priors::QMath.Marginals,
        ρ_states::Vector{<:States.AbstractDensityMatrix},
        Π::Observables.AbstractPOVM
    ) :: Float64

The probability of correctly distinguishing quantum states with the specifed POVM.
"""
function success_probability(priors::QMath.Marginals, ρ_states::Vector{<:States.AbstractDensityMatrix}, Π::Observables.AbstractPOVM) :: Float64
    sum(priors .* tr.(ρ_states .* Π))
end
success_probability(
    priors::AbstractVector{<:Real},
    ρ_states::Vector{<:AbstractMatrix{T}} where T <: Number,
    Π :: AbstractVector{Matrix{S}} where S <: Number
) = begin
    success_probability(
        (priors isa QMath.Marginals) ? priors : QMath.Marginals(priors),
        (ρ_states isa Vector{<:States.AbstractDensityMatrix}) ? ρ_states : States.DensityMatrix.(ρ_states),
        (Π isa Observables.AbstractPOVM) ? Π : Observables.POVM(Π)
    )
end

"""
    error_probability(
        priors::QMath.Marginals,
        ρ_states::Vector{<:States.AbstractDensityMatrix},
        Π::Observables.AbstractPOVM
    ) :: Float64

The probability of incorrectly distinguishing quantum states with the specifed POVM.
"""
function error_probability(priors::QMath.Marginals, ρ_states::Vector{<:States.AbstractDensityMatrix}, Π::Observables.AbstractPOVM) :: Float64
    1 .- success_probability(priors, ρ_states, Π)
end
error_probability(
    priors::AbstractVector{<:Real},
    ρ_states::Vector{<:AbstractMatrix{T}} where T <: Number,
    Π :: AbstractVector{Matrix{S}} where S <: Number
) = begin
    error_probability(
        (priors isa QMath.Marginals) ? priors : QMath.Marginals(priors),
        (ρ_states isa Vector{<:States.AbstractDensityMatrix}) ? ρ_states : States.DensityMatrix.(ρ_states),
        (Π isa Observables.AbstractPOVM) ? Π : Observables.POVM(Π)
    )
end

end
