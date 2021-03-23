export shannon_entropy, von_neumann_entropy, joint_entropy, conditional_entropy
export holevo_bound, holevo_information
export mutual_information
export success_probability, error_probability

"""
    shannon_entropy( probabilities :: Probabilities ) :: Float64

    shannon_entropy( probabilities :: Vector{<:Real}) :: Float64

The classical entropy of a probability distribution.
"""
function shannon_entropy( probabilities :: Probabilities ) :: Float64
    entropy = -1*sum(map(
        (p) -> isapprox(p, 0, atol=1e-7) ? 0 : p*log2(p),
        probabilities
    ))

    entropy
end
shannon_entropy(probs::Vector{<:Real}) :: Float64 = shannon_entropy(Probabilities(probs))

"""
    von_neumann_entropy( ρ :: AbstractState ) :: Float64

    von_neumann_entropy( ρ :: Matrix{<:Number} ) :: Float64

The von neumann entropy of a density matrix.
"""
function von_neumann_entropy(ρ::AbstractState) :: Float64
    λs = Probabilities(eigvals(ρ))
    entropy = shannon_entropy(λs)

    entropy
end
von_neumann_entropy(ρ::Matrix{<:Number}) :: Float64 = von_neumann_entropy(State(ρ))

"""
    holevo_bound(
        priors :: Probabilities,
        ρ_states :: Vector{<:AbstractState}
    ) :: Float64

    holevo_bound(
        priors :: Vector{<:Real},
        ρ_states :: Vector{Matrix{<:Number}}
    ) :: Float64

Computes the upper bound of a quantum channel's information capacity. The information
shared through a quantum channel cannot exceed a classical channel of the same dimension.
"""
function holevo_bound(priors::Probabilities, ρ_states::Vector{<:AbstractState}) :: Float64
    ρ = mixed_state(priors, ρ_states)
    ρ_ent = von_neumann_entropy(ρ)

    ρ_states_ent = von_neumann_entropy.(ρ_states)
    ρ_sum_ent = sum( priors .* ρ_states_ent )

    bound = ρ_ent - ρ_sum_ent

    bound
end
holevo_bound(priors::AbstractVector{<:Real}, ρ_states::Vector{<:AbstractMatrix{T}} where T <: Number ) :: Float64 = begin
    holevo_bound(
        (priors isa Probabilities) ? priors : Probabilities(priors),
        (ρ_states isa Vector{<:AbstractState}) ? ρ_states : State.(ρ_states)
    )
end

"""
    holevo_information(
        priors :: Probabilities,
        ρ_states :: Vector{<:AbstractState},
        Π :: POVM
    ) :: Float64

Computes the holevo (mutual) information shared through a quantum channel.
"""
function holevo_information(
    priors::Probabilities,
    ρ_states::Vector{<:AbstractState},
    Π::POVM
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

    conditionals = Conditionals(joint_probs .* (1 ./ priors)')

    h_information = mutual_information(priors, conditionals)

    h_information
end
holevo_information(
    priors::AbstractVector{<:Real},
    ρ_states::Vector{<:AbstractMatrix{T}} where T <: Number,
    Π :: AbstractVector{Matrix{S}} where S <: Number
) = begin
    holevo_information(
        (priors isa Probabilities) ? priors : Probabilities(priors),
        (ρ_states isa Vector{<:AbstractState}) ? ρ_states : State.(ρ_states),
        (Π isa POVM) ? Π : POVM(Π)
    )
end

"""
    joint_entropy(priors :: Probabilities, conditionals :: Conditionals) :: Float64

Returns the entropy for the union of pdf ``P(x,y)``.
"""
function joint_entropy(priors::Probabilities, conditionals::Conditionals) :: Float64
    joint_probabilities = joint_probabilities(priors, conditionals)

    shannon_entropy(Probabilities(joint_probabilities[:]))
end
joint_entropy(
    priors::AbstractVector{<:Number},
    conditionals::AbstractMatrix{<:Number}
) :: Float64 = begin
    joint_entropy(
        (priors isa Probabilities) ? priors : Probabilities(priors),
        (conditionals isa Conditionals) ? conditionals : Conditionals(conditionals)
    )
end

"""
    conditional_entropy(priors::Probabilities, conditionals::Conditionals) :: Float64

Returns the conditional entropy for the system with specified `priors` and `conditionals`.
"""
function conditional_entropy(priors::Probabilities, conditionals::Conditionals) :: Float64
    joint_probabilities = joint_probabilities(priors, conditionals)

    (num_row, num_col) = size(joint_probabilities)

    conditional_entropy = shannon_entropy(Probabilities(joint_probabilities[:])) + sum(sum(map( row -> map( col ->
        joint_probabilities[row,col]*( isapprox(priors[col], 0, atol=1e-7) ? 0 : log2(priors[col]) ),
    1:num_col), 1:num_row)))

    conditional_entropy
end
conditional_entropy(
    priors::AbstractVector{<:Number},
    conditionals::AbstractMatrix{<:Number}
) :: Float64 = begin
    conditional_entropy(
        (priors isa Probabilities) ? priors : Probabilities(priors),
        (conditionals isa Conditionals) ? conditionals : Conditionals(conditionals)
    )
end

"""
    mutual_information(
        priors :: Probabilities,
        conditionals :: Conditionals
    ) :: Float64

The entropy of the overlap between p(x) and p(y). The information shared from y to x.
"""
function mutual_information(priors::Probabilities, conditionals::Conditionals) :: Float64

    prior_ent = shannon_entropy(priors)

    outcome_probs = outcome_probabilities(priors,conditionals)
    outcome_ent = shannon_entropy( Probabilities(outcome_probs) )
    joint_ent = joint_entropy(priors, conditionals)

    mutual_information = prior_ent + outcome_ent - joint_ent

    mutual_information
end
mutual_information(
    priors::AbstractVector{<:Number},
    conditionals::AbstractMatrix{<:Number}
) :: Float64 = begin
    mutual_information(
        (priors isa Probabilities) ? priors : Probabilities(priors),
        (conditionals isa Conditionals) ? conditionals : Conditionals(conditionals)
    )
end

"""
    success_probability(
        priors::Probabilities,
        ρ_states::Vector{<:AbstractState},
        Π::POVM
    ) :: Float64

The probability of correctly distinguishing quantum states with the specifed POVM.
"""
function success_probability(priors::Probabilities, ρ_states::Vector{<:AbstractState}, Π::POVM) :: Float64
    sum(priors .* tr.(ρ_states .* Π))
end
success_probability(
    priors::AbstractVector{<:Real},
    ρ_states::Vector{<:AbstractMatrix{T}} where T <: Number,
    Π :: AbstractVector{Matrix{S}} where S <: Number
) = begin
    success_probability(
        (priors isa Probabilities) ? priors : Probabilities(priors),
        (ρ_states isa Vector{<:AbstractState}) ? ρ_states : State.(ρ_states),
        (Π isa POVM) ? Π : POVM(Π)
    )
end

"""
    error_probability(
        priors::Probabilities,
        ρ_states::Vector{<:AbstractState},
        Π::POVM
    ) :: Float64

The probability of incorrectly distinguishing quantum states with the specifed POVM.
"""
function error_probability(priors::Probabilities, ρ_states::Vector{<:AbstractState}, Π::POVM) :: Float64
    1 .- success_probability(priors, ρ_states, Π)
end
error_probability(
    priors::AbstractVector{<:Real},
    ρ_states::Vector{<:AbstractMatrix{T}} where T <: Number,
    Π :: AbstractVector{Matrix{S}} where S <: Number
) = begin
    error_probability(
        (priors isa Probabilities) ? priors : Probabilities(priors),
        (ρ_states isa Vector{<:AbstractState}) ? ρ_states : State.(ρ_states),
        (Π isa POVM) ? Π : POVM(Π)
    )
end
