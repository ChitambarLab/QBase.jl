"""
A QBase.jl submodule providing general purpose mathematics useful to quantum mechanics.

# Exports:

*Matrices:*
- [`partial_trace`](@ref) - The partial trace matrix operation.
- [`computational_basis_vectors`](@ref) - A general orthonormal set of vectors.

*Combinatorics:*
- [`stirling2`](@ref) - Counts the ways to partition `n` items into `k` sets.
- [`stirling2_partitions`](@ref) - Enumerates the stirling2 partitions.
- [`stirling2_matrices`](@ref) - Generates matrices representing the stirling2 partitions.
- [`permutations`] - Passthrough for `Combinatorics.permutations` ([JuliaMath/Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl)).
- [`permutation_matrices`](@ref) - Generates the set of permutation matrices.
- [`combinations`] - Passthrough for `Combinatorics.combinations` ([JuliaMath/Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl)).
- [`n_choose_k_matrices`](@ref) - Generates matrices representing n-choose-k combinations.
- [`base_n_val`](@ref) - Converts a base-n number into the base-10 value.

*Probability:*
- [`is_probability_distribution`](@ref) - Checks if a vector is probability distribution.
- [`is_conditional_distribution`](@ref) - Checks if a matrix is a conditional probability distribution.
- [`Marginals`](@ref) - Base type for marginal probability distributions.
- [`Conditionals`](@ref) - Base tyep for conditional probability distributions.
"""
module QMath

using LinearAlgebra

include("./QMath/validation.jl")
include("./QMath/matrices.jl")
include("./QMath/probability.jl")
include("./QMath/combinatorics.jl")

end
