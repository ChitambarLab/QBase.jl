"""
A QBase.jl submodule providing general purpose mathematics useful to quantum mechanics.

# Exports:

*Matrices:*
- [`partial_trace`](@ref) - The partial trace matrix operation.
- [`computational_basis_vectors`](@ref) - A general orthonormal set of vectors.
- [`commutes`](@ref) - Checks if two matrices commute.
- [`is_hermitian`](@ref) - Checks if a matrix is hermitian.
- [`is_positive_semidefinite`](@ref) - Checks if a matrix is positive semi-definite.
- [`is_square`](@ref) - Checks if a matrix is square.

*Combinatorics:*
- [`stirling2`](@ref) - Counts the ways to partition `n` items into `k` sets.
- [`stirling2_partitions`](@ref) - Enumerates the stirling2 partitions.
- [`stirling2_matrices`](@ref) - Generates matrices representing the stirling2 partitions.
- [`permutations`](https://github.com/JuliaMath/Combinatorics.jl) - Passthrough for `Combinatorics.permutations`.
- [`permutation_matrices`](@ref) - Generates the set of permutation matrices.
- [`combinations`](https://github.com/JuliaMath/Combinatorics.jl) - Passthrough for `Combinatorics.combinations`.
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

include("./QMath/matrices.jl")
include("./QMath/probability.jl")
include("./QMath/combinatorics.jl")

end