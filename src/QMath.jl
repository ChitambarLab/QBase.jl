"""
A QBase.jl submodule providing general purpose mathematics useful to quantum mechanics.

# Exports:

*Matrices:*
- [`partial_trace`](@ref) - The partial trace matrix operation.
- [`computational_basis_vectors`](@ref) - A general orthonormal set of vectors.
"""
module QMath

using LinearAlgebra

include("./QMath/validation.jl")

export partial_trace, computational_basis_vectors
include("./QMath/matrices.jl")


include("./QMath/probability.jl")
include("./QMath/enumerate.jl")

end
