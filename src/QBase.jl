"""
A collection of methods and submodules useful for computation of quantum sytems.
"""
module QBase

using LinearAlgebra
using Base.Iterators: flatten
using Combinatorics: permutations, combinations
using RandomMatrices: Haar

# Methods extended by QBase.jl
import Base: convert, *, kron, adjoint, show, sqrt
import LinearAlgebra: rank, eigvals

const ATOL = 1e-7

# math utilities
include("math/matrices.jl")
include("math/combinatorics.jl")

# type exports
include("types/probabilities.jl")
include("types/brakets.jl")
include("types/operators.jl")
include("types/states.jl")
include("types/unitaries.jl")
include("types/measurements.jl")

# TODO: Unitary
# TODO: KrausOperator
# TODO: AbstractChannel
# TODO: AbstractSuperOperator
# TODO: AbstractObservable

include("constructors/brakets.jl")
include("constructors/states.jl")
include("constructors/unitaries.jl")
include("constructors/measurements.jl")

# main methods
include("evolve.jl")
include("measure.jl")
include("channels.jl")
include("information.jl")

end
