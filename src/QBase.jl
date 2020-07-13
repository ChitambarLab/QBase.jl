"""
The main module for QBase.jl provides a collection of submodules which are useful
for computations of quantum sytems.

# Exports

*Modules:*
- [`QMath`](@ref) - Mathematics useful for modeling quantum operations.
- [`States`](@ref) - Types and constructors for representing quantum states.
- [`Observables`](@ref) - Types and constructors for representing measureable quantities.
- [`Unitaries`](@ref) - Types and constructors for representing unitary operators.
- [`Information`](@ref) - Functions for computing information-theoretic quantities.
"""
module QBase

using LinearAlgebra

# submodules are exported
export QMath, Unitaries, States, Observables, Information

# include modules
include("./QMath.jl")
using .QMath

include("./Unitaries.jl")
using .Unitaries

include("./States.jl")
using .States

include("./Observables.jl")
using .Observables

include("./Information.jl")
using .Information


export evolve

"""
    evolve(
        U::Unitaries.AbstractUnitary,
        ρ::States.AbstractDensityMatrix
    ) :: DensityMatrix

Apply a unitary evolution `U` to density matrix `ρ`.

    evolve(
        U::Unitaries.AbstractUnitary,
        ψ::States.AbstractKet
    ) :: Ket

Apply a unitary evolution `U` to a state ket `ψ`.
"""
evolve(U::Unitaries.AbstractUnitary, ρ::States.AbstractDensityMatrix) :: DensityMatrix = DensityMatrix(U*ρ*U')
evolve(U::Unitaries.AbstractUnitary, ψ::States.AbstractKet) :: Ket = U*ψ

end
