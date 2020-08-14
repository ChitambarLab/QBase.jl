"""
A collection of methods and submodules useful for computation of quantum sytems.

# Exports

*Methods:*
- [`evolve`](@ref) - Performs unitary evolution of quantum states.
- [`measurement_probs`](@ref) - Outcome probabilities of quantum measurement.

*Modules:*
- [`States`](@ref) - Types and constructors for representing quantum states.
- [`Observables`](@ref) - Types and constructors for representing measureable quantities.
- [`Unitaries`](@ref) - Types and constructors for representing unitary operators.
- [`Information`](@ref) - Functions for computing information-theoretic quantities.
- [`QMath`](@ref) - Mathematics useful for modeling quantum operations.
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

export evolve, measurement_probs

"""
Apply a unitary evolution `U` to density matrix `ρ`.

    evolve(
        U::Unitaries.AbstractUnitary,
        ρ::States.AbstractDensityMatrix
    ) :: DensityMatrix

Apply a unitary evolution `U` to a state ket `ψ`.

    evolve(
        U::Unitaries.AbstractUnitary,
        ψ::States.AbstractKet
    ) :: Ket
"""
evolve(U::Unitaries.AbstractUnitary, ρ::States.AbstractDensityMatrix) :: DensityMatrix = DensityMatrix(U*ρ*U')
evolve(U::Unitaries.AbstractUnitary, ψ::States.AbstractKet) :: Ket = U*ψ

"""
Computes the outcome probabilities for a quantum measurement. The conditional
probabilities are determined by the Born rule, ``P(i|j) = \\text{Tr}[\\Pi_i \\rho_j]``,
where ``\\Pi_j`` is a POVM element and ``\\rho_j`` is a density matrix.

Measurement of a single `Ket` or  `DensityMatrix`:

    measurement_probs(
        Π :: Observables.AbstractPOVM,
        ρ :: States.AbstractDensityMatrix
    ) :: QMath.Conditionals

    measurement_probs(
        Π :: Observables.AbstractPOVM,
        ψ :: States.AbstractKet
    ) :: QMath.Conditionals

Measurement of an ensemble of `Ket` or `DensityMatrix` types:

    measurement_probs(
        Π :: Observables.AbstractPOVM,
        ρ_states :: Vector{<:States.AbstractDensityMatrix}
    ) :: QMath.Conditionals

    measurement_probs(
        Π :: Observables.AbstractPOVM,
        ψ_kets :: Vector{<:States.ABstractKet}
    ) :: QMath.Conditionals
"""
function measurement_probs(
    Π :: Observables.AbstractPOVM,
    ρ :: States.AbstractDensityMatrix
) :: QMath.Marginals
    QMath.Marginals( real.(map(Π_el -> tr(Π_el * ρ), Π)) )
end
function measurement_probs(
    Π :: Observables.AbstractPOVM,
    ρ_states :: Vector{<:States.AbstractDensityMatrix}
) :: QMath.Conditionals
    QMath.Conditionals( real.(tr.(transpose(ρ_states * transpose(Π))) ))
end
function measurement_probs(
    Π :: Observables.AbstractPOVM,
    ψ :: States.AbstractKet
) :: QMath.Marginals
    QMath.Marginals( real.(map( Π_el -> ψ' * Π_el * ψ, Π)) )
end
function measurement_probs(
    Π :: Observables.AbstractPOVM,
    ψ_kets :: Vector{<:States.AbstractKet}
) :: QMath.Conditionals
    measurement_probs(Π, States.pure_state.(ψ_kets))
end

end
