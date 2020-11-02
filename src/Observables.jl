"""
Observables describe measurable quantities of a quantum system. Quantum
measurement is a probabilistic process. The measurement outcomes are not definite,
but described by conditional probabilities.

The `QBase.Observables` submodule provides types and constructors for representing
quantum observables.
"""
module Observables

using LinearAlgebra
using ..QMath
using ..States

# validation
export is_povm

# types
export AbstractPOVM, POVM, QubitPOVM

# constructors
export planar_symmetric_qubit_povm
export mirror_symmetric_qubit_3povm, asymmetric_qubit_3povm, sqrt_povm, sic_qubit_povm, trine_qubit_povm

# methods
export kraus_operators, naimark_dilation

"""
    AbstractPOVM <: AbstractMatrix{Complex{Float64}}

The abstract type representing positive-operator valued measures (POVMs). An
`AbstractPOVM` cannot be instantiated, it serves as a supertype from which concrete
types are derived.
"""
abstract type AbstractPOVM <: AbstractVector{Matrix{Complex{Float64}}} end
Base.size(P::AbstractPOVM) = size(P.Π)
Base.getindex(P::AbstractPOVM, id::Int) = getindex(P.Π, id...)
Base.setindex!(P::AbstractPOVM, v, id::Int) = (P.Π[id...] = v)

"""
    is_povm( Π :: Vector ) :: Bool

Returns `true` if `Π` is a POVM. The following constraints must be satisfied:
* Each POVM element is hermitian
* Each POVM element positive semi-definite
* The POVM is complete: `sum(Π) == I`
"""
is_povm(Π::AbstractPOVM) :: Bool = true
function is_povm(Π::Vector) :: Bool
    dim = size(Π[1])[1]

    all_positive_semidefinite = all(QMath.is_positive_semidefinite, Π)
    all_hermitian = all(QMath.is_hermitian, Π)
    is_complete = isapprox(sum(Π), Matrix{Complex{Float64}}(I,dim,dim), atol=10e-5)

    (all_hermitian & all_positive_semidefinite & is_complete)
end

"""
    POVM( Π :: Vector{Matrix} ) <: AbstractPOVM

Positve-operator valued measures (POVMs) represent a general quantum measurement.
Each POVM-element is a hermitian, positive-semidefinite matrix. The sum of all
POVM-elements yields the identity matrix. The constructor, `POVM(Π)` throws a
`DomainError` if the provided array of matrices, `Π` is not a valid POVM.
"""
struct POVM <: AbstractPOVM
    Π::Vector{Matrix{Complex{Float64}}}
    POVM(Π) = is_povm(Π) ? new(Π) : throw(DomainError(Π, "povm Π is invalid"))
end

"""
    QubitPOVM( Π :: Vector{Matrix} ) <: AbstractPOVM

A general qubit measurement. A `DomainError` is thrown if `Π` does not contain
2x2 elements or if it is not a valid POVM.
"""
struct QubitPOVM <: AbstractPOVM
    Π::POVM
    QubitPOVM(Π) = (size(Π[1]) == (2,2)) ? new(POVM(Π)) : throw(DomainError(Π, "POVM Π does not contain 2x2 elements"))
end

"""
    mirror_symmetric_qubit_3povm( θ :: Real ) :: QubitPOVM

Constructs a `QubitPOVM` aligned with the three mirror symmetric qubit states.
The first measurement is aligned with the ``|0\\rangle`` state and the remaining
two are symmetric across the z-axis.

A `DomainError` is thrown if argument `θ` ∉ [π/4,π/2].
"""
function mirror_symmetric_qubit_3povm(θ::Real)::QubitPOVM
    if !( π/4 <= θ <= π/2)
        throw(DomainError(θ, "angle must be in range [π/4,π/2]"))
    end

    # POVM weights
    γ1 = 1 - cot(θ)^2
    γ2 = 0.5*(1 + cot(θ)^2)

    mirror_symmetric_qubits = States.mirror_symmetric_qubits(θ)

    π1 = γ1*mirror_symmetric_qubits[1]
    π2 = γ2*mirror_symmetric_qubits[2]
    π3 = γ2*mirror_symmetric_qubits[3]

    QubitPOVM([π1, π2, π3])
end

"""
    asymmetric_qubit_3povm( θ1::Real, θ2::Real ) :: QubitPOVM

Constructs a general non-orthogonal 3-element `QubitPOVM`.

Inputs:
* `θ1 ∈ (0,π/2] or [-π/2,0)`, the angle element 2 makes with ``|0\\rangle`` state.
* `θ2 ∈ [-π/2,0) or (0,π/2]`, the angle element 3 makes with ``|0\\rangle`` state.

A `DomainError` is thrown if `θ1` and `θ2` are not in the valid ranges. Note that
one angle must be positive and the other negative.
"""
function asymmetric_qubit_3povm(θ1::Real,θ2::Real)::QubitPOVM

    if !(((0 < θ1 <= π/2) && (-π/2 <= θ2 < 0)) || ((0 < θ2 <= π/2) && (-π/2 <= θ1 < 0)))
        throw(DomainError((θ1,θ2), "angles (θ1, θ2) are not in the valid range"))
    end

    # povm scaling terms derived analytically
    γ2 = 1/(sin(θ1)^2 - sin(θ1)*cos(θ1)*tan(θ2))
    γ3 = -γ2*(sin(θ1)*cos(θ1))/(sin(θ2)*cos(θ2))
    γ1 = 1 - γ2*cos(θ1)^2 - γ3*cos(θ2)^2

    π1 = γ1*[1 0;0 0]
    π2 = γ2*[cos(θ1)^2 cos(θ1)*sin(θ1); cos(θ1)*sin(θ1) sin(θ1)^2]
    π3 = γ3*[cos(θ2)^2 cos(θ2)*sin(θ2); cos(θ2)*sin(θ2) sin(θ2)^2]

    QubitPOVM([π1, π2, π3])
end

"""
    planar_symmetric_qubit_povm( n :: Int64 ) :: QubitPOVM

Constructs an `n`-element `QubitPOVM` from the [`planar_symmetric_qubit_states`](@ref).
Each state is multipled by a factor of `2/n` to satisfy the completeness relation.

A `DomainError` is thrown if `n ≥ 2` is not satisfied.
"""
function planar_symmetric_qubit_povm(n :: Int64) :: QubitPOVM
    QubitPOVM(2/n*States.planar_symmetric_qubits(n))
end

"""
    sqrt_povm(priors :: QMath.Marginals, states :: Vector{<:States.AbstractDensityMatrix}) :: POVM

Returns the "pretty good" square-root povm for the given density operator
states and priors.
"""
function sqrt_povm(priors :: QMath.Marginals, states :: Vector{<:States.AbstractDensityMatrix}) :: POVM

    ρ_mix = sum(priors .* states)
    ρ_sqrt = sqrt(inv(ρ_mix))

    Π = map(
        (ρ) -> ρ_sqrt * ρ * ρ_sqrt,
        (priors .* states)
    )

    POVM(Π)
end

"""
    sic_qubit_povm :: QubitPOVM

The POVM with elements parallel to the symmetric informationally complete (SIC) qubits.
"""
const sic_qubit_povm = QubitPOVM(0.5 * States.sic_qubits)

"""
    trine_qubit_povm :: QubitPOVM

The POVM with elements parallel to the trine qubit states.
"""
const trine_qubit_povm = QubitPOVM(2/3 * States.trine_qubits)

"""
    kraus_operators(Π::AbstractPOVM) :: Array{Array{Complex{Float64},2},1}

Returns the Kraus operators for the provided POVM. In general, the Kraus operators
form a continuum and are non-unique. In this method, the construction

``k_i = \\sqrt{\\Pi_i)}\\otimes |i\\rangle``
"""
function kraus_operators(Π::AbstractPOVM) :: Array{Array{Complex{Float64},2},1}
    map(i -> kron(sqrt(Π[i]), Matrix(1I,3,3)[:,i]), 1:length(Π))
end

"""
    naimark_dilation( Π::AbstractPOVM )

Returns the dilated projectors which are equivalent to the provided POVM. During
measurement, the measured state must be tensored with the ancilla.
"""
function naimark_dilation(Π::AbstractPOVM)
    n = length(Π)
    d = size(Π[1])[1]

    ancilla = zeros(n,n)
    ancilla[1,1] = 1

    k = kraus_operators(Π)
    V = sum(k)

    orth_comp = nullspace(V')

    U_cols = []
    null_id = 1
    for i in 1:n*d
        if floor((i-1)/n) == (i-1)/n
            push!(U_cols,V[:,convert(Int,(i-1)/n + 1)])
        else
            push!(U_cols, orth_comp[:,null_id])
            null_id = null_id + 1
        end

    end

    U = hcat(U_cols...)

    p_ancilla = map( i -> begin
        p = zeros(n,n)
        p[i,i] = 1

        return p
    end, 1:n )

    projectors = map(i -> U'*kron(Matrix(1I,d,d),p_ancilla[i])*U, 1:n)

    Dict(
        "projectors" => projectors,
        "ancilla" => ancilla
    )
end

end
