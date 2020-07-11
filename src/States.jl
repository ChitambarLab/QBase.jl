"""
The `QBase.States` submodule provides:

* `abstract` and `concrete` types for quantum state representations.
* A catalog of constructors for instantiating quantum states.
"""
module States

using ..QMath
using ..Unitaries

using LinearAlgebra

# Types
export AbstractKet, AbstractDensityMatrix, Ket, QubitKet, DensityMatrix, Qubit

# Validation
export is_ket, is_density_matrix

# State Constructors
export bloch_qubit_ket, bloch_qubit, pure_state, pure_qubit, mixed_state, mixed_qubit

export bell_kets, basis_kets, basis_states

# 3 State Array Constructors
export mirror_symmetric_qubit_kets, mirror_symmetric_qubits, trine_qubit_kets, trine_qubits

# 4 State Array Constructors
export bb84_qubits, sic_qubits

"""
    AbstractKet <: AbstractVector{Complex{Float64}}

The abstract type representing a quantum state ket. Since kets are contained within
a complex-valued hilbert space, they are appropriately `AbstractVector{Complex{Float64}}`
subtypes. An `AbstractKet` cannot be instantiated, it serves as a supertype from
which ket types are defined.
"""
abstract type AbstractKet <: AbstractVector{Complex{Float64}} end
Base.size(S::AbstractKet) = size(S.ψ)
Base.getindex(S::AbstractKet, I::Int) = getindex(S.ψ, I...)
Base.setindex!(S::AbstractKet, v, I::Int) = (S.ψ[I...] = v)

"""
    is_ket( ψ :: Vector ) :: Bool

Returns `true` if vector `ψ` is a valid ket representation of a quamtum state:

* `ψ` is a real or complex-valued vector.
* `ψ` is normalized with respect to the bra-ket inner prodcut (`ψ'ψ == 0`).
"""
function is_ket(ψ::Vector)::Bool
    norm(ψ) ≈ 1.0
end

"""
    Ket( ψ :: Vector{Complex{Float64}} ) <: AbstractKet

A ket representation of a general quantum state. When given invalid input, the
constructor, `Ket(ψ)`, throws:
* `DomainError` -- If `ψ` is not normalized.
* `MethodError` -- If `ψ` is not a column vector (`[a,b,c]` or `[a;b;c]`)
"""
struct Ket <: AbstractKet
    ψ :: Vector{Complex{Float64}}
    Ket(ψ) = is_ket(ψ) ? new(ψ) : throw(DomainError(ψ, "vector ψ is not normalized"))
end

"""
    QubitKet( ψ :: Vector{Complex{Float64}} ) <: AbstractKet

A ket representation of a 2-dimensional quantum state. When given invalid input,
the constructor, `QubitKet(ψ)`, throws:
* `DomainError` -- If `ψ` is not normalized.
* `MethodError` -- If `ψ` is not a column vector (`[a,b]` or `[a;b]`).
"""
struct QubitKet <: AbstractKet
    ψ :: Ket
    QubitKet(ψ) = (size(ψ) == (2,)) ? new(Ket(ψ)) : throw(DomainError(ψ, "vector ψ is not length 2"))
end

"""
    is_density_matrix( ρ :: Matrix ) :: Bool

Returns true if input `\rho` is:
* Hermitian
* Positive Semi-Definite
* Trace[ρ] = 1 (normalization)
"""
function is_density_matrix(ρ::Matrix)::Bool
    if !(QMath.is_square_matrix(ρ))
        throw(DomainError(ρ, "provided matrix ρ is not square"))
    end

    is_hermitian = QMath.is_hermitian(ρ)
    is_pos_sd = QMath.is_positive_semidefinite(ρ)
    is_trace_one = (tr(ρ) ≈ 1)

    (is_hermitian & is_pos_sd & is_trace_one)
end

"""
    AbstractDensityMatrix <: AbstractMatrix{Complex{Float64}}

The abstract type representing all density matrices.
"""
abstract type AbstractDensityMatrix <: AbstractMatrix{Complex{Float64}} end
Base.size(S::AbstractDensityMatrix) = size(S.ρ)
Base.getindex(S::AbstractDensityMatrix, I::Vararg{Int,2}) = getindex(S.ρ, I...)
Base.setindex!(S::AbstractDensityMatrix, v, I::Vararg{Int,2}) = (S.ρ[I...] = v)

"""
    DensityMatrix( ρ :: Matrix{Complex{Float64}} ) <: AbstractDensityMatrix

The density matrix representation of a quantum state. The constructor, `DensityMatrix(ρ)`
throws a `DomainError` if `is_density_matrix(ρ)` is `false`.

Base methods extended to use the `DensityMatrix` type:
* `QMath.partial_trace` - Returns `DensityMatrix` if supplied with one.
* `Base.kron` - The kronecker product of two density matrices is a `DensityMatrix`.
"""
struct DensityMatrix <: AbstractDensityMatrix
    ρ :: Matrix{Complex{Float64}}
    DensityMatrix(ρ) = is_density_matrix(ρ) ? new(ρ) : throw(DomainError(ρ, "matrix ρ is invalid"))
end
Base.kron(ρ1::AbstractDensityMatrix, ρ2::AbstractDensityMatrix) = DensityMatrix(kron(ρ1.ρ,ρ2.ρ))
QMath.partial_trace(ρ::DensityMatrix, subsystem_dims::Vector{Int64}, subsystem_id::Int64) = begin
    DensityMatrix(QMath.partial_trace(ρ.ρ, subsystem_dims, subsystem_id))
end

"""
    Qubit( ρ :: Matrix{Complex{Float64}} ) <: AbstractDensityMatrix

The 2x2 density matrix representation of a qubit. The constructor, `Qubit(ρ)`
throws a `DomainError` if `is_density_matrix(ρ)` is false or if `ρ` is not 2x2 in
dimension.
"""
struct Qubit <: AbstractDensityMatrix
    ρ :: Matrix{Complex{Float64}}
    Qubit(ρ) = (size(ρ) == (2,2)) ? new(DensityMatrix(ρ)) : throw(DomainError(ρ, "matrix is not 2x2"))
end

"""
    pure_state( ψ :: AbstractKet ) :: Qubit

A state is considered "pure" if it is rank-one. A rank-one density matrix
is constructed by  taking the outer-product of a ket state.

The method alternatively accepts a `Vector` input.

    pure_state( ψ :: Vector ) :: Qubit

A `DomainError` is thrown if `ψ` is not a valid ket.`
"""
function pure_state(ψ::AbstractKet)::DensityMatrix
    DensityMatrix(ψ*ψ')
end

function pure_state(ψ::Vector)::DensityMatrix
    if !(is_ket(ψ))
        throw(DomainError(ψ, "ψ is not a valid wavefunction"))
    end

    DensityMatrix(ψ*ψ')
end

"""
    basis_kets( dim :: Int64 ) :: Vector{Ket}

The computation basis vectors for a hilbert space of dimension, `dim`.
"""
function basis_kets(dim::Int64)::Vector{Ket}
    Ket.(QMath.computational_basis_vectors(dim))
end

"""
    basis_states( dim :: Int64 ) :: Vector{DensityMatrix}

The density matrices for the computational basis of dimension, `dim`.
"""
function basis_states(dim::Int64)::Vector{DensityMatrix}
    pure_state.(QMath.computational_basis_vectors(dim))
end

@doc raw"""
    bell_kets :: Vector{Ket}

The Bell basis kets, ordered as ``\{|\Phi^+\rangle, |\Phi^-\rangle, |\Psi^+\rangle, |\Psi^-\rangle \}``, where

```math
\begin{matrix}
    |\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle), &
    |\Phi^-\rangle = \frac{1}{\sqrt{2}}(|00\rangle - |11\rangle), \\
    |\Psi^+\rangle = \frac{1}{\sqrt{2}}(|01\rangle + |10\rangle), &
    |\Psi^-\rangle = \frac{1}{\sqrt{2}}(|01\rangle - |10\rangle). \\
\end{matrix}
```
"""
const bell_kets = Ket.([
    1/sqrt(2)*[1,0,0,1],
    1/sqrt(2)*[1,0,0,-1],
    1/sqrt(2)*[0,1,1,0],
    1/sqrt(2)*[0,1,-1,0]
])

"""
    pure_qubit( ψ :: AbstractKet ) :: Qubit

A qubit is considered "pure" if it is rank-one. A rank-one density matrix
is constructed by  taking the outer-product of a ket state.

The method alternatively accepts a `Vector` input.

    pure_qubit( ψ :: Vector ) :: Qubit

A `DomainError` is thrown if `ψ` is not a valid ket.`
"""
function pure_qubit(ψ::AbstractKet)::Qubit
    Qubit(ψ*ψ')
end

function pure_qubit(ψ::Vector)::Qubit
    if !(is_ket(ψ))
        throw(DomainError(ψ, "ψ is not a valid wavefunction"))
    end

    Qubit(ψ*ψ')
end

"""
    mixed_state( priors :: QMath.Marginals, ρ_states :: Vector{<:AbstractDensityMatrix} ) :: DensityMatrix

Constructs the statistical mixture (weighted average) of quantum states. The
method accepts states as type `DensityMatrix` or subbtypes `AbstractDensityMatrix`.
"""
function mixed_state( priors::QMath.Marginals, ρ_set::Vector{<:AbstractDensityMatrix} )::DensityMatrix
    DensityMatrix(sum(priors .* ρ_set))
end

"""
    mixed_qubit( priors :: QMath.Marginals, ρ_states :: Vector{Qubit} ) :: Qubit

Constructs the statistical mixture (weighted average) of qubits.
"""
function mixed_qubit(priors::QMath.Marginals, ρ_set::Vector{Qubit})::Qubit
    Qubit(sum(priors .* ρ_set))
end

"""
    trine_qubit_kets :: Vector{QubitKet}

The triplet of kets representing three quantum states separated by equal angles in
the z-x plane of bloch sphere.

```jldoctest
julia> QBase.trine_qubit_kets == [[1.0, 0], [0.5, sqrt(3)/2], [0.5, -sqrt(3)/2]]
true
```
"""
const trine_qubit_kets = QubitKet.([
    [1.0, 0], [0.5, sqrt(3)/2], [0.5, -sqrt(3)/2]
])

"""
    mirror_symmetric_qubit_kets( θ :: Real ) :: Vector{QubitKet}

Returns the triplet of qubit kets in the z-x plane of bloch sphere. The first ket
is  ``|0\\rangle`` and the other two are symmetric across the z-axis.

Input:

* `θ ∈ [0,π/2]`: the hilbert space angle between |0> and symmetric kets.
"""
function mirror_symmetric_qubit_kets(θ::Real)::Vector{QubitKet}
    ψ1 = QubitKet([1;0])
    ψ2 = bloch_qubit_ket(2*θ,0)
    ψ3 = bloch_qubit_ket(2*θ,π)

    [ψ1,ψ2,ψ3]
end

"""
    trine_qubits :: Vector{Qubit}

Returns the qubit trine states in density matrix form.
"""
const trine_qubits = pure_qubit.(trine_qubit_kets)

"""
    mirror_symmetric_qubits(θ) :: Vector{Qubit}

Returns a set of 3 mirror symmetric qubit density matrices. The first
state is ``|0\\rangle\\langle 0|`` the other two are symmetric about the  ``|0\\rangle`` axis.

Input:
* `θ ∈ [0,π/2]`: the hilbert space angle between ``|0\\rangle`` and ``|\\psi_{2/3}\\rangle``.
"""
function mirror_symmetric_qubits(θ)::Vector{Qubit}
    pure_qubit.(mirror_symmetric_qubit_kets(θ))
end

"""
    bloch_qubit_ket( θ :: Real, ϕ :: Real ) :: QubitKet

Returns the qubit ket for the specified spherical coordinate on the surface of
bloch sphere, `(r=1, θ, ϕ)`.

Inputs:

* `θ ∈ [0, π]`: polar angle (w.r.t z-axis)
* `ϕ ∈ [0, 2π]`: azimuthal angle (x-y plane)

A `DomainError` is thrown if inputs `θ` and/or `ϕ` do are not within the valid range.
"""
function bloch_qubit_ket(θ::Real, ϕ::Real) :: QubitKet
    if !(0 <= θ <= π)
        throw(DomainError(θ, "polar bloch-angle (θ) must be in domain [0,π]"))
    elseif !(0 <= ϕ <= 2π)
        throw(DomainError(ϕ, "azimuthal angle (ϕ) must be in domain [0,2π]"))
    end

    QubitKet([cos(θ/2);exp(im*ϕ)*sin(θ/2)])
end

"""
Returns the qubit density matrix for quantum state described by a coordinate on
bloch sphere.

Spherical Coordinates:

States on the surface of bloch sphere may be described by spherical coordinates.

    bloch_qubit(θ::Real, ϕ::Real) :: Qubit

* `θ ∈ [0, π]`: polar angle (w.r.t z-axis).
* `ϕ ∈ [0, 2π]`: azimuthal angle (x-y plane)

Cartesian Coordinates:

States within the volume of bloch sphere may be described in cartesian coordinates.

    bloch_qubit(x::Real, y::Real, z::Real) :: Qubit

 * where `x`, `y`, and `z` are constrained to the unit sphere, `0 <= norm([x,y,z]) <= 1`.

 A `DomainError` is thrown if the coordinates `(x,y,z)` do  not adhere to constraints.
"""
function bloch_qubit(θ::Real,ϕ::Real) :: Qubit
    pure_qubit(bloch_qubit_ket(θ,ϕ))
end

function bloch_qubit(x,y,z) :: Qubit
    r = [x,y,z]

    if !(0 <= norm(r) <= 1)
        throw(DomainError(r, "norm([x,y,z]) is not in valid range"))
    end

    ρ = 0.5*([1 0;0 1] + sum(r .* paulis))

    Qubit(ρ)
end

"""
    sic_qubits :: Vector{Qubit}

The quadruplet of symmetric informationally complete (SIC) qubits. The qubits
are the vertices of a tetrahedron inscribed on bloch sphere.

```jldoctest
julia> QBase.sic_qubits
4-element Array{QBase.States.Qubit,1}:
 [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]
 [0.33333333333333337 + 0.0im 0.47140452079103173 + 0.0im; 0.47140452079103173 + 0.0im 0.6666666666666666 + 0.0im]
 [0.33333333333333337 + 0.0im -0.2357022603955158 - 0.4082482904638631im; -0.2357022603955158 + 0.4082482904638631im 0.6666666666666666 + 0.0im]
 [0.33333333333333337 - 0.0im -0.2357022603955158 + 0.4082482904638631im; -0.2357022603955158 - 0.4082482904638631im 0.6666666666666666 - 0.0im]
```
"""
const sic_qubits = [
    bloch_qubit([0,0,1]...),
    bloch_qubit([sqrt(2)*2,0,-1]./3...),
    bloch_qubit([-2/sqrt(2),sqrt(3)*sqrt(2),-1]./3...),
    bloch_qubit([-2/sqrt(2),-sqrt(3)*sqrt(2),-1]./3...)
]

"""
    bb84_qubits :: Vector{Qubit}

The quadruplet of qubits used in the BB84 Quantum Key Distribution protocol. The
states are ``|0\\rangle\\langle 0|``, ``|1\\rangle\\langle 1|``, ``|+\\rangle\\langle +|``, and ``|- \\rangle\\langle -|``.

```jldoctest
julia> QBase.bb84_qubits
4-element Array{QBase.States.Qubit,1}:
 [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]
 [0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]
 [0.5 + 0.0im 0.5 + 0.0im; 0.5 + 0.0im 0.5 + 0.0im]
 [0.5 + 0.0im -0.5 + 0.0im; -0.5 + 0.0im 0.5 + 0.0im]
```
"""
const bb84_qubits = Qubit.([
    [1 0;0 0], [0 0;0 1], [0.5 0.5;0.5 0.5], [0.5 -0.5;-0.5 0.5]
])

end
