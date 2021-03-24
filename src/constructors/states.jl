# State Constructors
export pure_state, mixed_state
export bloch_qubit_state
export computational_basis_states
export bell_states, generalized_bell_states

# State ensemble constructors
export planar_symmetric_qubit_states, mirror_symmetric_qubit_states
export trine_qubit_states, bb84_qubit_states, sic_qubit_states

"""
    pure_state( ψ :: AbstractKet; atol=ATOL :: Float64 ) :: State

A state is considered "pure" if it is rank-one. A rank-one density matrix
is constructed by taking the outer-product of a `Ket` or `Bra`.
The method alternatively accepts a `Vector` input.

    pure_state( ψ :: Vector ) :: State
"""
function pure_state(ket::AbstractKet; atol=ATOL :: Float64) :: State
    State(ket.ψ*ket.ψ', atol=atol)
end
function pure_state(bra::AbstractBra; atol=ATOL :: Float64) :: State
    State(bra.ψ'*bra.ψ, atol=atol)
end
function pure_state(ψ::Vector{<:Number}; atol=ATOL :: Float64) :: State
    State(ψ*ψ', atol=atol)
end
function pure_state(
    ψ::Adjoint{T,Vector{T}} where T <: Number; atol=ATOL :: Float64
) :: State
    State(ψ'*ψ, atol=atol)
end

"""
    mixed_state(
        priors :: AbstractVector{<:Number},
        states :: Vector{<:AbstractMatrix{<:Number}}
    ) :: State

Constructs the statistical mixture (weighted average) of quantum states.
"""
function mixed_state(
    priors::AbstractVector{<:Number},
    states::Vector{<:AbstractMatrix{<:Number}};
    atol=ATOL :: Float64
) :: State
    if !is_probability_distribution(priors)
        throw(DomainError(priors, "Input `priors` is not a valid probability distribution."))
    elseif findfirst(!is_density_matrix, states) !== nothing
        throw(DomainError(states, "Input `states` is not a set of valid density matrices."))
    end

    State(sum(priors .* states))
end

"""
    computational_basis_states( dim :: Int64 ) :: Vector{State}

The density matrices for the computational basis of dimension, `dim`.
"""
function computational_basis_states(dim::Int64)::Vector{State}
    pure_state.(computational_basis_vectors(dim))
end

"""
    bell_states() :: Vector{State{Float64}}

The Bell basis density matrices. See [`bell_kets`](@ref) for more details.
"""
bell_states() :: Vector{State{Float64}} = pure_state.(bell_kets())

"""
    generalized_bell_states( dim :: Int64 ) :: Vector{State{Complex{Float64}}}

The density matrix representation of the generalized Bell basis. See  [`bell_kets`](@ref)
for more details.

A `DomainError` is thrown if `dim ≥ 2` is not satisfied.
"""
generalized_bell_states(dim :: Int64) :: Vector{State{Complex{Float64}}} = pure_state.(generalized_bell_kets(dim))

"""
    planar_symmetric_qubit_states( n :: Int64 ) :: Vector{State{Float64}}

Constructs a set of pure `State`s oriented symmetrically in the x-z-plane.
See [`planar_symmetric_qubit_kets`](@ref) for details.
"""
function planar_symmetric_qubit_states(n :: Int64) :: Vector{State{Float64}}
    pure_state.(planar_symmetric_qubit_kets(n))
end

"""
    mirror_symmetric_qubit_states( θ ::  Real ) :: Vector{State{Float64}}

Returns a set of 3 mirror symmetric qubit density matrices. The first
state is ``|0\\rangle\\langle 0|`` the other two are symmetric about the  ``|0\\rangle`` axis.

Input:
* `θ ∈ [0,π/2]`: the hilbert space angle between ``|0\\rangle`` and ``|\\psi_{2/3}\\rangle``.
"""
function mirror_symmetric_qubit_states(θ :: Real)::Vector{State{Float64}}
    pure_state.(mirror_symmetric_qubit_kets(θ))
end


"""
Returns the qubit density matrix for quantum state described by a coordinate on
bloch sphere.

Spherical Coordinates:

States on the surface of bloch sphere may be described by spherical coordinates.

    bloch_qubit_state(θ::Real, ϕ::Real) :: State{Complex{Float64}}

* `θ ∈ [0, π]`: polar angle (w.r.t z-axis).
* `ϕ ∈ [0, 2π]`: azimuthal angle (x-y plane)

Cartesian Coordinates:

States within the volume of bloch sphere may be described in cartesian coordinates.

    bloch_qubit_state(x::Real, y::Real, z::Real) :: State{Complex{Float64}}

 * where `x`, `y`, and `z` are constrained to the unit sphere, `0 <= norm([x,y,z]) <= 1`.

 A `DomainError` is thrown if the coordinates `(x,y,z)` do  not adhere to constraints.
"""
function bloch_qubit_state(θ::Real,ϕ::Real) :: State{Complex{Float64}}
    pure_state(bloch_qubit_ket(θ,ϕ))
end
function bloch_qubit_state(θ::Real) :: State{Float64}
    pure_state(bloch_qubit_ket(θ))
end
function bloch_qubit_state(x :: Real, y :: Real, z :: Real) :: State{Complex{Float64}}
    r = [x,y,z]

    if !(0 <= norm(r) <= 1)
        throw(DomainError(r, "norm([x,y,z]) is not in valid range"))
    end

    paulis_xyz = [[0 1;1 0],[0 -im;im 0],[1 0;0 -1]]

    ρ = 0.5*([1 0;0 1] + sum(r .* paulis_xyz))

    State(ρ)
end

"""
    trine_qubit_states() :: Vector{State{Float64}}

Returns the qubit trine states in density matrix form.
"""
trine_qubit_states() :: Vector{State{Float64}} = pure_state.(trine_qubit_kets())

"""
    sic_qubit_states :: Vector{State{Complex{Float64}}}

The quadruplet of symmetric informationally complete (SIC) qubits. The qubits
are the vertices of a tetrahedron inscribed on bloch sphere.

```jldoctest
julia> sic_qubits
4-element Array{QBase.State{Complex{Float64}},1}:
 [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]
 [0.33333333333333337 + 0.0im 0.47140452079103173 + 0.0im; 0.47140452079103173 + 0.0im 0.6666666666666666 + 0.0im]
 [0.33333333333333337 + 0.0im -0.2357022603955158 - 0.4082482904638631im; -0.2357022603955158 + 0.4082482904638631im 0.6666666666666666 + 0.0im]
 [0.33333333333333337 - 0.0im -0.2357022603955158 + 0.4082482904638631im; -0.2357022603955158 - 0.4082482904638631im 0.6666666666666666 - 0.0im]
```
"""
sic_qubit_states() :: Vector{State{Complex{Float64}}} = pure_state.(sic_qubit_kets())

"""
    bb84_qubit_states :: Vector{State}

The quadruplet of qubits used in the BB84 Quantum Key Distribution protocol. The
states are ``|0\\rangle\\langle 0|``, ``|1\\rangle\\langle 1|``, ``|+\\rangle\\langle +|``, and ``|- \\rangle\\langle -|``.

```jldoctest
julia> bb84_qubits
4-element Array{QBase.State{Float64},1}:
 [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]
 [0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]
 [0.5 + 0.0im 0.5 + 0.0im; 0.5 + 0.0im 0.5 + 0.0im]
 [0.5 + 0.0im -0.5 + 0.0im; -0.5 + 0.0im 0.5 + 0.0im]
```
"""
bb84_qubit_states() :: Vector{State{Float64}} = pure_state.(bb84_qubit_kets())
