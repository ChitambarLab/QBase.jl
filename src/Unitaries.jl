"""
Quantum states evolve under unitary transformations. The `QBase.Unitaries`
submodule provides:
* Types and Constructors for unitary operators.
"""
module Unitaries

using ..QMath
using LinearAlgebra

# validation
export is_unitary

# types
export Unitary, QubitUnitary

# constants
export paulis, σx, σy, σz

# contructors
export qubit_rotation

"""
    AbstractUnitary <: AbstractMatrix{Complex{Float64}}

The abstract type representing unitary operators. An `AbstractUnitary` cannot be
instantiated, it serves as a supertype from which concrete types are derived.
"""
abstract type AbstractUnitary <: AbstractMatrix{Complex{Float64}} end
Base.size(unitary::AbstractUnitary) = size(unitary.U)
Base.getindex(unitary::AbstractUnitary, id::Vararg{Int,2}) = getindex(unitary.U, id...)
Base.setindex!(unitary::AbstractUnitary, v, id::Int) = (unitary.U[id...] = v)

"""
    is_unitary( U :: Matrix ) :: Bool

Returns `true` if matrix `U` is unitary. The hermitian adjoint of a unitary matrix
is its inverse:
* `U'U == I` where `I` is the identity matrix.

A unitary matrix must be square. A `DomainError` is thrown if input `U` is not square.
"""
function is_unitary(U::Matrix) :: Bool
    if !(QMath.is_square_matrix(U))
        throw(DomainError(U, "provided matrix U is not square"))
    end

    dim = size(U)[1]

    U'*U ≈ Matrix{Complex{Float64}}(I,dim,dim)
end

"""
    Unitary( U :: Matrix ) <: AbstractUnitary

Unitary matrices represent the physical evolution of quantum states. The Constructor,
`Unitary(U)`, throws a `DomainError` if the provided matrix, `U` is not unitary.
"""
struct Unitary <: AbstractUnitary
    U :: Matrix{Complex{Float64}}
    Unitary(U) = is_unitary(U) ? new(U) : throw(DomainError(U, "matrix U is not unitary"))
end

"""
    QubitUnitary( U :: Matrix ) <: AbstractUnitary

Constructs a 2x2 unitary for qubit evolution. Throws a `DomainError` if input
`U` is not of dimension 2x2 or if `U` is not unitary.
"""
struct QubitUnitary <: AbstractUnitary
    U :: Unitary
    QubitUnitary(U) = (size(U) == (2,2)) ? new(Unitary(U)) : throw(DomainError(U, "matrix U is not 2x2"))
end

"""
    σx :: QubitUnitary

Pauli-X unitary:

```jldoctest
julia> QBase.σx
2×2 QBase.Unitaries.QubitUnitary:
 0.0+0.0im  1.0+0.0im
 1.0+0.0im  0.0+0.0im
```
"""
const σx = QubitUnitary([0 1;1 0])

"""
    σy :: QubitUnitary

Pauli-Y unitary:

```jldoctest
julia> QBase.σy
2×2 QBase.Unitaries.QubitUnitary:
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im
```
"""
const σy = QubitUnitary([0 -im;im 0])

"""
    σz :: QubitUnitary

Pauli-Z unitary:

```
julia> QBase.σz
2×2 QBase.Unitaries.QubitUnitary:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im
```
"""
const σz = QubitUnitary([1 0;0 -1])

"""
    paulis :: Vector{QubitUnitary}

Returns a vector containing the three qubit pauli matrices, `[σx, σy, σz]`.
"""
const paulis = [σx, σy, σz]

"""
    qubit_rotation( θ :: Real; axis="x" :: String ) :: QubitUnitary

Returns a unitary which performs a qubit rotation along bloch sphere. `θ` designates
the angle of rotation and `axis` (`"x"`, `"y"`, `"z"`) designates the cartesian
axis about which the qubit is rotated.
"""
function qubit_rotation(θ :: Real; axis="x" :: String)::QubitUnitary

    pauli = [1 0; 0 1]
    if axis == "x"
        pauli = σx
    elseif axis == "y"
        pauli = σy
    elseif axis == "z"
        pauli = σz
    end

    QubitUnitary(exp(-im * pauli * θ/2))
end

end
