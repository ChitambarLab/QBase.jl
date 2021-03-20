export AbstractUnitary, Unitary, is_unitary

"""
    AbstractUnitary <: AbstractMatrix{Complex{Float64}}

The abstract type representing unitary operators. An `AbstractUnitary` cannot be
instantiated, it serves as a supertype from which concrete types are derived.
"""
abstract type AbstractUnitary{T<:Number} <: AbstractMatrix{Number} end
Base.size(unitary::AbstractUnitary) = size(unitary.U)
Base.getindex(unitary::AbstractUnitary, id::Vararg{Int,2}) = getindex(unitary.U, id...)
Base.setindex!(unitary::AbstractUnitary, val, id::Vararg{Int,2}) = (unitary.U[id...] = val)

"""
    is_unitary( U :: Matrix; atol=ATOL :: Float64 ) :: Bool

Returns `true` if matrix `U` is unitary. The hermitian adjoint of a unitary matrix
is its inverse:
* `U' * U == I` where `I` is the identity matrix.
"""
function is_unitary(U::Matrix; atol=ATOL :: Float64) :: Bool
    if !isequal(size(U)...)
        return false
    end
    dim = size(U,1)

    isapprox(U'*U, Matrix(I,dim,dim), atol=atol)
end
is_unitary(U :: AbstractUnitary) = true

"""
    Unitary( U :: Matrix ) <: AbstractUnitary

Unitary matrices represent the physical evolution of quantum states. The Constructor,
`Unitary(U)`, throws a `DomainError` if the provided matrix, `U` is not unitary.
"""
struct Unitary{T} <: AbstractUnitary{T}
    U :: Matrix{T}
    atol :: Float64
    Unitary(
        U :: Matrix{T};
        atol=ATOL :: Float64
    ) where T <: Number = is_unitary(U, atol=atol) ? new{T}(U, atol) : throw(DomainError(U, "matrix U is not unitary"))
end

"""
*(unitaries :: Vararg{AbstractUnitary}; atol=ATOL :: Float64) :: AbstractUnitary
"""
*(unitaries :: Vararg{AbstractUnitary}; atol=ATOL :: Float64) = Unitary(*(map(unitary -> unitary.U, unitaries)...), atol=atol)

"""
kron(unitaries :: Vararg{AbstractUnitary}; atol=ATOL :: Float64) :: AbstractUnitary
"""
kron(unitaries :: Vararg{AbstractUnitary}; atol=ATOL :: Float64) = Unitary(kron(map(unitary -> unitary.U, unitaries)...), atol=atol)
