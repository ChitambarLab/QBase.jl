export Unitary, is_unitary

"""
    is_unitary( U :: Matrix; atol=ATOL :: Float64 ) :: Bool

Returns `true` if matrix `U` is unitary. The hermitian adjoint of a unitary matrix
is its inverse:
* `U' * U == I` where `I` is the identity matrix.
"""
function is_unitary(U::AbstractMatrix{<:Number}; atol=ATOL :: Float64) :: Bool
    if !isequal(size(U)...)
        return false
    end
    dim = size(U,1)

    isapprox(U'*U, Matrix(I,dim,dim), atol=atol)
end

"""
    Unitary( U :: AbstractMatrix ) <: Operator{T}

Unitary matrices represent the physical evolution of quantum states. The Constructor,
`Unitary(U)`, throws a `DomainError` if the provided matrix, `U` is not unitary.
"""
struct Unitary{T<:Number} <: Operator{T}
    M :: Matrix{T}
    atol :: Float64
    Unitary(
        U :: AbstractMatrix{<:Number};
        atol=ATOL :: Float64
    ) = is_unitary(U, atol=atol) ? new{eltype(U)}(U, atol) : throw(DomainError(U, "matrix U is not unitary"))
end
is_unitary(U :: Unitary) = true

"""
*(unitaries :: Vararg{Unitary}; atol=ATOL :: Float64) :: Unitary
"""
*(unitaries :: Vararg{Unitary}; atol=ATOL :: Float64) :: Unitary = Unitary(*(map(U -> U.M, unitaries)...), atol=atol)

"""
kron(unitaries :: Vararg{Unitary}; atol=ATOL :: Float64) :: Unitary
"""
kron(unitaries :: Vararg{Unitary}; atol=ATOL :: Float64) :: Unitary = Unitary(kron(map(U -> U.M, unitaries)...), atol=atol)

adjoint(U :: Unitary) :: Unitary = Unitary(U.M', atol=U.atol)
