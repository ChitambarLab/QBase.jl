export Operator

"""
    abstract type Operator{T<:Number} <: AbstractMatrix{Number} end

Operators are matrix objects which describe the dynamics of a quantum system.
"""
abstract type Operator{T<:Number} <: AbstractMatrix{Number} end
Base.size(O :: Operator) = size(O.M)
Base.getindex(O :: Operator, id :: Vararg{Int,2}) = getindex(O.M, id...)
Base.setindex!(O :: Operator, val, id :: Vararg{Int,2}) = (O.M[id...] = val)

*(operators :: Vararg{Operator}) :: Matrix = *(map(O -> O.M, operators)...)
kron(operators :: Vararg{Operator}) :: Matrix =  kron(map(O -> O.M, operators)...)
rank(O :: Operator) :: Int64 = rank(O.M, atol=O.atol)
sqrt(O :: Operator) :: Matrix = sqrt(O.M)

*(O :: Operator, ket :: Ket) :: Vector = O.M * ket.ψ
*(bra :: Bra, O :: Operator) :: Adjoint{T, Vector{T}} where T <: Number = bra.ψ * O.M
*(bra :: Bra, O :: Operator, ket :: Ket) :: Number = bra.ψ * O.M * ket.ψ

*(O :: Operator, M :: Matrix{<:Number}) :: Matrix = O.M * M
*(M :: Matrix{<:Number}, O :: Operator) :: Matrix = M * O.M
