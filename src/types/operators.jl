export Operator

"""
    abstract type Operator{T} <: AbstractMatrix{T} end

Operators are matrix objects which describe the dynamics of a quantum system.
"""
abstract type Operator{T} <: AbstractMatrix{T} end
Base.size(O :: Operator) = size(O.M)
Base.getindex(O :: Operator, id :: Vararg{Int,2}) = getindex(O.M, id...)
Base.setindex!(O :: Operator, val, id :: Vararg{Int,2}) = (O.M[id...] = val)

*(operators :: Vararg{Operator}) :: Matrix = *(map(O -> O.M, operators)...)
kron(operators :: Vararg{Operator}) :: Matrix =  kron(map(O -> O.M, operators)...)
rank(O :: Operator) = rank(O.M, atol=O.atol)
sqrt(O :: Operator) :: Matrix = sqrt(O.M)


*(O :: Operator, ket :: AbstractKet) :: Vector = O.M * ket.ψ
*(bra :: AbstractBra, O :: Operator) :: Adjoint{T, Vector{T}} where T <: Number = bra.ψ * O.M
*(bra :: AbstractBra, O :: Operator, ket :: AbstractKet) :: Number = bra.ψ * O.M * ket.ψ
