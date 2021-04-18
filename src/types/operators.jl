export Operator

"""
    abstract type Operator{T<:Number} <: AbstractMatrix{Number} end

A matrix representing a linear operator that acts upon a complex-valued Hilbert space.
The `Operator` is an abstract type that parents all linear operator in quantum mechanics.
These matrix types can represent quantum states, evolution, and measurements, each with
their individual constraint.
These constraints are place upon the children of th `Operator` abstract type.
"""
abstract type Operator{T<:Number} <: AbstractMatrix{Number} end
Base.size(O :: Operator) = size(O.M)
Base.getindex(O :: Operator, id :: Vararg{Int,2}) = getindex(O.M, id...)
Base.setindex!(O :: Operator, val, id :: Vararg{Int,2}) = (O.M[id...] = val)

"""
Operator Multiplication:

```julia
*(operators :: Vararg{Operator}) :: Matrix
```
"""
*(operators :: Vararg{Operator}) :: Matrix = *(map(O -> O.M, operators)...)

"""
Kronecker Product:

```julia
kron(operators :: Vararg{Operator}) :: Matrix
```
"""
kron(operators :: Vararg{Operator}) :: Matrix =  kron(map(O -> O.M, operators)...)

"""
Matrix Rank:

```julia
rank(O :: Operator) :: Int64
```
"""
rank(O :: Operator) :: Int64 = rank(O.M, atol=O.atol)

"""
Matrix Square Root:

```julia
sqrt(O :: Operator) :: Matrix
```
"""
sqrt(O :: Operator) :: Matrix = sqrt(O.M)

"""
[`Operator`](@ref) types can multiply [`Bra`](@ref) and [`Ket`](@ref) types.

``O|\\psi\\rangle = | \\psi' \\rangle``:

```julia
*(O :: Operator, ket :: Ket) :: Vector
```

``\\langle \\psi |O = \\langle \\psi'|``:

```julia
*(bra :: Bra, O :: Operator) :: Adjoint{T, Vector{T}} where T <: Number
```

Inner product, ``\\langle \\psi |O|\\sigma\\rangle = \\langle \\psi|\\sigma'\\rangle``:

```julia
*(bra :: Bra, O :: Operator, ket :: Ket) :: Number
```
"""
*(O :: Operator, ket :: Ket) :: Vector = O.M * ket.ψ
*(bra :: Bra, O :: Operator) :: Adjoint{T, Vector{T}} where T <: Number = bra.ψ * O.M
*(bra :: Bra, O :: Operator, ket :: Ket) :: Number = bra.ψ * O.M * ket.ψ

*(O :: Operator, M :: Matrix{<:Number}) :: Matrix = O.M * M
*(M :: Matrix{<:Number}, O :: Operator) :: Matrix = M * O.M
