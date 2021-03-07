using LinearAlgebra

import Base: convert, *, kron, adjoint

export AbstractKet, AbstractBra
export Ket, Bra

const ATOL = 1e-7

function is_normalized(V :: Vector{<:Number}; atol=ATOL :: Float64)
    isapprox( norm(V,2), 1, atol=atol )
end
function is_normalized(V :: Adjoint{T, Vector{T}}; atol=ATOL :: Float64) where T <: Number
    isapprox( norm(V',2), 1, atol=atol )
end

abstract type AbstractKet{T<:Number} <: AbstractVector{Number} end
Base.size(K::AbstractKet) = size(K.ψ)
Base.getindex(K::AbstractKet, I::Int) = getindex(K.ψ, I...)
Base.setindex!(K::AbstractKet, v, I::Int) = (K.ψ[I...] = v)

abstract type AbstractBra{T<:Number} <: AbstractMatrix{Number} end
Base.size(B::AbstractBra{<:Number}) = size(B.ψ)
Base.getindex(B::AbstractBra{<:Number}, I::Vararg{Int,2}) = getindex(B.ψ, I...)
Base.setindex!(B::AbstractBra{<:Number}, v, I::Vararg{Int,2}) = (B.ψ[I...] = v)

struct Ket{T} <: AbstractKet{T}
    ψ :: Vector{T}
    atol :: Float64
    Ket(
        ψ :: Vector{<:Number};
        atol=ATOL :: Float64
    ) = begin
        println("uouofduofdu")
        is_normalized(ψ, atol=atol) ? new{eltype(ψ)}(ψ, atol) : throw(DomainError(ψ, "ψ is not normalized."))
    end
    Ket(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_normalized(ψ, atol=atol) ? new{eltype(ψ)}(
            ψ[:], atol
        ) : throw(DomainError(ψ, "ψ is not normalized."))
    Ket(ψ :: AbstractBra{T}) where T <: Number = new{T}(ψ.ψ'[:], ψ.atol)
end

struct Bra{T} <: AbstractBra{T}
    ψ :: Matrix{T}
    atol :: Float64
    Bra(
        ψ :: Matrix{<:Number}; atol=ATOL :: Float64
    ) = is_normalized(ψ[:], atol=atol) ? new{eltype(ψ)}(
        ψ, atol
    ) : throw(DomainError(ψ, "ψ is not normalized."))
    Bra(
        ψ :: Vector{<:Number}; atol=ATOL :: Float64
    ) = is_normalized(ψ, atol=atol) ? new{eltype(ψ)}(
            transpose!(zeros(T, 1, length(ψ)), ψ), atol
        ) : throw(DomainError(ψ, "ψ is not normalized."))
    Bra(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_normalized(ψ, atol=atol) ? new{eltype(ψ)}(
            Matrix{T}(ψ), atol
        ) : throw(DomainError(ψ, "ψ is not normalized."))
    Bra(ψ :: AbstractKet{<:Number}) = new{T}(ψ.ψ', ψ.atol)
end

*(K :: AbstractKet{<:Number}, B :: AbstractBra{<:Number}) = K.ψ * B.ψ
*(B :: AbstractBra{<:Number}, K :: AbstractKet{<:Number}) = B.ψ * K.ψ
adjoint(K :: AbstractKet{<:Number}) = Bra(K)
adjoint(B :: AbstractBra{<:Number}) = Ket(B)
kron(kets :: Vararg{AbstractKet{<:Number}}) = Ket(kron(map(ket -> ket.ψ, kets)...))
kron(bras :: Vararg{AbstractBra{<:Number}}) = Bra(kron(map(bra -> bra.ψ, bras)...))


# rig adjoint operator to make bra

# AbstractOperator < Matrix
    # AbstractUnitary
    # AbstractKrausOperator

# AbstractChannel

# AbstractSuperOperator

# AbstractObservable

# AbstractMeasurement
