export is_wave_vector
export AbstractKet, AbstractBra, Ket, Bra

abstract type AbstractKet{T<:Number} <: AbstractVector{Number} end
Base.size(ket::AbstractKet) = size(ket.ψ)
Base.getindex(ket::AbstractKet, I::Int) = getindex(ket.ψ, I...)
Base.setindex!(ket::AbstractKet, val, I::Int) = (ket.ψ[I...] = val)

abstract type AbstractBra{T<:Number} <: AbstractMatrix{Number} end
Base.size(bra::AbstractBra) = size(bra.ψ)
Base.getindex(bra::AbstractBra, I::Vararg{Int,2}) = getindex(bra.ψ, I...)
Base.setindex!(bra::AbstractBra, val, I::Vararg{Int,2}) = (bra.ψ[I...] = val)

function is_wave_vector(ψ; atol=ATOL :: Float64)
    isapprox( norm(ψ,2), 1, atol=atol )
end
function is_wave_vector(
    ψ :: Union{Matrix{<:Number},Adjoint{T,Matrix{T}}};
    atol=ATOL :: Float64
) where T <: Number
    min(size(ψ)...) == 1 ? is_wave_vector(ψ[:], atol=atol) : false
end
is_wave_vector(ψ :: AbstractKet{<:Number}) = true
is_wave_vector(ψ :: AbstractBra{<:Number}) = true

_wave_vector_error(ψ) = throw(DomainError(ψ, "ψ must be a vector that satisfies `norm(ψ,2) != 1`."))

struct Ket{T} <: AbstractKet{T}
    ψ :: Vector{T}
    atol :: Float64
    Ket(
        ψ :: Vector{<:Number}; atol=ATOL :: Float64
    ) = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(ψ, atol) : _wave_vector_error(ψ)
    Ket(
        ψ :: Matrix{<:Number}; atol=ATOL :: Float64
    ) = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(ψ[:], atol) : _wave_vector_error(ψ)
    Ket(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
        ψ[:], atol
    ) : _wave_vector_error(ψ)
    Ket(
        ψ :: Adjoint{T,Matrix{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
        ψ[:], atol
    ) : _wave_vector_error(ψ)
    Ket(ψ :: AbstractBra{T}) where T <: Number = new{T}(ψ.ψ'[:], ψ.atol)
end

struct Bra{T} <: AbstractBra{T}
    ψ :: Matrix{T}
    atol :: Float64
    Bra(
        ψ :: Matrix{<:Number}; atol=ATOL :: Float64
    ) = is_wave_vector(ψ[:], atol=atol) ? new{eltype(ψ)}(
        ψ, atol
    ) : _wave_vector_error(ψ)
    Bra(
        ψ :: Adjoint{T, Matrix{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
        ψ, atol
    ) : _wave_vector_error(ψ)
    Bra(
        ψ :: Vector{<:Number}; atol=ATOL :: Float64
    ) = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
            transpose!(zeros(eltype(ψ), 1, length(ψ)), ψ), atol
        ) : _wave_vector_error(ψ)
    Bra(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
            Matrix{T}(ψ), atol
        ) : _wave_vector_error(ψ)
    Bra(ψ :: AbstractKet{T}) where T <: Number = new{T}(ψ.ψ', ψ.atol)
end

*(ket :: AbstractKet{<:Number}, bra :: AbstractBra{<:Number}) = ket.ψ * bra.ψ
*(bra :: AbstractBra{<:Number}, ket :: AbstractKet{<:Number}) :: Number = (bra.ψ * ket.ψ)[1]

adjoint(ket :: AbstractKet{<:Number}) = Bra(ket)
adjoint(bra :: AbstractBra{<:Number}) = Ket(bra)

kron(kets :: Vararg{AbstractKet{<:Number}}; atol=ATOL) = Ket(kron(map(ket -> ket.ψ, kets)...), atol=atol)
kron(bras :: Vararg{AbstractBra{<:Number}}; atol=ATOL) = Bra(kron(map(bra -> bra.ψ, bras)...), atol=atol)
