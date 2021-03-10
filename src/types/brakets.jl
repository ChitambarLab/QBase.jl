export is_wave_vector
export AbstractKet, AbstractBra, Ket, Bra

"""
    AbstractKet{T<:Number} <: AbstractVector{Number}

The abstract type representing a quantum state ket. A ket can be thought of as a
column vector in a complex-valued Hilbert space. The elements of a `Ket` can be any
subtype of `Number`.
An `AbstractKet` cannot be instantiated, it serves as a supertype from
which ket types are defined.
"""
abstract type AbstractKet{T<:Number} <: AbstractVector{Number} end
Base.size(ket::AbstractKet) = size(ket.ψ)
Base.getindex(ket::AbstractKet, I::Int) = getindex(ket.ψ, I...)
Base.setindex!(ket::AbstractKet, val, I::Int) = (ket.ψ[I...] = val)

"""
    AbstractBra{T<:Number} <: AbstractVector{Number}

The abstract type representing a quantum state ket. A ket can be thought of as a
column vector in a complex-valued Hilbert space. The elements of a `Bra` can be any
subtype of `Number`.
An `AbstractBra` cannot be instantiated, it serves as a supertype from
which ket types are defined.
"""
abstract type AbstractBra{T<:Number} <: AbstractMatrix{Number} end
Base.size(bra::AbstractBra) = size(bra.ψ)
Base.getindex(bra::AbstractBra, I::Vararg{Int,2}) = getindex(bra.ψ, I...)
Base.setindex!(bra::AbstractBra, val, I::Vararg{Int,2}) = (bra.ψ[I...] = val)

"""
    is_wave_vector( ψ :: Vector; atol=ATOL :: Real) :: Bool

Returns `true` if vector `ψ` is a valid ket representation of a quantum state:

* `ψ` is a real or complex-valued vector.
* `ψ` is normalized with respect to the bra-ket inner prodcut (`ψ' * ψ == 0`).
"""
function is_wave_vector(ψ; atol=ATOL :: Real)
    isapprox( norm(ψ,2), 1, atol=atol )
end
function is_wave_vector(
    ψ :: Union{Matrix{<:Number},Adjoint{T,Matrix{T}}};
    atol=ATOL :: Real
) where T <: Number
    min(size(ψ)...) == 1 ? is_wave_vector(ψ[:], atol=atol) : false
end
is_wave_vector(ψ :: AbstractKet{<:Number}) = true
is_wave_vector(ψ :: AbstractBra{<:Number}) = true

"""
    _wave_vector_error(ψ)

Throws the `DomainError`, `"ψ must be a vector that satisfies \`norm(ψ,2) != 1\`."`.
"""
_wave_vector_error(ψ) = throw(DomainError(ψ, "ψ must be a vector that satisfies `norm(ψ,2) != 1`."))

"""
    Ket{T} <: AbstractKet{T} where T <: Number

A ket representation of a general quantum state. When given invalid input, the
constructor, `Ket(ψ)`, throws:
* `DomainError` - If `ψ` is not normalized.
"""
struct Ket{T} <: AbstractKet{T}
    ψ :: Vector{T}
    atol :: Float64
    Ket(
        ψ :: Vector{<:Number}; atol=ATOL :: Real
    ) = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(ψ, atol) : _wave_vector_error(ψ)
    Ket(
        ψ :: Matrix{<:Number}; atol=ATOL :: Real
    ) = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(ψ[:], atol) : _wave_vector_error(ψ)
    Ket(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Real
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
        ψ[:], atol
    ) : _wave_vector_error(ψ)
    Ket(
        ψ :: Adjoint{T,Matrix{T}}; atol=ATOL :: Real
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
        ψ[:], atol
    ) : _wave_vector_error(ψ)
    Ket(ψ :: AbstractBra{T}) where T <: Number = new{T}(ψ.ψ'[:], ψ.atol)
end

"""
    Bra{T} <: AbstractBra{T} where T <: Number

A bra representation of a general quantum state. When given invalid input, the
constructor, `Bra(ψ)`, throws:
* `DomainError` - If `ψ` is not normalized.
"""
struct Bra{T} <: AbstractBra{T}
    ψ :: Matrix{T}
    atol :: Float64
    Bra(
        ψ :: Matrix{<:Number}; atol=ATOL :: Real
    ) = is_wave_vector(ψ[:], atol=atol) ? new{eltype(ψ)}(
        ψ, atol
    ) : _wave_vector_error(ψ)
    Bra(
        ψ :: Adjoint{T, Matrix{T}}; atol=ATOL :: Real
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
        ψ, atol
    ) : _wave_vector_error(ψ)
    Bra(
        ψ :: Vector{<:Number}; atol=ATOL :: Real
    ) = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
            transpose!(zeros(eltype(ψ), 1, length(ψ)), ψ), atol
        ) : _wave_vector_error(ψ)
    Bra(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Real
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
            Matrix{T}(ψ), atol
        ) : _wave_vector_error(ψ)
    Bra(ψ :: AbstractKet{T}) where T <: Number = new{T}(ψ.ψ', ψ.atol)
end

"""
    *(ket :: AbstractKet{<:Number}, bra :: AbstractBra{<:Number}) :: Matrix{<:Number}
    *(bra :: AbstractBra{<:Number}, ket :: AbstractKet{<:Number}) :: Number

"""
*(ket :: AbstractKet{<:Number}, bra :: AbstractBra{<:Number}) :: Matrix{<:Number} = ket.ψ * bra.ψ
*(bra :: AbstractBra{<:Number}, ket :: AbstractKet{<:Number}) :: Number = (bra.ψ * ket.ψ)[1]

"""
    adjoint(ket :: AbstractKet{<:Number}) :: AbstractBra
    adjoint(bra :: AbstractBra{<:Number}) :: AbstractKet
"""
adjoint(ket :: AbstractKet{<:Number}) = Bra(ket)
adjoint(bra :: AbstractBra{<:Number}) = Ket(bra)

"""
    kron(kets :: Vararg{AbstractKet{<:Number}}; atol=ATOL) :: AbstractKet
    kron(bras :: Vararg{AbstractBra{<:Number}}; atol=ATOL) :: AbstractBra
"""
kron(kets :: Vararg{AbstractKet{<:Number}}; atol=ATOL) = Ket(kron(map(ket -> ket.ψ, kets)...), atol=atol)
kron(bras :: Vararg{AbstractBra{<:Number}}; atol=ATOL) = Bra(kron(map(bra -> bra.ψ, bras)...), atol=atol)
