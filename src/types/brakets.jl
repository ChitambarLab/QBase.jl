export is_wave_vector
export Ket, Bra

"""
    is_wave_vector( ψ :: Vector; atol=ATOL :: Float64) :: Bool

Returns `true` if vector `ψ` is a valid ket representation of a quantum state:

* `ψ` is a real or complex-valued vector.
* `ψ` is normalized with respect to the bra-ket inner prodcut (`ψ' * ψ == 0`).
"""
function is_wave_vector(ψ :: Vector{<:Number}; atol=ATOL :: Float64)
    isapprox( norm(ψ,2), 1, atol=atol )
end
function is_wave_vector(ψ :: Adjoint{T, Vector{T}}; atol=ATOL :: Float64) where T <: Number
    isapprox( norm(ψ,2), 1, atol=atol )
end

"""
    _wave_vector_error(ψ)

Throws the `DomainError`, `"ψ must be a vector that satisfies \`norm(ψ,2) != 1\`."`.
"""
_wave_vector_error(ψ) = throw(DomainError(ψ, "ψ must be a vector that satisfies `norm(ψ,2) != 1`."))

"""
    Ket(ψ :: Vector{T}; atol=ATOL :: Float64) where T <: Number

A normalized column vector on a complex-valued Hilbert space. If called with
an adjoint vector `ψ'`, the adjoint will be taken on construction.

    Ket(ψ :: Adjoint{Vector{T}}; atol=ATOL :: Float64) where T <: Number

Similarly a `Bra` can be converted into a `Ket` via the adjoint.

    Ket(ψ :: Bra{T}; atol=ATOL :: Float64) where T <: Number

When given invalid input, the constructor, `Ket(ψ)`, throws:
* `DomainError` - If `ψ` is not normalized.
"""
struct Ket{T<:Number} <: AbstractVector{Number}
    ψ :: Vector{T}
    atol :: Float64
    Ket(
        ψ :: Vector{T}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{T}(ψ, atol) : _wave_vector_error(ψ)
    Ket(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{T}(
        ψ', atol
    ) : _wave_vector_error(ψ)
end
Base.size(ket::Ket) = size(ket.ψ)
Base.getindex(ket::Ket, I::Int) = getindex(ket.ψ, I...)
Base.setindex!(ket::Ket, val, I::Int) = (ket.ψ[I...] = val)
is_wave_vector(::Ket) = true

"""
    Bra(ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Float64) where T <: Number

A row vector on a complex valued HIlbet space. Bras are dual to Kets via the adjoint
If called with an `Vector{T}` type, the adjoint will automatically be taken `ψ'`.

    Bra(ψ :: Vector{T}; atol=ATOL :: Float64) where T <: Number

Similarly a `Ket` can be converted into a `Bra` via the adjoint.

    Bra(ψ :: Ket{T}; atol=ATOL :: Float64) where T <: Number

When given invalid input, the constructor, `Ket(ψ)`, throws:
* `DomainError` - If `ψ` is not normalized.
"""
struct Bra{T<:Number} <: AbstractMatrix{Number}
    ψ :: Adjoint{T,Vector{T}}
    atol :: Float64
    Bra(
        ψ :: Vector{T}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{T}(ψ', atol) : _wave_vector_error(ψ)
    Bra(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{T}(ψ, atol) : _wave_vector_error(ψ)
end
Base.size(bra::Bra) = size(bra.ψ)
Base.getindex(bra::Bra, I::Vararg{Int,2}) = getindex(bra.ψ, I...)
Base.setindex!(bra::Bra, val, I::Vararg{Int,2}) = (bra.ψ[I...] = val)
is_wave_vector(::Bra) = true

Ket(bra :: Bra) = Ket(bra.ψ', atol=bra.atol)
Bra(ket :: Ket) = Bra(ket.ψ', atol=ket.atol)


"""
    *(ket :: Ket{<:Number}, bra :: Bra{<:Number}) :: Matrix{<:Number}
    *(bra :: Bra{<:Number}, ket :: Ket{<:Number}) :: Number

"""
*(ket :: Ket{<:Number}, bra :: Bra{<:Number}) :: Matrix{<:Number} = ket.ψ * bra.ψ
*(bra :: Bra{<:Number}, ket :: Ket{<:Number}) :: Number = (bra.ψ * ket.ψ)[1]

"""
    adjoint(ket :: Ket{<:Number}) :: Bra
    adjoint(bra :: Bra{<:Number}) :: Ket
"""
adjoint(ket :: Ket{<:Number}) = Bra(ket)
adjoint(bra :: Bra{<:Number}) = Ket(bra)

"""
    kron(kets :: Vararg{Ket{<:Number}}; atol=ATOL) :: Ket
    kron(bras :: Vararg{Bra{<:Number}}; atol=ATOL) :: Bra
"""
kron(kets :: Vararg{Ket{<:Number}}; atol=ATOL) = Ket(kron(map(ket -> ket.ψ, kets)...), atol=atol)
kron(bras :: Vararg{Bra{<:Number}}; atol=ATOL) = Bra(kron(map(bra -> bra.ψ', bras)...), atol=atol)
