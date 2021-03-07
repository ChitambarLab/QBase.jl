export is_wave_vector, wave_vector_error
export AbstractKet, AbstractBra, Ket, Bra

abstract type AbstractKet{T<:Number} <: AbstractVector{Number} end
Base.size(ket::AbstractKet) = size(ket.ψ)
Base.getindex(ket::AbstractKet, I::Int) = getindex(ket.ψ, I...)
Base.setindex!(ket::AbstractKet, val, I::Int) = (ket.ψ[I...] = val)

abstract type AbstractBra{T<:Number} <: AbstractMatrix{Number} end
Base.size(bra::AbstractBra{<:Number}) = size(bra.ψ)
Base.getindex(bra::AbstractBra{<:Number}, I::Vararg{Int,2}) = getindex(bra.ψ, I...)
Base.setindex!(bra::AbstractBra{<:Number}, val, I::Vararg{Int,2}) = (bra.ψ[I...] = val)

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

wave_vector_error(ψ) = throw(DomainError(ψ, "ψ is not normalized `norm(ψ,2) != 1`."))

struct Ket{T} <: AbstractKet{T}
    ψ :: Vector{T}
    atol :: Float64
    Ket(
        ψ :: Vector{<:Number}; atol=ATOL :: Float64
    ) = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(ψ, atol) : wave_vector_error(ψ)
    Ket(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
            ψ[:], atol
        ) : wave_vector_error(ψ)
    Ket(ψ :: AbstractBra{T}) where T <: Number = new{T}(ψ.ψ'[:], ψ.atol)
end

struct Bra{T} <: AbstractBra{T}
    ψ :: Matrix{T}
    atol :: Float64
    Bra(
        ψ :: Matrix{<:Number}; atol=ATOL :: Float64
    ) = is_wave_vector(ψ[:], atol=atol) ? new{eltype(ψ)}(
        ψ, atol
    ) : wave_vector_error(ψ)
    Bra(
        ψ :: Vector{<:Number}; atol=ATOL :: Float64
    ) = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
            transpose!(zeros(T, 1, length(ψ)), ψ), atol
        ) : wave_vector_error(ψ)
    Bra(
        ψ :: Adjoint{T,Vector{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_wave_vector(ψ, atol=atol) ? new{eltype(ψ)}(
            Matrix{T}(ψ), atol
        ) : wave_vector_error(ψ)
    Bra(ψ :: AbstractKet{<:Number}) = new{T}(ψ.ψ', ψ.atol)
end

*(ket :: AbstractKet{<:Number}, bra :: AbstractBra{<:Number}) = ket.ψ * bra.ψ
*(bra :: AbstractBra{<:Number}, ket :: AbstractKet{<:Number}) = bra.ψ * ket.ψ

adjoint(ket :: AbstractKet{<:Number}) = Bra(ket, atol=ket.atol)
adjoint(bra :: AbstractBra{<:Number}) = Ket(bra, atol=bra.atol)

kron(kets :: Vararg{AbstractKet{<:Number}}; atol=Atol) = Ket(kron(map(ket -> ket.ψ, kets)...), atol=atol)
kron(bras :: Vararg{AbstractBra{<:Number}}; atol=Atol) = Bra(kron(map(bra -> bra.ψ, bras)...), atol=atol)

# abstract type AbstractState{T<:Number} <: AbstractMatrix{Number} end
# Base.size(state::AbstractState{<:Number}) = size(state.ρ)
# Base.getindex(state::AbstractState{<:Number}, I::Vararg{Int,2}) = getindex(state.ρ, I...)
# Base.setindex!(state::AbstractState{<:Number}, val, I::Vararg{Int,2}) = (state.ρ[I...] = val)
#
# struct State{T} <: AbstractState{T}
#     ρ :: Matrix{T}
#     atol :: Float64
#     State(
#         ρ :: Matrix{<:Number}; atol=ATOL :: Float64
#     ) = is_density_matrix(ρ, atol=atol) ? new{eltype(ρ)}(ρ, atol) : density_matrix_error(ρ)
# end
#
# kron(states :: Vararg{AbstractState{<:Number}};atol=ATOL) = State(
#     kron(map(state -> state.ρ, states)...),
# atol=atol)
#
# function is_density_matrix(ρ::Matrix; atol=ATOL)::Bool
#     is_hermitian = QMath.is_hermitian(ρ) # TODO: atol
#     is_pos_sd = QMath.is_positive_semidefinite(ρ) # TODO: atol
#     is_trace_one = isapprox(tr(ρ), 1, atol=atol)
#
#     (is_hermitian & is_pos_sd & is_trace_one)
# end


# AbstractOperator < Matrix
    # AbstractUnitary
    # AbstractKrausOperator

# AbstractChannel

# AbstractSuperOperator

# AbstractObservable

# AbstractMeasurement
