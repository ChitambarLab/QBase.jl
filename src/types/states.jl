export AbstractState, State
export is_density_matrix

"""
    AbstractState{T<:Number} <: AbstractMatrix{Number}

The abstract type representing all density matrices.
"""
abstract type AbstractState{T<:Number} <: AbstractMatrix{Number} end
Base.size(state::AbstractState{<:Number}) = size(state.ρ)
Base.getindex(state::AbstractState{<:Number}, I::Vararg{Int,2}) = getindex(state.ρ, I...)
Base.setindex!(state::AbstractState{<:Number}, val, I::Vararg{Int,2}) = (state.ρ[I...] = val)

"""
    is_density_matrix( ρ :: Matrix; atol=ATOL :: Float64 ) :: Bool

Returns true if input `ρ` is:
* Hermitian
* Positive Semi-Definite
* Trace[ρ] = 1 (normalization)
"""
function is_density_matrix(ρ::Matrix; atol=ATOL :: Float64)::Bool
    if !is_hermitian(ρ, atol=atol)
        return false
    elseif !isapprox(tr(ρ), 1, atol=atol)
        return false
    elseif !is_positive_semidefinite(ρ, atol=atol)
        return false
    end

    return true
end

"""
    _density_matrix_error(ρ; atol=ATOL :: Float64)

Throws a `DomainError` indicating why density matrix `ρ` is not a density matrix
the error messages are self-explanatory and are listed as follows:
* `"Density matrix `ρ` is not hermitian."`
* `"Density matrix `ρ` is not trace-one."`
* `"Density matrix `ρ` is not positive semi-definite."`
"""
function _density_matrix_error(ρ; atol=ATOL :: Float64)
    if !is_hermitian(ρ, atol=atol)
        throw(DomainError(ρ, "Density matrix `ρ` is not hermitian."))
    elseif !isapprox(tr(ρ), 1, atol=atol)
        throw(DomainError(ρ, "Density matrix `ρ` is not trace-one."))
    elseif !is_positive_semidefinite(ρ, atol=atol)
        throw(DomainError(ρ, "Density matrix `ρ` is not positive semi-definite."))
    end
end

"""
    State( ρ :: Matrix{Complex{Float64}} ) <: AbstractDensityMatrix

The density matrix representation of a quantum state. The constructor, `DensityMatrix(ρ)`
throws a `DomainError` if `is_density_matrix(ρ)` is `false`.

Base methods extended to use the `DensityMatrix` type:
* `partial_trace` - Returns `DensityMatrix` if supplied with one.
* `kron` - The kronecker product of two quantum states is a `State`.
"""
struct State{T} <: AbstractState{T}
    ρ :: Matrix{T}
    atol :: Float64
    State(
        ρ :: Matrix{<:Number}; atol=ATOL :: Float64
    ) = is_density_matrix(ρ, atol=atol) ? new{eltype(ρ)}(ρ, atol) : _density_matrix_error(ρ, atol=atol)
    State(
        ρ :: Adjoint{T,Matrix{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_density_matrix(ρ[:,:], atol=atol) ? new{eltype(ρ)}(ρ[:,:], atol) : _density_matrix_error(ρ[:,:], atol=atol)
end

"""
    kron(states :: Vararg{AbstractState}; atol=ATOL :: Float64)
"""
kron(states :: Vararg{AbstractState}; atol=ATOL :: Float64) = State(
    kron(map(state -> state.ρ, states)...),
atol=atol)

"""
    partial_trace(ρ::AbstractState, system::Vector{Int64}, id::Int64)
"""
partial_trace(ρ::AbstractState, system::Vector{Int64}, id::Int64; atol=ATOL) = begin
    State(partial_trace(ρ.ρ, system, id), atol=atol)
end
