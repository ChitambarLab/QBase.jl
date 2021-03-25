export State
export is_density_matrix, is_pure, is_mixed

"""
    is_density_matrix( ρ :: Matrix{<:Number}; atol=ATOL :: Float64 ) :: Bool

Returns true if input `ρ` is:
* Hermitian
* Positive Semi-Definite
* Trace[ρ] = 1 (normalization)
"""
function is_density_matrix(ρ::Matrix{<:Number}; atol=ATOL :: Float64)::Bool
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
    State( ρ :: Matrix{T<:Number} ) <: Operator{T}

The density matrix representation of a quantum state. The constructor, `DensityMatrix(ρ)`
throws a `DomainError` if `is_density_matrix(ρ)` is `false`.

Base methods extended to use the `DensityMatrix` type:
* `partial_trace` - Returns `DensityMatrix` if supplied with one.
* `kron` - The kronecker product of two quantum states is a `State`.
"""
struct State{T<:Number} <: Operator{T}
    M :: Matrix{T}
    atol :: Float64
    State(
        ρ :: Matrix{<:Number}; atol=ATOL :: Float64
    ) = is_density_matrix(ρ, atol=atol) ? new{eltype(ρ)}(ρ, atol) : _density_matrix_error(ρ, atol=atol)
    State(
        ρ :: Adjoint{T,Matrix{T}}; atol=ATOL :: Float64
    ) where T <: Number = is_density_matrix(ρ[:,:], atol=atol) ? new{eltype(ρ)}(ρ[:,:], atol) : _density_matrix_error(ρ[:,:], atol=atol)
end
is_density_matrix(::State) = true


"""
    kron(states :: Vararg{State}; atol=ATOL :: Float64)
"""
kron(states :: Vararg{State}; atol=ATOL :: Float64) = State(
    kron(map(ρ -> ρ.M, states)...),
atol=atol)

"""
    partial_trace(ρ::State, system::Vector{Int64}, id::Int64)
"""
partial_trace(ρ::State, system::Vector{Int64}, id::Int64; atol=ATOL :: Float64) = begin
    State(partial_trace(ρ.M, system, id), atol=atol)
end

"""
    eigvals(state :: State)
"""
function eigvals(ρ :: State)
    λs = eigvals(ρ.M)
    for λ in λs
        if !isapprox(imag(λ), 0, atol=ρ.atol)
            throw(DomainError(λ), "eigenvalue λ has a imaginary component outside the absolute tolerance.")
        end
    end

    real.(λs)
end

"""
    is_pure(ρ :: State) :: Bool

Returns `true` if `rho` is pure, i.e. `rank(ρ) == 1`. This method also can be used
on an matrix.
"""
function is_pure(ρ :: State) :: Bool
    rank(ρ) == 1
end
is_pure(ρ :: Matrix{<:Number}; atol=ATOL :: Float64) :: Bool = is_pure(State(ρ, atol=atol))

"""
    is_mixed(ρ :: State) :: Bool

Returns `true` if `rho` is mixed, i.e. `rank(ρ) > 1`. This method also can be used
on an matrix.
"""
function is_mixed(ρ :: State) :: Bool
    rank(ρ) > 1
end
is_mixed(ρ :: Matrix{<:Number}; atol=ATOL :: Float64) :: Bool = is_mixed(State(ρ, atol=atol))
