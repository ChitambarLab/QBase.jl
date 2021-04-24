export State
export is_density_matrix, is_pure, is_mixed

"""
    is_density_matrix( ρ :: Matrix; atol=ATOL :: Float64 ) :: Bool

Returns `true` if `ρ` is a valid density matrix. The following constraints must be
satisfied for all density matrices:
* Hermitian: `ρ' == ρ`
* Positive Semidefinite: `eigmin(ρ) ≥ 0`
* Trace-one: `tr(ρ) == 1`
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
    State( ρ :: Matrix{T<:Number} ) <: Operator{T}

The density matrix representation of a quantum state. The constructor, `State(ρ)`
throws a `DomainError` if [`is_density_matrix`](@ref) evaluates to `false`.
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
    kron(states :: Vararg{State}; atol=ATOL :: Float64) :: State

Performs the kronecker product on the supplied `states` to construct a new `State`: ``\\rho_A \\otimes \\rho_B = \\rho_{AB}.``

This method extends [`Base.kron`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#Base.kron).
"""
kron(states :: Vararg{State}; atol=ATOL :: Float64) :: State = State(
    kron(map(ρ -> ρ.M, states)...),
atol=atol)

"""
    partial_trace(ρ::State, system::Vector{Int64}, id::Int64) :: State

Takes the partial_trace of a `State` `ρ` to construct a new `State`: ``\\text{Tr}_B[\\rho_{AB}] = \\rho_A``.
"""
partial_trace(ρ::State, system::Vector{Int64}, id::Int64; atol=ATOL :: Float64) :: State = begin
    State(partial_trace(ρ.M, system, id), atol=atol)
end

"""
    eigvals(ρ :: State) :: Vector{Float64}

Computes the eigenvalues of `State` `ρ`. The eigenvalues of a density matrix
are real. A `DomainError` is thrown if the imaginary component of any eigenvalues
exceeds the absolute tolerance `ρ.atol`. This method extends [`LinearAlgebra.eigvals`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.eigvals).
"""
function eigvals(ρ :: State) :: Vector{Float64}
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

Returns `true` if `ρ` is pure, i.e. `rank(ρ) == 1` within the absolute tolerance `ρ.atol`.
This method can also be used on `Matrix` types.

```julia
is_pure(ρ :: Matrix; atol=ATOL :: Float64) :: Bool
```
"""
function is_pure(ρ :: State) :: Bool
    rank(ρ) == 1
end
is_pure(ρ :: Matrix; atol=ATOL :: Float64) :: Bool = is_pure(State(ρ, atol=atol))

"""
    is_mixed(ρ :: State) :: Bool

Returns `true` if `ρ` is mixed, i.e. `rank(ρ) > 1` within the absolute tolerance `ρ.atol`.
This method also can be used on `Matrix` types.

```julia
is_mixed(ρ :: Matrix; atol=ATOL :: Float64) :: Bool
```
"""
function is_mixed(ρ :: State) :: Bool
    rank(ρ) > 1
end
is_mixed(ρ :: Matrix; atol=ATOL :: Float64) :: Bool = is_mixed(State(ρ, atol=atol))
