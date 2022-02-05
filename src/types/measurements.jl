# validation
export is_povm, is_povm_element

# types
export Measurement
export POVMel, POVM, PVM

"""
    Measurement{T} <: AbstractVector{T}

The abstract type representing a quantum measurement. A measurement consists of a
complete set of possible outcomes.
"""
abstract type Measurement{T} <: AbstractVector{T} end
Base.size( Π :: Measurement) = size(Π.V)
Base.getindex( Π :: Measurement, id::Int) = getindex(Π.V, id...)
Base.setindex!( Π :: Measurement, val, id::Int) = (Π.V[id...] = val)

"""
    is_povm_element(M :: AbstractMatrix; atol=ATOL :: Float64) :: Bool

Returns `true` if matrix `M` satisfies the following constraints:
* `M` is Hermitian
* `M` is positive semi-definite
"""
function is_povm_element(M :: AbstractMatrix; atol=ATOL :: Float64) :: Bool
    if !is_hermitian(M, atol=atol)
        return false
    elseif !is_positive_semidefinite(M, atol=atol)
        return false
    end

    return true
end

"""
    is_povm( Π :: Vector; atol=ATOL :: Float64 ) :: Bool

Returns `true` if `Π` satisfies the following constraints
* The POVM is complete: `sum(Π) == I`
* Each POVM element is hermitian
* Each POVM element positive semi-definite
"""
function is_povm(Π::Vector{<:AbstractMatrix}; atol=ATOL :: Float64) :: Bool
    if !is_complete(Π, atol=atol)
        return false
    elseif !all(M -> is_povm_element(M, atol=atol), Π)
        return false
    end

    return true
end

"""
    POVMel( M :: AbstractMatrix{T<:Number}; atol=ATOL :: Float64) <: Operator{T}

A [`POVM`](@ref) element. A `DomainError` is thrown if matrix `M` is not hermitian
and positive semi-definite within absolute tolerance `atol`.
"""
struct POVMel{T<:Number} <: Operator{T}
    M :: Matrix{T}
    atol :: Float64
    POVMel(
        M :: Matrix{<:Number}; atol=ATOL :: Float64
    ) = is_povm_element(M,atol=atol) ? new{eltype(M)}(M,atol) : throw(DomainError(M, "POVM element M is invalid"))
    POVMel(
        ρ :: State; atol=ATOL :: Float64
    ) = is_povm_element(ρ, atol=atol) ? new{eltype(ρ.M)}(ρ.M,atol) : throw(DomainError(ρ.M, "POVM element M is invalid"))
end
is_povm_element(::POVMel) :: Bool = true
kron(Π_els :: Vararg{POVMel}) = POVMel(kron(map(Π_el -> Π_el.M, Π_els)...))

"""
    POVM( Π :: Vector{POVMel{T}} ) <: Measurement{T}
    POVM( Π :: Vector{Matrix{T}} ) <: Measurement{T}

Positve-operator valued measures (POVMs) represent a general quantum measurement.
Each POVM-element is a hermitian, positive-semidefinite matrix. The sum of all
POVM-elements yields the identity matrix. The constructor, `POVM(Π)` throws a
`DomainError` if the provided array of matrices, `Π` is not a valid POVM.
"""
struct POVM{T<:Number} <: Measurement{POVMel{T}}
    V :: Vector{POVMel{T}}
    atol :: Float64
    POVM(
        Π :: Vector{POVMel{T}},
        atol=ATOL :: Float64
    ) where T <: Number = is_povm(Π, atol=atol) ? new{T}(
            Π, atol
        ) : throw(DomainError(Π, "povm Π is invalid"))
    POVM(
        Π :: Vector{Matrix{T}};
        atol=ATOL :: Float64
    ) where T <: Number = is_povm(Π, atol=atol) ? new{T}(
            POVMel.(Π, atol=atol), atol
        ) : throw(DomainError(Π, "povm Π is invalid"))
    POVM(
        Π :: Vector{<:State};
        atol=ATOL :: Float64
    ) = is_povm(Π, atol=atol) ? new{eltype(Π[1].M)}(
            POVMel.(Π, atol=atol), atol
        ) : throw(DomainError(Π, "povm Π is invalid"))
end
is_povm(::POVM) :: Bool = true

# print out matrix forms when POVM types are displayed
function show(io::IO, mime::MIME{Symbol("text/plain")}, Π :: POVM)
    summary(io, Π)
    for i in 1:length(Π)
        print(io,"\nΠ[",i,"] : ")
        show(io, mime, Π[i])
    end
end

"""
    PVM( Π :: Vector{Vector{T}}; atol=ATOL :: Float64 ) <: Measurement{T}

The concret type for a projector-valued measure. The projectors are represented
as a set of orthonormal basis vectors

A `DomainError` is thrown if `Π` does not contain an orthonormal basis.
"""
struct PVM{T<:Number} <: Measurement{T}
    V :: Vector{Ket{T}}
    atol :: Float64
    PVM(
        Π :: Vector{Vector{T}};
        atol=ATOL :: Float64
    ) where T <: Number = is_orthonormal_basis(Π,  atol=atol) ? new{T}(Ket.(Π, atol=atol), atol) : throw(DomainError(Π, "povm Π is invalid"))
    PVM(
        Π :: Vector{Ket{T}};
        atol=ATOL :: Float64
    ) where  T <: Number = is_orthonormal_basis(Π,  atol=atol) ? new{T}(Π, atol) : throw(DomainError(Π, "povm Π is invalid"))
end
is_orthonormal_basis(::PVM) = true

# TODO: Observable
    # hermitian matrix M
    # M = ∑_j m_j * P_j
    # P_j are rank one projectors
    # m_j are real eigenvalues
    # P_j are the orthonormal eigenvectors of M
# TODO: POVM and PVM conversions
# TODO: kron extensions for measurements
