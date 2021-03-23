# validation
export is_povm

# types
export Measurement, POVM, PVM #AbstractPOVM, POVM

abstract type Measurement{T} <: AbstractVector{T} end
Base.size( m :: Measurement) = size(m.Π)
Base.getindex( m :: Measurement, id::Int) = getindex(m.Π, id...)
Base.setindex!( m :: Measurement, val, id::Int) = (m.Π[id...] = val)


"""
    AbstractPOVM <: AbstractMatrix{Complex{Float64}}

The abstract type representing positive-operator valued measures (POVMs). An
`AbstractPOVM` cannot be instantiated, it serves as a supertype from which concrete
types are derived.
"""
# abstract type AbstractPOVM{T<:Number} <: AbstractMeasurement end


"""
    is_povm( Π :: Vector; atol=ATOL :: Float64 ) :: Bool

Returns `true` if `Π` is a POVM. The following constraints must be satisfied:
* Each POVM element is hermitian
* Each POVM element positive semi-definite
* The POVM is complete: `sum(Π) == I`
"""
function is_povm(Π::Vector{Matrix}; atol=ATOL :: Float64) :: Bool
    dim = size(Π[1])[1]

    if !isapprox(sum(Π), Matrix{Complex{Float64}}(I,dim,dim), atol=atol)
        return false
    elseif !all(Πx -> is_hermitian(Πx, atol=atol), Π)
        return false
    elseif !all(Πx -> is_positive_semidefinite(Πx, atol=atol), Π)
        return false
    end

    return true
end

function is_pvm(Π::Vector{Vector}; atol=ATOL :: Float64) :: Bool

end

"""
    POVM( Π :: Vector{Matrix} ) <: Measurement

Positve-operator valued measures (POVMs) represent a general quantum measurement.
Each POVM-element is a hermitian, positive-semidefinite matrix. The sum of all
POVM-elements yields the identity matrix. The constructor, `POVM(Π)` throws a
`DomainError` if the provided array of matrices, `Π` is not a valid POVM.
"""
struct POVM{T} <: Measurement{T}
    Π :: Vector{Matrix{T}}
    atol :: Float64
    POVM(
        Π :: Vector{Matrix{T}};
        atol=ATOL :: Float64
    ) where T <: Number = is_povm(Π, atol=atol) ? new{T}(Π, atol) : throw(DomainError(Π, "povm Π is invalid"))
end

function show(io::IO, mime::MIME{Symbol("text/plain")}, povm :: POVM)
    summary(io, povm)
    print(io, "\n\natol : ")
    show(io,mime, povm.atol)
    println(io, "\nlength : ", length(povm))
    for i in 1:length(povm.Π)
        print(io,"\nΠ[",i,"] : ")
        show(io, mime, povm.Π[i])
    end
end
is_povm(Π::POVM) :: Bool = true

struct PVM{T} <: Measurement{T}
    Π :: Vector{Vector{T}}
    atol :: Float64
    PVM(
        Π :: Vector{Vector{T}};
        atol=ATOL :: Float64
    ) where T <: Number = true ? new{T}(Π, atol) : throw(DomainError(Π, "povm Π is invalid"))
end










#
# abstract type AbstractPVM{T<:Number} <: AbstractPOVM end
#
# abstract type AbstractObservable{T<:Number} <: AbstractMatrix{Number} end


# Observable
    # hermitian matrix M
    # M = ∑_j m_j * P_j
    # P_j are rank one projectors
    # m_j are real eigenvalues
    # P_j are the orthonormal eigenvectors of M

# PVM projection-valued measure
    # orthonormal and complete
    # should be represented by vector
    # measurements will be much faster that way

# POVM
    # Hermitian positive semi-definite elements that sum to the identity
    # vector of matrices

# Measurement Operator
    # matrix that is element of

# a PVM is also a POVM
    # can convert PVM into POVM

# PVM can be converted into an observable and vice-versa
