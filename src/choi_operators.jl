export is_choi_matrix, choi_matrix, ChoiOp, choi_evolve

"""
    is_choi_matrix(
        Λ :: AbstractMatrix,
        dims :: Vector{Int};
        atol=ATOL :: Float
    ) :: Bool

Returns `true` if matrix `Λ` is a Choi operator representation of a quantum
channel.

The requirements on `Λ` are ([https://arxiv.org/abs/1111.6950](https://arxiv.org/abs/1111.6950)):
* Is Hermitian, ``\\Lambda = \\Lambda^{\\dagger}``.
* Is trace-preserving, ``\\text{Tr}_A[\\Lambda^{AB}]=\\mathbb{I}_B``.
* Is completely-positive, ``\\Lambda \\geq 0``.
"""
function is_choi_matrix(
    Λ :: AbstractMatrix,
    dims :: Vector{Int};
    atol=ATOL::Float64
) :: Bool
    if !is_hermitian(Λ, atol=atol)
        @warn "The Choi matrix `Λ` is not Hermitian-preserving.`"
        return false
    elseif !is_positive_semidefinite(Λ, atol=atol)
        @warn "The Choi matrix `Λ` is not completely-positive."
        return false
    elseif !isapprox(partial_trace(Λ, dims, 1), I, atol=atol)
        @warn "The Choi matrix `Λ` is not trace-preserving."
        return false
    end

    return true
end

"""
    choi_matrix(𝒩 :: Function, dims :: Vector{Int}) :: Matrix{ComplexF64}

Returns the Choi operator of a channel. The Choi matrix is constructed as

```math
        \\Lambda^{AB} = \\sum_{i,j \\in [d_A]} E^A_{i,j} \\otimes \\mathcal{N}^B(E_{i,j}) ,
```

where ``[d_A]`` is the finite alphabet indexing the input space and ``E_{i,j}``
is a square matrix of dimension ``d_A`` with a ``1`` in the ``(i,j)`` entry
and a ``0`` everywhere else. The input ``\\Lambda`` is the output dimension.

The input function `𝒩` is called as `𝒩(X)` for arbitrary input `X`.
Channel functions with multiple parameters can be considered by declaring
a function
`𝒩_xy(ρ) = 𝒩(ρ,x,y)` for fixed `(x,y)` and then call, `choi(𝒩_xy, dim_in, dim_out)`.
"""
function choi_matrix(𝒩 :: Function, dims :: Vector{Int}) :: Matrix{ComplexF64}
    dim_in, dim_out = dims

    eab_matrix = zeros(ComplexF64, dim_in , dim_in)
    Λ = zeros(ComplexF64, dim_in * dim_out, dim_in * dim_out)
    for i in 1:dim_in
        row_ids = ((i-1) * dim_out + 1):((i-1) * dim_out + dim_out)
        for j in 1:dim_out
            col_ids = ((j-1) * dim_out + 1):((j-1) * dim_out + dim_out)

            eab_matrix[i,j] = 1
            Λ[row_ids,col_ids] += 𝒩(eab_matrix)
            eab_matrix[i,j] = 0
        end
    end

    return Λ
end

"""
    ChoiOp( Λ :: AbstractMatrix, dims :: Vector{Int} ) :: ChoiOp{<:Number}
    ChoiOp( 𝒩 :: Function, dims :: Vector{Int} ) :: ChoiOp{ComplexF64}

Constructs the Choi operator representation of a quantum channel.
If either a function `𝒩` or set of kraus operators is provided as input, the
Choi matrix is constructed with the [`choi_matrix`](@ref) method.

The `ChoiOp` type contains the fields:
* `M :: Matrix{<:Number}` - The Choi matrix.
* `dims :: Vector{Int}` - The channel's input/output dimension `[dim_in, dim_out]`.

A `DomainError` is thrown if [`is_choi_matrix`](@ref) returns `false`.
"""
struct ChoiOp{T<:Number} <: Operator{T}
    M :: Matrix{T}
    dims :: Vector{Int}
    atol :: Float64
    ChoiOp(
        Λ :: AbstractMatrix{T},
        dims :: Vector{Int};
        atol=ATOL :: Float64
    ) where T <: Number = is_choi_matrix(Λ, dims, atol=atol) ? new{T}(Λ, dims, atol) : throw(
        DomainError(Λ, "The Choi operator is not a valid quantum channel.")
    )
    ChoiOp(
        𝒩 :: Function,
        dims :: Vector{Int};
        atol=ATOL :: Float64
    ) = ChoiOp( choi_matrix(𝒩, dims), dims, atol=atol)
end

# print out matrix forms when Choi types are displayed
function show(io::IO, mime::MIME{Symbol("text/plain")}, Λ :: ChoiOp)
    summary(io, Λ)
    print(io, "\ndims : ", Λ.dims)
    print(io, "\natol : ", Λ.atol)
    print(io, "\nM : ")
    show(io, mime, Λ.M)
end

"""
    choi_evolve(Λ :: Matrix{<:Number}, ρ :: Matrix{<:Number}) :: Matrix
    choi_evolve(
        Λ :: Matrix{<:Number}, ρ :: Matrix{<:Number}, dims :: Vector{Int}
    ) :: Matrix

Applies the Choi operator `Λ` to the density operator `ρ`. The output of this
quantum channel evaluated as,

```math
    \\rho'_B = \\text{Tr}_A[\\Lambda_{AB}(\\rho_A^{T}\\otimes \\mathbb{I}_B],
```

where the [`partial_trace`](@ref) is take with respect to the input system.
"""
function choi_evolve(Λ :: Matrix{<:Number}, ρ :: Matrix{<:Number}) :: Matrix
    dim_in = size(ρ, 1)
    dim_out = size(Λ, 1) ÷ dim_in

    return choi_evolve(Λ, ρ, [dim_in, dim_out])
end
function choi_evolve(
    Λ :: Matrix{<:Number}, ρ :: Matrix{<:Number}, dims :: Vector{Int}
) :: Matrix
    ρI = kron(transpose(ρ), Matrix{Int}(I, dims[2], dims[2]))
    return partial_trace(Λ * ρI, dims, 1)
end
