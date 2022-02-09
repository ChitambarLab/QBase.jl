export is_kraus_channel, KrausChannel, kraus_evolve

"""
    is_kraus_channel(
        kraus_ops :: Vector{<:AbstractMatrix};
        atol=ATOL :: Float64
    ) :: Bool

Returns `true` if the `kraus_ops` constitute a quantum channel.
The requirements are:
* Completeness, ``\\sum_{i} K^{\\dagger}_i K_i = \\mathbb{I}``.
"""
function is_kraus_channel(
    kraus_ops :: Vector{<:AbstractMatrix};
    atol=ATOL :: Float64
) :: Bool
    if !isapprox(sum(K -> K'K, kraus_ops), I, atol=atol)
        @warn "`kraus_ops` do not sum to identity."
        return false
    end

    return true
end

"""
    KrausChannel(
        kraus_ops :: Vector{Matrix{T}};
        atol=ATOL :: Float64
    ) :: KrausChannel{T}

The Kraus operator representation of a quantum channel.
The channel evolves a density operator ``\\rho`` as,

```math
    \\mathcal{N}(\\rho) = \\sum_{i} K_{i} \\rho K_{i}^{\\dagger},
```

where ``\\{K_i \\in L(A,B)\\}_i`` are a set of linear operators known as
Kraus operators.

A `DomainError` is thrown if `kraus_ops` do not constitute a quantum channel.
See [`is_kraus_channel`](@ref) for details.
"""
struct KrausChannel{T<:Number}
    kraus_ops :: Vector{Matrix{T}}
    atol :: Float64
    KrausChannel(
        kraus_ops :: Vector{Matrix{T}};
        atol=ATOL :: Float64
    ) where T <: Number = is_kraus_channel(
        kraus_ops, atol=atol
    ) ? new{T}(kraus_ops, atol) : throw(
        DomainError(kraus_ops, "`kraus_ops` do not represent a channel.")
    )
end

# print out matrix forms when Channel types are displayed
function show(io::IO, mime::MIME{Symbol("text/plain")}, 𝒩 :: KrausChannel)
    summary(io, 𝒩)
    print(io, "\natol : ", 𝒩.atol)
    print(io, "\nkraus_ops : ")
    show(io, mime, 𝒩.kraus_ops)
end

"""
    kraus_evolve(
        kraus_ops :: Vector{<:AbstractMatrix},
        ρ :: Matrix
    ) :: Matrix

Applies the Kraus operators `kraus_ops` to evolve the density operator `\rho`.
The evolved state ``\\rho'`` is constructed as,

```math
    \\rho' = \\sum_{i} K_{i} \\rho K_{i}^{\\dagger}
```

where ``K_i`` is a Kraus operator.
"""
kraus_evolve(
    kraus_ops :: Vector{<:AbstractMatrix}, ρ :: Matrix
) :: Matrix = sum(K -> K * ρ * K', kraus_ops)
