"""
The `Channels` submodule provides:

* A catalog of functions which apply common quantum channels to quantum states.
"""
module Channels

using ..States
using LinearAlgebra

export depolarizing, erasure

"""
    replacer(
        ρ :: AbstractDensityMatrix,
        σ :: AbstractDensityMatrix,
        μ :: Real
    ) :: DensityMatrix


The replacer channel replaces the quantum state `ρ` with quantum state `σ` with
probability `μ`.
The replacer channel is expressed:

```math
    \\mathcal{R}_{\\mu}(\\rho) = \\mu \\rho + (1-\\mu)\\sigma
```

A `DomainError` is thrown if
* `μ` does not satisfy `1 ≥ μ ≥ 0`
* `ρ` and `σ` are not the same size
"""
function replacer(ρ :: AbstractDensityMatrix, σ :: AbstractDensityMatrix, μ :: Real) :: DensityMatrix
    if !(1 ≥ μ ≥ 0)
        throw(DomainError(μ, "Input `μ` must satisfy `1 ≥ μ ≥ 0`"))
    elseif size(ρ) != size(σ)
        throw(DomainError((ρ, σ), "Inputs `ρ` and `σ` must have the same size."))
    end

    DensityMatrix(μ*ρ + (1-μ)*σ)
end

@doc raw"""
    depolarizing( ρ :: AbstractDensityMatrix, μ :: Real ) :: DensityMatrix

The depolarizing channel mixes uniform classical noise into a quantum state `ρ`.
The argument `μ` describes the amount of noise mixed into the quantum states.
For a quantum state ``\rho``, the depolarizing channel is expressed:

```math
    \mathcal{D}_{\mu}(\rho) = \mu \rho + \frac{(1 - \mu)}{d}\mathbb{I}_d
```

A `DomainError` is thrown if `μ` does not satisfy `1 ≥ μ ≥ 0`.
"""
function depolarizing(ρ :: AbstractDensityMatrix, μ :: Real) :: DensityMatrix
    d = size(ρ, 1)
    σ = DensityMatrix(Matrix{Int64}(I, (d,d))/d)
    replacer(ρ, σ, μ)
end

@doc raw"""
    erasure( ρ :: AbstractDensityMatrix, μ :: Real ) :: DensityMatrix

The erasure channel mixes a quantum state `ρ` with an error flag ``|F\rangle``
orthogonal to the Hilbert space of `ρ`.
The argument `μ` describes the probability that `ρ` is replaced with the error flag.
For a quantum state ``\rho``, the erasure channel is expressed:

```math
    \mathcal{E}_{\mu}(\rho) = \mu \rho + (1 - \mu) |F \rangle \langle F|
```

Note that the erasure channel increases the dimension of the Hilbert space by 1.

A `DomainError` is thrown if `μ` does not satisfy `1 ≥ μ ≥ 0`.
"""
function erasure(ρ :: AbstractDensityMatrix, μ :: Real) :: DensityMatrix
    d = size(ρ, 1)

    ρ_ext = zeros(Complex{Float64}, d+1, d+1)
    ρ_ext[1:d,1:d] = ρ

    σ = zeros(Complex{Float64}, d+1, d+1)
    σ[d+1,d+1] = 1

    replacer(
        DensityMatrix(ρ_ext),
        DensityMatrix(σ),
        μ
    )
end

end
