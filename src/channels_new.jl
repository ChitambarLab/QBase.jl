export replacer_channel, depolarizing_channel, erasure_channel

"""
    replacer_channel(
        ρ :: AbstractState,
        σ :: AbstractState,
        μ :: Real
    ) :: AbstractState


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
function replacer_channel(ρ :: AbstractState, σ :: AbstractState, μ :: Real) :: AbstractState
    if !(1 ≥ μ ≥ 0)
        throw(DomainError(μ, "Input `μ` must satisfy `1 ≥ μ ≥ 0`"))
    elseif size(ρ) != size(σ)
        throw(DomainError((ρ, σ), "Inputs `ρ` and `σ` must have the same size."))
    end

    State(μ*ρ + (1-μ)*σ)
end

"""
    depolarizing_channel( ρ :: AbstractState, μ :: Real ) :: AbstractState

The depolarizing channel mixes uniform classical noise into a quantum state `ρ`.
The argument `μ` describes the amount of noise mixed into the quantum states.
For a quantum state ``\\rho``, the depolarizing channel is expressed:

```math
    \\mathcal{D}_{\\mu}(\\rho) = \\mu \\rho + \\frac{(1 - \\mu)}{d}\\mathbb{I}_d
```

A `DomainError` is thrown if `μ` does not satisfy `1 ≥ μ ≥ 0`.
"""
function depolarizing_channel(ρ :: AbstractState, μ :: Real) :: AbstractState
    d = size(ρ, 1)
    σ = State(Matrix{Int64}(I, (d,d))/d)
    replacer_channel(ρ, σ, μ)
end

"""
    erasure_channel( ρ :: AbstractState, μ :: Real ) :: State

The erasure channel mixes a quantum state `ρ` with an error flag ``|F\\rangle``
orthogonal to the Hilbert space of `ρ`.
The argument `μ` describes the probability that `ρ` is replaced with the error flag.
For a quantum state ``\\rho``, the erasure channel is expressed:

```math
    \\mathcal{E}_{\\mu}(\\rho) = \\mu \\rho + (1 - \\mu) |F \\rangle \\langle F|
```

Note that the erasure channel increases the dimension of the Hilbert space by 1.

A `DomainError` is thrown if `μ` does not satisfy `1 ≥ μ ≥ 0`.
"""
function erasure_channel(ρ :: AbstractState, μ :: Real) :: AbstractState
    d = size(ρ, 1)

    ρ_ext = zeros(Complex{Float64}, d+1, d+1)
    ρ_ext[1:d,1:d] = ρ

    σ = zeros(Complex{Float64}, d+1, d+1)
    σ[d+1,d+1] = 1

    replacer_channel(
        State(ρ_ext),
        State(σ),
        μ
    )
end
