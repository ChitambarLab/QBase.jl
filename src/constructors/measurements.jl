# constructors
export mirror_symmetric_qubit_3povm, asymmetric_qubit_3povm
export trine_qubit_povm, sic_qubit_povm
export planar_symmetric_qubit_povm

export sqrt_povm
export naimark_dilation

"""
    mirror_symmetric_qubit_3povm( θ :: Real ) :: POVM{Float64}

Constructs a `POVM` aligned with the three mirror symmetric qubit states.
The first measurement is aligned with the ``|0\\rangle`` state and the remaining
two are symmetric across the z-axis.

A `DomainError` is thrown if argument `θ` ∉ [π/4,π/2].
"""
function mirror_symmetric_qubit_3povm(θ::Real) :: POVM{Float64}
    if !( π/4 <= θ <= π/2)
        throw(DomainError(θ, "angle must be in range [π/4,π/2]"))
    end

    # POVM weights
    γ1 = 1 - cot(θ)^2
    γ2 = 0.5*(1 + cot(θ)^2)

    mirror_symmetric_qubits = mirror_symmetric_qubit_states(θ)

    π1 = γ1*mirror_symmetric_qubits[1]
    π2 = γ2*mirror_symmetric_qubits[2]
    π3 = γ2*mirror_symmetric_qubits[3]

    POVM([π1, π2, π3])
end

"""
    asymmetric_qubit_3povm( θ1::Real, θ2::Real ) :: POVM

Constructs a general non-orthogonal 3-element `POVM`.

Inputs:
* `θ1 ∈ (0,π/2] or [-π/2,0)`, the angle element 2 makes with ``|0\\rangle`` state.
* `θ2 ∈ [-π/2,0) or (0,π/2]`, the angle element 3 makes with ``|0\\rangle`` state.

A `DomainError` is thrown if `θ1` and `θ2` are not in the valid ranges. Note that
one angle must be positive and the other negative.
"""
function asymmetric_qubit_3povm(θ1::Real,θ2::Real) :: POVM{Float64}

    if !(((0 < θ1 <= π/2) && (-π/2 <= θ2 < 0)) || ((0 < θ2 <= π/2) && (-π/2 <= θ1 < 0)))
        throw(DomainError((θ1,θ2), "angles (θ1, θ2) are not in the valid range"))
    end

    # povm scaling terms derived analytically
    γ2 = 1/(sin(θ1)^2 - sin(θ1)*cos(θ1)*tan(θ2))
    γ3 = -γ2*(sin(θ1)*cos(θ1))/(sin(θ2)*cos(θ2))
    γ1 = 1 - γ2*cos(θ1)^2 - γ3*cos(θ2)^2

    π1 = γ1*[1 0;0 0]
    π2 = γ2*[cos(θ1)^2 cos(θ1)*sin(θ1); cos(θ1)*sin(θ1) sin(θ1)^2]
    π3 = γ3*[cos(θ2)^2 cos(θ2)*sin(θ2); cos(θ2)*sin(θ2) sin(θ2)^2]

    POVM([π1, π2, π3])
end

"""
    planar_symmetric_qubit_povm( n :: Int64 ) :: POVM{Float64}

Constructs an `n`-element `POVM{Float64}` from the [`planar_symmetric_qubits`](@ref).
Each state is multipled by a factor of `2/n` to satisfy the completeness relation.

A `DomainError` is thrown if `n ≥ 2` is not satisfied.
"""
function planar_symmetric_qubit_povm(n :: Int64) :: POVM{Float64}
    POVM(2/n*planar_symmetric_qubit_states(n))
end

"""
    trine_qubit_povm() :: POVM{Float64}

The `POVM` with elements parallel to the trine qubit states (see [`trine_qubit_kets`](@ref)).
"""
trine_qubit_povm() :: POVM{Float64} = POVM(2/3 * trine_qubit_states())

"""
    sic_qubit_povm() :: POVM{Complex{Float64}}

The `POVM` with elements parallel to the symmetric informationally complete (SIC) qubit states
(see [`sic_qubit_kets`](@ref)).
"""
sic_qubit_povm() :: POVM{Complex{Float64}} = POVM(0.5 * sic_qubit_states())

"""
    sqrt_povm(priors :: Probabilities, states :: Vector{<:State}) :: POVM

Returns the "pretty good" square-root povm for the given density operator
states and priors.
"""
function sqrt_povm(priors :: Probabilities, states :: Vector{<:State}) :: POVM
    ρ_mix = sum(priors .* states)
    ρ_sqrt = sqrt(inv(ρ_mix))

    Π = map(
        (ρ) -> ρ_sqrt * ρ * ρ_sqrt,
        (priors .* states)
    )

    POVM(Π)
end


"""
    _naimark_kraus_operators(Π::POVM) :: Array{Array{Complex{Float64},2},1}

Returns the Kraus operators for the provided POVM. In general, the Kraus operators
form a continuum and are non-unique. This method applies the construction

```math
k_i = \\sqrt{\\Pi_i}\\otimes |i\\rangle
```
"""
function _naimark_kraus_operators(Π::POVM) :: Array{Array{Complex{Float64},2},1}
    num_el = length(Π)

    map(i -> kron(sqrt(Π[i]), Matrix(I,num_el,num_el)[:,i]), 1:num_el)
end

"""
    naimark_dilation( Π::POVM )

Returns the dilated projectors which are equivalent to the provided POVM. During
measurement, the measured state must be tensored with the ancilla.
"""
function naimark_dilation(Π::POVM)
    n = length(Π)
    d = size(Π[1])[1]

    ancilla = zeros(n,n)
    ancilla[1,1] = 1

    k = _naimark_kraus_operators(Π)
    V = sum(k)

    orth_comp = nullspace(V')

    U_cols = []
    null_id = 1
    for i in 1:n*d
        if floor((i-1)/n) == (i-1)/n
            push!(U_cols,V[:,convert(Int,(i-1)/n + 1)])
        else
            push!(U_cols, orth_comp[:,null_id])
            null_id = null_id + 1
        end

    end

    U = hcat(U_cols...)

    p_ancilla = map( i -> begin
        p = zeros(n,n)
        p[i,i] = 1

        return p
    end, 1:n )

    projectors = map(i -> U'*kron(Matrix(1I,d,d),p_ancilla[i])*U, 1:n)

    Dict(
        "projectors" => projectors,
        "ancilla" => ancilla
    )
end
