export measure

"""
Computes the outcome probabilities for a quantum measurement. The conditional
probabilities are determined by the Born rule, ``P(i|j) = \\text{Tr}[\\Pi_i \\rho_j]``,
where ``\\Pi_j`` is a POVM element and ``\\rho_j`` is a density matrix.

Measurement of a single `Ket` or  `DensityMatrix`:

    measure(
        Π :: POVM,
        ρ :: State
    ) :: Conditionals

    measure(
        Π :: POVM,
        ψ :: AbstractKet
    ) :: Conditionals

Measurement of an ensemble of `Ket` or `DensityMatrix` types:

    measure(
        Π :: POVM,
        ρ_states :: Vector{<:State}
    ) :: Conditionals

    measure(
        Π :: POVM,
        ψ_kets :: Vector{<:ABstractKet}
    ) :: Conditionals
"""
function measure(
    Π :: POVM,
    ρ :: State
) :: Probabilities
    Probabilities( real.(map(Π_el -> tr(Π_el * ρ), Π)) )
end
function measure(
    Π :: POVM,
    states :: Vector{<:State}
) :: Conditionals
    Conditionals( real.(tr.(transpose(states * transpose(Π))) ))
end
function measure(
    Π :: POVM,
    ψ :: AbstractKet
) :: Probabilities
    Probabilities( real.(map( Π_el -> ψ' * Π_el * ψ, Π)) )
end
function measure(
    Π :: POVM,
    ψ_kets :: Vector{<:AbstractKet}
) :: Conditionals
    measure(Π, pure_state.(ψ_kets))
end
