export measure

"""
Computes the outcome probabilities for a quantum measurement. The conditional
probabilities are determined by the Born rule, ``P(i|j) = \\text{Tr}[\\Pi_i \\rho_j]``,
where ``\\Pi_j`` is a POVM element and ``\\rho_j`` is a density matrix.

Measurement of a single `Ket` or  `State`:

    measure(Π :: POVM, ρ :: State) :: Probabilities
    measure(Π :: POVM, ψ :: Ket) :: Probabilities
    measure(Π :: PVM, ρ :: State) :: Probabilities
    measure(Π :: PVM, ψ :: Ket) :: Probabilities


Measurement of an ensemble of `Ket` or `State` types:

    measure(Π :: POVM, ρ_states :: Vector{<:State}) :: Conditionals
    measure(Π :: POVM, ψ_kets :: Vector{<:Ket}) :: Conditionals
    measure(Π :: PVM, ρ_states :: Vector{<:State}) :: Conditionals
    measure(Π :: PVM, ψ_kets :: Vector{<:Ket}) :: Conditionals
"""
function measure(
    Π :: POVM,
    ρ :: State
) :: Probabilities
    Probabilities( real.(map(Π_el -> tr(Π_el * ρ), Π)) )
end
function measure(
    Π :: PVM,
    ρ :: State
) :: Probabilities
    Probabilities( real.(map(Π_el -> Π_el' * ρ * Π_el, Π)))
end
function measure(
    Π :: POVM,
    states :: Vector{<:State}
) :: Conditionals
    Conditionals( real.(tr.(transpose(states * transpose(Π))) ))
end
function measure(
    Π :: PVM,
    states :: Vector{<:State}
) :: Conditionals
    Conditionals(
        real.(hcat(map( ρ -> map(Π_el -> Π_el' * ρ * Π_el, Π), states)...))
    )
end
function measure(
    Π :: POVM,
    ψ :: Ket
) :: Probabilities
    Probabilities( real.(map( Π_el -> ψ' * Π_el * ψ, Π)) )
end
function measure(
    Π :: PVM,
    ψ :: Ket
) :: Probabilities
    Probabilities( map( Π_el -> norm(Π_el' * ψ)^2, Π) )
end
function measure(
    Π :: POVM,
    ψ_kets :: Vector{<:Ket}
) :: Conditionals
    measure(Π, pure_state.(ψ_kets))
end
function measure(
    Π :: PVM,
    ψ_kets :: Vector{<:Ket}
) :: Conditionals
    Conditionals( hcat(map( ψ -> map(Π_el -> norm(Π_el' * ψ)^2, Π), ψ_kets)...) )
end
