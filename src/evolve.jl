export evolve

"""
    *( U :: Unitary, ket :: AbstractKet ) :: AbstractKet
    *( bra :: AbstractBra, U :: Unitary ) :: AbstractBra
    *()
"""
*(U :: Unitary, ket :: AbstractKet) :: AbstractKet = Ket(U.M * ket.ψ, atol=max(U.atol, ket.atol))
*(bra :: AbstractBra, U :: Unitary) :: AbstractBra = Bra(bra.ψ * U.M, atol=max(U.atol, bra.atol))
*(U1 :: Unitary, ρ :: State, U2 ::  Unitary) :: Matrix =  U1.M*ρ.M*U2.M

"""
Apply a unitary evolution `U` to density matrix `ρ`.

    evolve(
        U :: Unitary,
        ρ :: State
    ) :: DensityMatrix
"""
evolve(U :: Unitary, ρ :: State) :: State = State(
    U*ρ*U', atol=max(U.atol, ρ.atol)
)
evolve(U :: AbstractMatrix, ρ :: AbstractMatrix; atol=ATOL :: Float64) :: State = State(U*ρ*U', atol=atol)
"""
    evolve(
        U ::  Unitary,
        ψ :: AbstractKet
    ) :: AbstractKet

Apply a unitary evolution `U` to a state ket `ψ`.
"""
evolve(U::Unitary, ψ::AbstractKet) :: AbstractKet = U*ψ
evolve(U :: AbstractMatrix, ψ :: AbstractVector; atol=ATOL :: Float64) :: AbstractKet = Ket(U*ψ, atol=atol)
