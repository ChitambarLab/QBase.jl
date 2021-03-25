export evolve

"""
    *( U :: Unitary, ket :: Ket ) :: Ket
    *( bra :: Bra, U :: Unitary ) :: Bra
    *()
"""
*(U :: Unitary, ket :: Ket) :: Ket = Ket(U.M * ket.ψ, atol=max(U.atol, ket.atol))
*(bra :: Bra, U :: Unitary) :: Bra = Bra(bra.ψ * U.M, atol=max(U.atol, bra.atol))
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
        ψ :: Ket
    ) :: Ket

Apply a unitary evolution `U` to a state ket `ψ`.
"""
evolve(U::Unitary, ψ::Ket) :: Ket = U*ψ
evolve(U :: AbstractMatrix, ψ :: AbstractVector; atol=ATOL :: Float64) :: Ket = Ket(U*ψ, atol=atol)
