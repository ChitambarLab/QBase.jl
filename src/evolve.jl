export evolve

"""
    *(U :: Unitary, ket :: Ket) :: Ket
    *(bra :: Bra, U :: Unitary) :: Bra
    *(U1 :: Unitary, ρ :: State, U2 ::  Unitary) :: Matrix
"""
*(U :: Unitary, ket :: Ket) :: Ket = Ket(U.M * ket.ψ, atol=max(U.atol, ket.atol))
*(bra :: Bra, U :: Unitary) :: Bra = Bra(bra.ψ * U.M, atol=max(U.atol, bra.atol))
*(U1 :: Unitary, ρ :: State, U2 ::  Unitary) :: Matrix =  U1.M*ρ.M*U2.M

"""
Apply a unitary evolution `U` to density matrix `ρ`: ``\\rho' = U\\rho U^{\\dagger}``.

    evolve(U :: Unitary, ρ :: State) :: State
"""
evolve(U :: Unitary, ρ :: State) :: State = State(
    U*ρ*U', atol=max(U.atol, ρ.atol)
)
evolve(U :: AbstractMatrix, ρ :: AbstractMatrix; atol=ATOL :: Float64) :: State = State(U*ρ*U', atol=atol)
"""
Apply a unitary evolution `U` to a state ket `ψ`: ``|\\psi'\\rangle = U|\\psi\\rangle``.

    evolve(U ::  Unitary, ψ :: Ket) :: Ket
"""
evolve(U::Unitary, ψ::Ket) :: Ket = U*ψ
evolve(U :: AbstractMatrix, ψ :: AbstractVector; atol=ATOL :: Float64) :: Ket = Ket(U*ψ, atol=atol)

"""
Evolve a `State` `ρ` by `Λ` the [`ChoiOp`](@ref) representation of a channel.

    evolve(Λ ::  ChoiOp, ρ :: State) :: State
    evolve(Λ :: ChoiOp, ρ :: AbstractMatrix) :: Matrix

See the [`choi_evolve`](@ref) method for details.
"""
evolve(Λ::ChoiOp, ρ::State) :: State = State(choi_evolve(Λ.M, ρ.M, Λ.dims))
evolve(Λ::ChoiOp, ρ::AbstractMatrix) :: Matrix = choi_evolve(Λ.M, ρ, Λ.dims)

"""
Evolve a `State` `ρ` by the [`KrausChannel`](@ref) `𝒩`.

    evolve(𝒩 :: KrausChannel, ρ :: State) :: State
    evolve(𝒩 :: KrausChannel, ρ :: AbstractMatrix) :: Matrix

See the [`kraus_evolve`](@ref) method for details.
"""
evolve(𝒩::KrausChannel, ρ::State) :: State = State(kraus_evolve(𝒩.kraus_ops, ρ.M))
evolve(𝒩::KrausChannel, ρ::AbstractMatrix) :: Matrix = kraus_evolve(𝒩.kraus_ops, ρ)
