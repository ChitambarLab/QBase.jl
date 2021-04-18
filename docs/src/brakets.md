```@meta
CurrentModule = QBase
```
# Bras and Kets

In bra-ket notation, vectors on a complex-valued Hilbert space are used to represent
quantum states and operations.
A column vector ``|\psi\rangle`` is referred to as a **ket** whereas a row vector
``\langle \psi|`` is referred to as a **bra**.
Bras and kets are dual to each other via the adjoint operation.
That is, ``\langle\psi| = |\psi\rangle^{\dagger}`` and ``|\psi\rangle = \langle\psi|^{\dagger}``.

Quantum systems behave probabilistically under observation, hence quantum mechanics
is used to construct probability distributions which describe the behavior of quantum
systems.
For this reason, bras and kets must be normalized such that ``\langle\psi|\psi\rangle = 1`` holds true.
Here ``\langle  \cdot | \cdot \rangle`` denotes the inner product (dot product) between bra and ket.
These constraints are checked with the following method.

```@docs
is_braket
```

## Bra-Ket Types

```@docs
Ket
Bra
```

## Bra-Ket Algebra

```@docs
*(ket :: Ket{<:Number}, bra :: Bra{<:Number})
adjoint(ket :: Ket{<:Number})
kron(kets :: Vararg{Ket{<:Number}}; atol=ATOL)
```

## Ket State Constructors

QBase.jl provides a catalog for constructing various `Ket`s.
To construct similar `Bra`s you must manually convert the `Ket` to a `Bra` using
either [`adjoint(::Ket)`](@ref) or [`Bra(::Ket)`](@ref).

### Ket Singlets
```@docs
bloch_qubit_ket
```

### Ket Ensembles
```@docs
computational_basis_kets
bell_kets
generalized_bell_kets
mirror_symmetric_qubit_kets
planar_symmetric_qubit_kets
trine_qubit_kets
sic_qubit_kets
bb84_qubit_kets
```
