```@meta
CurrentModule = QBase.States
```

# QBase.States

```@docs
States
```

In discrete and finite quantum systems, states can be represented computationally
with vectors (Bra-Ket Representation) or matrices (Density Matrix Representation).

## Bra-Ket Representation

A quantum state may be represented by a vector on a complex-valued Hilbert space.
In Bra-Ket notation, a quantum state is represented by a ket denoted ``|\psi\rangle``.
A ket can simply be understood as a column vector. Each ket has an associated dual
vector known as a bra. A bra is may be represented by a row vector and can be
constructed from a ket via Hermitian adjoint, ``\langle\psi| = (|\psi\rangle)^{\dagger}``.

It is essential that a quantum states are normalized, ``\langle\psi|\psi\rangle = 1``,
where ``\langle \; | \; \rangle`` denotes the inner product (dot product) between bra and ket.

The `QBase.States` module provides a validation method, `is_ket()` for checking
whether a vector satisfies the requirements for being a quantum state ket.

```@docs
is_ket
```

### Ket Types

```@docs
AbstractKet
```

The `QBase.States` module provides two concrete subtypes of `AbstractKet`:

```@docs
Ket
QubitKet
```

### Ket Constructors

```@docs
basis_kets
```

#### Singlet States
```@docs
bloch_qubit_ket
```


#### Triplet States
```@docs
trine_qubit_kets
mirror_symmetric_qubit_kets
```

## Density Matrix Representation

Quantum states can be represented in matrix form.

```@docs
is_density_matrix
AbstractDensityMatrix
DensityMatrix
Qubit
```

## Density Matrix Constructors

The density matrix ``\rho`` can be constructed from ket ``|\psi\rangle`` by taking
the outer-product of the ket with itself, ``|\psi\rangle\langle\psi| = \rho``.

```@docs
pure_state
pure_qubit
basis_states
```

The rank of the density matrix can be greater than 1. If a density matrix has
rank greater than one, we call the state mixed.

```@docs
mixed_state
mixed_qubit
bloch_qubit
```

### Triplets
```@docs
trine_qubits
mirror_symmetric_qubits
```

### Quadruplets
```@docs
sic_qubits
bb84_qubits
```
