```@meta
CurrentModule = QBase
```
# Base Types

A goal of QBase.jl is to provide a general framework for representing finite quantum
systems and their dynamics.
Quantum mechanics is simply an application of linear algebra where particular
constraints are applied to the vectors and matrices involved in the representation
of quantum systems.
Hence, there are two core data structures used to  represent quantum systems
1. **Bras and Kets:**  Row and column vectors defined on a complex-valued Hilbert space.
2. **Operators:** Matrices defined on a complex-valued Hilbert space.

Furthermore, quantum systems behave probabilistically under observation. Thus,
QBase.jl also provides data structures for the various probability distributions.

Details regarding the definitions and constraints for each of these data structure
are provided below.
Ideally, the constraints on quantum objects should be met exactly, however, numerical
errors are innate in the computations involved.
Therefore, each type has an absolute tolerance parameter `atol` which specifies how
much error is allowed before the quantum object is deemed invalid.
By default, the `atol=1e-7` and is stored in the constant `QBase.ATOL`.
This tolerance is sufficient for most tasks, however, it can easily be relaxed or
tightened as needed.

## Bras and Kets

In Bra-Ket notation, the complex-valued column vector ``|\psi\rangle`` is referred
to as a ket.
Each ket has an associated dual ``\langle \psi |`` which is referred to as a bra.
A bra is represented by a complex-valued row vector and can be
constructed from a ket via the Hermitian adjoint, ``\langle\psi| = (|\psi\rangle)^{\dagger}``.

Bras and kets are commonly used to represent a quantum state where it is important
that the normalization ``\langle\psi|\psi\rangle = 1`` hold true.
Here ``\langle \; | \; \rangle`` denotes the inner product (dot product) between bra and ket.

```@docs
is_braket
Ket
Bra
```

The following methods have been extended to return bras and kets.

```@docs
*(ket :: Ket{<:Number}, bra :: Bra{<:Number})
adjoint(ket :: Ket{<:Number})
kron(kets :: Vararg{Ket{<:Number}}; atol=ATOL)
```

## Operators

```@docs
Operator
```

## Probabilities

```@docs
ProbabilityDistribution
is_probability_distribution
Probabilities
ConditionalDistribution
is_conditional_distribution
Conditionals
JointProbabilityDistribution
JointProbabilities
outcome_probabilities
```
