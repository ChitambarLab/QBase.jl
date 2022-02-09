```@meta
CurrentModule = QBase
```
# Evolution

```@docs
evolve
*(U :: Unitary, ket :: Ket)
```

## Unitaries

```@docs
Unitary
is_unitary
```

### Unitary Operations

```@docs
*(unitaries :: Vararg{Unitary}; atol=ATOL :: Float64)
kron(unitaries :: Vararg{Unitary}; atol=ATOL :: Float64)
adjoint(U :: Unitary)
```

### Unitary Constructors

```@docs
σI
σx
σy
σz
qubit_rotation
random_unitary
```

## Channels

```@docs
replacer_channel
depolarizing_channel
erasure_channel
```

## Kraus Channels

```@docs
KrausChannel
is_kraus_channel
kraus_evolve
```

## Choi Operators

```@docs
ChoiOp
is_choi_matrix
choi_matrix
choi_evolve
```
