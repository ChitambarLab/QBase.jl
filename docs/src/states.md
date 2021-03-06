```@meta
CurrentModule = QBase
```
# States

A quantum state is represented by a density operator.
This matrix describes the wave function of the quantum system.

```@docs
is_density_matrix
State
```

## State Operations

```@docs
kron(states :: Vararg{State}; atol=ATOL :: Float64)
partial_trace(ρ::State, system::Vector{Int64}, id::Int64; atol=ATOL :: Float64)
eigvals(ρ :: State)
is_pure
is_mixed
```

## State Constructors

### Singlet States

```@docs
pure_state
mixed_state
bloch_qubit_state
```

### Ensemble States

```@docs
computational_basis_states
bell_states
generalized_bell_states
planar_symmetric_qubit_states
trine_qubit_states
mirror_symmetric_qubit_states
sic_qubit_states
bb84_qubit_states
```
