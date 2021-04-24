# Tutorial

## Quickstart

### Add the QBase.jl Package:

```julia
julia> using Pkg; Pkg.add("QBase")
```

### Import the QBase module

```@example tutorial
using QBase
```

## Algebra with Bras and Kets

See the [Bras and Kets](@ref) section for details.

### Create a [`Ket`](@ref)
```@example tutorial
ψ = Ket([1,0])
```

### Create a [`Bra`](@ref) by taking the adjoint.
```@example tutorial
ψ'
```

### Inner product between `Bra` and `Ket`
```@example tutorial
ψ'ψ
```

### Outer product between `Ket` and `Bra`
```@example tutorial
ψ*ψ'
```

## Quantum States

See [States](@ref) section for details.

### Constructing States
```@example tutorial
ρ_0 = State( [1 0;0 0] )
```
```@example tutorial
ρ_1 = State( [0 0;0 1] )
```

### Create Product States
```@example tutorial
kron( ρ_0, ρ_1 )
```

### Create a Bloch Qubit State
```@example tutorial
bloch_qubit_state( π/3, 7π/6 )
```

### Create a State from a Ket

```@example tutorial
pure_state(ψ)
```

### Create a Mixed State
```@example tutorial
mixed_state( [0.6, 0.4], [ρ_0, ρ_1] )
```

## Evolution of Quantum States

See [Evolution](@ref) section for details.

### Pauli Constants
```@example tutorial
[σI, σx, σy, σz]
```

### Unitary Evolution
```@example tutorial
evolve(σx, ρ_0)
```

```@example tutorial
evolve(σx, ψ)
```

### Channel Evolution


```@example tutorial
depolarizing_channel(ρ_0, 0.5)
```

## Measurement of Quantum States

See [Measurements](@ref) section for details.

### Constructing POVM Measurements
```@example tutorial
# measurement in z-basis
Π_Z = POVM([ [1 0;0 0], [0 0;0 1] ])
```

```@example tutorial
# measurement in  x-basis
Π_X = POVM([ [0.5 0.5;0.5 0.5], [0.5 -0.5;-0.5 0.5] ])
```

### Measurement Probabilities
Outcome [`Probabilities`](@ref) are obtained with the [`measure`](@ref) method.

```@example tutorial
measure( Π_Z, ρ_0 )
```

```@example tutorial
measure( Π_X, ρ_0 )
```

```@example tutorial
measure( Π_Z, [ρ_0, ρ_1] )
```

### Quantum Information

See the [Information Theory](@ref) section for a complete list of methods.

#### Von Neumann Entropy of Bell state

```@example tutorial
# maximally entangled bipartite qubit state
ρ_bell = bell_states()[1]

von_neumann_entropy( ρ_bell )
```
```@example tutorial
# maximally mixed two-qubit state
ρ_bell_mix = mixed_state( [1,1,1,1]/4, bell_states() )

von_neumann_entropy( ρ_bell_mix )
```

## Advanced Examples

### Absolute Tolerance

In some cases, it may be desirable to relax or tighten the tolerated numerical error.
This example demonstrates  how to pass the `atol` parameter to the `State` constructor.
This parameter can be used for any type defined in this project.
The absolute tolerance should be used with caution as it may delegitimize computed results.

```@example tutorial
ϵ = 1e-5    # introducing a small error

# state creation fails due to non-hermiticity
try
    State([1 ϵ;-ϵ 0])
catch err
    println( err )
end
```

```@example tutorial
# the atol parameter overrides the non-hermiticity failure
State([1 ϵ;-ϵ 0], atol=1e-4)
```

### Random States

```@example tutorial
# Create a 5x5 random unitary
U_rand = random_unitary( 5 )
```

```@example tutorial
# evolve the |0> Ket
ψ_rand = U_rand * Ket( [1, 0, 0, 0, 0] )
```

```@example tutorial
# create a pure state from the random ket
ρ_rand = pure_state(ψ_rand)
```

### Noisy Quantum Channel

Computing measurement probabilities in a noiseless channel:

```@example tutorial
# create bipartite Bell states
ρ_bell_states = bell_states()

# create a measurement in the same Bell basis
Π_bell_basis = POVM( ρ_bell_states )

# compute the ideal measurement probabilities
measure( Π_bell_basis, ρ_bell_states )
```

Computing measurement probabilities with a noisy channel:

```@example tutorial
# create a depolarizing channel which mixes in 50% white noise
𝒩(ρ) = depolarizing_channel(ρ, 0.5)

# Add noise to quantum states
ρ_noisy_bell_states = 𝒩.(ρ_bell_states)

# compute the noisy measurement probabilities
measure( Π_bell_basis, ρ_noisy_bell_states  )
```

### Reduced Density States

The reduced density matrix of a Bell state is a maximally mixed state.

```@example tutorial
# bipartite entangled qutrit states
ρ = generalized_bell_states( 3 )[1]

# tracing out the first qutrit subsystem
# creates a reduced density matrix State
partial_trace(ρ, [3,3], 1)
```
