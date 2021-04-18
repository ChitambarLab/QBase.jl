# Ket Constructors
# qubit constructors
export bloch_qubit_ket
export mirror_symmetric_qubit_kets, planar_symmetric_qubit_kets
export trine_qubit_kets, bb84_qubit_kets, sic_qubit_kets

# 2-qubit constructors
export bell_kets

# general constructors
export computational_basis_kets
export generalized_bell_kets

"""
    computational_basis_kets( dim :: Int64 ) :: Vector{Ket{Int64}}

The computational basis vectors for the Hilbert space of dimension, `dim`.

```math
\\Psi = \\{ |j\\rangle \\}_{j=0}^{(d-1)}
```
"""
function computational_basis_kets(dim::Int64) :: Vector{Ket{Int64}}
    Ket.(computational_basis_vectors(dim))
end

"""
    bloch_qubit_ket( θ :: Real, ϕ :: Real ) :: Ket{Complex{Float64}}

Returns the qubit ket for the specified spherical coordinate on the surface of
bloch sphere, `(r=1, θ, ϕ)`:

```math
|\\psi\\rangle = \\cos(\\theta/2)|0\\rangle + e^{i\\phi}\\sin(\\theta/2)|1\\rangle
```

If the ket does not have a phase then a real-valued `Ket` is constructed as:

```julia
bloch_qubit_ket( θ :: Real ) :: Ket{Float64}
```

Inputs:

* `θ ∈ [0, 2π]`: polar angle (w.r.t z-axis)
* `ϕ ∈ [0, 2π]`: azimuthal angle (x-y plane)

A `DomainError` is thrown if inputs `θ` and/or `ϕ` do are not within the valid range.
"""
function bloch_qubit_ket(θ::Real, ϕ::Real) :: Ket{Complex{Float64}}
    if !(0 <= θ <= 2π)
        throw(DomainError(θ, "polar bloch-angle (θ) must be in domain [0,2π]"))
    elseif !(0 <= ϕ <= 2π)
        throw(DomainError(ϕ, "azimuthal angle (ϕ) must be in domain [0,2π]"))
    end

    Ket([cos(θ/2),exp(im*ϕ)*sin(θ/2)])
end

function bloch_qubit_ket(θ::Real) :: Ket{Float64}
    if !(-2π <= θ <= 2π)
        throw(DomainError(θ, "polar bloch-angle (θ) must be in domain [0,2π]"))
    end

    Ket([cos(θ/2),sin(θ/2)])
end

"""
    bell_kets() :: Vector{Ket{Float64}}

The Bell basis kets for maximally entangled bipartite qubit systems.
These kets are ordered as ``\\{|\\Phi^+\\rangle, |\\Phi^-\\rangle, |\\Psi^+\\rangle, |\\Psi^-\\rangle \\}``, where

```math
\\begin{matrix}
    |\\Phi^+\\rangle = \\frac{1}{\\sqrt{2}}(|00\\rangle + |11\\rangle), &
|\\Phi^-\\rangle = \\frac{1}{\\sqrt{2}}(|00\\rangle - |11\\rangle), \\\\
    |\\Psi^+\\rangle = \\frac{1}{\\sqrt{2}}(|01\\rangle + |10\\rangle), &
    |\\Psi^-\\rangle = \\frac{1}{\\sqrt{2}}(|01\\rangle - |10\\rangle). \\\\
\\end{matrix}
```
"""
bell_kets() :: Vector{Ket{Float64}} = Ket.([
    1/sqrt(2)*[1,0,0,1],
    1/sqrt(2)*[1,0,0,-1],
    1/sqrt(2)*[0,1,1,0],
    1/sqrt(2)*[0,1,-1,0]
])

"""
    generalized_bell_kets( dim :: Int64 ) :: Vector{Ket{Float64}}

The Bell basis for entangled bipartite quantum states each of dimension `dim`. Each state is constructed by

```math
|\\Psi^p_c\\rangle = \\frac{1}{\\sqrt{d}}\\sum_{j=0}^{d-1} e^{i2\\pi pj/d}|j\\rangle |\\mod(j+c,d)\\rangle
```

where ``p,c\\in \\{0,\\cdots, (d-1)\\}`` and ``d`` is `dim`. When iterated, ``c`` is the major
index and ``p`` is the minor index.

A `DomainError` is thrown if `dim ≥ 2` is not satisfied.
"""
function generalized_bell_kets(dim :: Int64) :: Vector{Ket{Complex{Float64}}}
    if !(dim ≥ 2)
        throw(DomainError(dim, "Hilbert space dimension must satisfy `dim ≥ 2`"))
    end

    basis = computational_basis_vectors(dim)

    kets = map(
        c -> map(
            p -> Ket(sum(
                j -> exp(im*2*π*p*j/dim)*kron(basis[j+1],basis[mod(j + c, dim) + 1]),
                0:dim-1
            )/sqrt(dim)),
        0:dim-1),
    0:dim-1)

    collect(flatten(kets))
end

"""
    mirror_symmetric_qubit_kets( θ :: Real ) :: Vector{Ket{Float64}}

Returns the triplet of qubit kets in the x-z plane of bloch sphere. The first ket
is  ``|0\\rangle`` and the other two are symmetric across the z-axis,
``|\\pm\\rangle = \\cos(\\theta)|0 \\rangle \\pm \\sin(\\theta)|1\\rangle``.

Input:

* `θ ∈ [0,π/2]`: the hilbert space angle between ``|0\\rangle`` and symmetric kets.
"""
function mirror_symmetric_qubit_kets(θ::Real) :: Vector{Ket{Float64}}
    if !(0 <= θ <= π/2)
        throw(DomainError(θ, "Hilbert space angle θ must satisfy `0 <= θ <= π/2`"))
    end

    ψ1 = Ket([1.,0.])
    ψ2 = bloch_qubit_ket(2*θ)
    ψ3 = bloch_qubit_ket(-2*θ)

    [ψ1,ψ2,ψ3]
end

"""
    planar_symmetric_qubit_kets( n :: Int64 ) :: Vector{Ket{Float64}}

Constructs a set of `n ``Ket`s oriented symmetrically in the x-z-plane.
Each ket is separated by a bloch angle of `2π/n`. The `Ket`s are constructed
with the form:

```math
|\\psi_j \\rangle = \\cos(j \\pi/n) | 0\\rangle + \\sin(j \\pi/n)|1\\rangle
```

where ``j \\in \\{0,\\cdots, (n-1)\\}``.

A `DomainError` is thrown if `n < 2`.
"""
function planar_symmetric_qubit_kets(n :: Int64) :: Vector{Ket{Float64}}
    if !(n ≥ 2)
        throw(DomainError(n, "n must satisfy `n ≥ 2`"))
    end

    map(θ -> (θ ≤ π) ? bloch_qubit_ket(θ) : bloch_qubit_ket(θ - 2π), 0:2π/n:2π*(n-1)/n+π/(2*n))
end

"""
    trine_qubit_kets() :: Vector{Ket{Float64}}

The triplet of `Ket`s separated by equal angles in the x-z plane of bloch sphere.

```math
    |\\psi_1\\rangle = |0\\rangle, \\quad
    |\\psi_2\\rangle = \\frac{1}{2}|0\\rangle + \\frac{\\sqrt{3}}{2}|1\\rangle, \\quad
    |\\psi_3\\rangle = \\frac{1}{2}|0\\rangle - \\frac{\\sqrt{3}}{2}|1\\rangle
```

```jldoctest
julia> trine_qubit_kets() == [[1.0, 0], [0.5, sqrt(3)/2], [0.5, -sqrt(3)/2]]
true
```
"""
trine_qubit_kets() :: Vector{Ket{Float64}} = Ket.([
    [1.0, 0], [0.5, sqrt(3)/2], [0.5, -sqrt(3)/2]
])

"""
    sic_qubit_kets() :: Vector{Ket{Complex{Float64}}}

The quadruplet of symmetric informationally complete (SIC) qubits. This set of qubits
correspond to the vertices of a tetrahedron inscribed on bloch sphere.

```math
\\begin{matrix}
    |\\psi_1\\rangle = |0\\rangle, & \\quad |\\psi_2\\rangle = \\frac{1}{\\sqrt{3}}|0\\rangle + \\sqrt{\\frac{2}{3}}|1\\rangle, \\\\
    |\\psi_3\\rangle = \\frac{1}{\\sqrt{3}}|0\\rangle + \\sqrt{\\frac{2}{3}}\\exp{i 2\\pi/3}|1\\rangle, & \\quad |\\psi_4\\rangle = \\frac{1}{\\sqrt{3}}|0\\rangle + \\sqrt{\\frac{2}{3}}\\exp{i 4\\pi/3}|1\\rangle
\\end{matrix}
```

```jldoctest
julia> sic_qubit_kets() == [
    [1., 0im],
    [1/sqrt(3), sqrt(2/3)+ 0im],
    [1/sqrt(3), sqrt(2/3)*exp(im*2π/3)],
    [1/sqrt(3), sqrt(2/3)*exp(im*4π/3)]
]
true
```
"""
sic_qubit_kets() :: Vector{Ket{Complex{Float64}}} = Ket.([
    [1., 0im], [1/sqrt(3), sqrt(2/3)+ 0im],
    [1/sqrt(3), sqrt(2/3)*exp(im*2π/3)], [1/sqrt(3), sqrt(2/3)*exp(im*4π/3)]
])

"""
    bb84_qubit_kets :: Vector{Ket{Float64}}

The quadruplet of qubit kets used in the BB84 Quantum Key Distribution protocol. The
states are ``|0\\rangle``, ``|+\\rangle``, ``|1\\rangle``, and ``|- \\rangle``.

```jldoctest
julia> bb84_qubit_kets() == [
    [1,0], [1,1]/sqrt(2), [0,1], [1,-1]/sqrt(2)
]
true
```
"""
bb84_qubit_kets() :: Vector{Ket{Float64}} = Ket.([
    [1,0], [1,1]/sqrt(2), [0,1], [1,-1]/sqrt(2)
])
