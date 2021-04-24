# constants
export σI, σx, σy, σz

# contructors
export qubit_rotation, random_unitary

"""
    σI :: Unitary{Int64}

Identity unitary:

```jldoctest
julia> σI
2×2 Unitary{Int64}:
 1  0
 0  1
```
"""
const σI = Unitary([1 0;0 1])

"""
    σx :: Unitary{Int64}

Pauli-X unitary:

```jldoctest
julia> σx
2×2 Unitary{Int64}:
 0  1
 1  0
```
"""
const σx = Unitary([0 1;1 0])

"""
    σy :: Unitary{Complex{Int64}}

Pauli-Y unitary:

```jldoctest
julia> σy
2×2 Unitary{Complex{Int64}}:
 0+0im  0-1im
 0+1im  0+0im
```
"""
const σy = Unitary([0 -im;im 0])

"""
    σz :: Unitary{Int64}

Pauli-Z unitary:

```
julia> σz
2×2 Unitary:
 1   0
 0  -1
```
"""
const σz = Unitary([1 0;0 -1])

"""
    qubit_rotation( θ :: Real; axis="x" :: String ) :: Unitary{Complex{Float64}}

Returns a unitary which performs a qubit rotation along bloch sphere. `θ` designates
the angle of rotation and `axis` (`"x"`, `"y"`, `"z"`) designates the cartesian
axis about which the qubit is rotated.
"""
function qubit_rotation(θ :: Real; axis="x" :: String)::Unitary{Complex{Float64}}

    pauli = [1 0; 0 1]
    if axis == "x"
        pauli = σx
    elseif axis == "y"
        pauli = σy
    elseif axis == "z"
        pauli = σz
    end

    Unitary(exp(-im * pauli * θ/2))
end

"""
    random_unitary( d :: Int64 ) :: Unitary{Complex{Float64}}

Constructs a `d x d` random unitary matrix according to the Haar measure.
"""
function random_unitary(d :: Int64) :: Unitary{Complex{Float64}}
    generator = Haar(2) # Specifies a complex hermitian unitary matrix
    Unitary(rand(generator, d))
end
