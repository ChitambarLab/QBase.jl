var documenterSearchIndex = {"docs":
[{"location":"Unitaries/#","page":"Unitaries","title":"Unitaries","text":"CurrentModule = QBase.Unitaries","category":"page"},{"location":"Unitaries/#QBase.Unitaries-1","page":"Unitaries","title":"QBase.Unitaries","text":"","category":"section"},{"location":"Unitaries/#","page":"Unitaries","title":"Unitaries","text":"Unitaries\nis_unitary\nAbstractUnitary\nUnitary\nQubitUnitary","category":"page"},{"location":"Unitaries/#QBase.Unitaries","page":"Unitaries","title":"QBase.Unitaries","text":"Quantum states evolve under unitary transformations. The QBase.Unitaries submodule provides:\n\nTypes and Constructors for unitary operators.\n\n\n\n\n\n","category":"module"},{"location":"Unitaries/#QBase.Unitaries.is_unitary","page":"Unitaries","title":"QBase.Unitaries.is_unitary","text":"is_unitary( U :: Matrix ) :: Bool\n\nReturns true if matrix U is unitary. The hermitian adjoint of a unitary matrix is its inverse:\n\nU'U == I where I is the identity matrix.\n\nA unitary matrix must be square. A DomainError is thrown if input U is not square.\n\n\n\n\n\n","category":"function"},{"location":"Unitaries/#QBase.Unitaries.AbstractUnitary","page":"Unitaries","title":"QBase.Unitaries.AbstractUnitary","text":"AbstractUnitary <: AbstractMatrix{Complex{Float64}}\n\nThe abstract type representing unitary operators. An AbstractUnitary cannot be instantiated, it serves as a supertype from which concrete types are derived.\n\n\n\n\n\n","category":"type"},{"location":"Unitaries/#QBase.Unitaries.Unitary","page":"Unitaries","title":"QBase.Unitaries.Unitary","text":"Unitary( U :: Matrix ) <: AbstractUnitary\n\nUnitary matrices represent the physical evolution of quantum states. The Constructor, Unitary(U), throws a DomainError if the provided matrix, U is not unitary.\n\n\n\n\n\n","category":"type"},{"location":"Unitaries/#QBase.Unitaries.QubitUnitary","page":"Unitaries","title":"QBase.Unitaries.QubitUnitary","text":"QubitUnitary( U :: Matrix ) <: AbstractUnitary\n\nConstructs a 2x2 unitary for qubit evolution. Throws a DomainError if input U is not of dimension 2x2 or if U is not unitary.\n\n\n\n\n\n","category":"type"},{"location":"Unitaries/#","page":"Unitaries","title":"Unitaries","text":"paulis\nσx\nσy\nσz","category":"page"},{"location":"Unitaries/#QBase.Unitaries.paulis","page":"Unitaries","title":"QBase.Unitaries.paulis","text":"paulis :: Vector{QubitUnitary}\n\nReturns a vector containing the three qubit pauli matrices, [σx, σy, σz].\n\n\n\n\n\n","category":"constant"},{"location":"Unitaries/#QBase.Unitaries.σx","page":"Unitaries","title":"QBase.Unitaries.σx","text":"σx :: QubitUnitary\n\nPauli-X unitary:\n\njulia> QBase.σx\n2×2 QBase.Unitaries.QubitUnitary:\n 0.0+0.0im  1.0+0.0im\n 1.0+0.0im  0.0+0.0im\n\n\n\n\n\n","category":"constant"},{"location":"Unitaries/#QBase.Unitaries.σy","page":"Unitaries","title":"QBase.Unitaries.σy","text":"σy :: QubitUnitary\n\nPauli-Y unitary:\n\njulia> QBase.σy\n2×2 QBase.Unitaries.QubitUnitary:\n 0.0+0.0im  0.0-1.0im\n 0.0+1.0im  0.0+0.0im\n\n\n\n\n\n","category":"constant"},{"location":"Unitaries/#QBase.Unitaries.σz","page":"Unitaries","title":"QBase.Unitaries.σz","text":"σz :: QubitUnitary\n\nPauli-Z unitary:\n\njulia> QBase.σz\n2×2 QBase.Unitaries.QubitUnitary:\n 1.0+0.0im   0.0+0.0im\n 0.0+0.0im  -1.0+0.0im\n\n\n\n\n\n","category":"constant"},{"location":"Unitaries/#","page":"Unitaries","title":"Unitaries","text":"qubit_rotation","category":"page"},{"location":"Unitaries/#QBase.Unitaries.qubit_rotation","page":"Unitaries","title":"QBase.Unitaries.qubit_rotation","text":"qubit_rotation( θ :: Real; axis=\"x\" :: String ) :: QubitUnitary\n\nReturns a unitary which performs a qubit rotation along bloch sphere. θ designates the angle of rotation and axis (\"x\", \"y\", \"z\") designates the cartesian axis about which the qubit is rotated.\n\n\n\n\n\n","category":"function"},{"location":"States/#","page":"States","title":"States","text":"CurrentModule = QBase.States","category":"page"},{"location":"States/#QBase.States-1","page":"States","title":"QBase.States","text":"","category":"section"},{"location":"States/#","page":"States","title":"States","text":"States","category":"page"},{"location":"States/#QBase.States","page":"States","title":"QBase.States","text":"The QBase.States submodule provides:\n\nabstract and concrete types for quantum state representations.\nA catalog of constructors for instantiating quantum states.\n\n\n\n\n\n","category":"module"},{"location":"States/#","page":"States","title":"States","text":"In discrete and finite quantum systems, states can be represented computationally with vectors (Bra-Ket Representation) or matrices (Density Matrix Representation).","category":"page"},{"location":"States/#Bra-Ket-Representation-1","page":"States","title":"Bra-Ket Representation","text":"","category":"section"},{"location":"States/#","page":"States","title":"States","text":"A quantum state may be represented by a vector on a complex-valued Hilbert space. In Bra-Ket notation, a quantum state is represented by a ket denoted psirangle. A ket can simply be understood as a column vector. Each ket has an associated dual vector known as a bra. A bra is may be represented by a row vector and can be constructed from a ket via Hermitian adjoint, langlepsi = (psirangle)^dagger.","category":"page"},{"location":"States/#","page":"States","title":"States","text":"It is essential that a quantum states are normalized, langlepsipsirangle = 1, where langle    rangle denotes the inner product (dot product) between bra and ket.","category":"page"},{"location":"States/#","page":"States","title":"States","text":"The QBase.States module provides a validation method, is_ket() for checking whether a vector satisfies the requirements for being a quantum state ket.","category":"page"},{"location":"States/#","page":"States","title":"States","text":"is_ket","category":"page"},{"location":"States/#QBase.States.is_ket","page":"States","title":"QBase.States.is_ket","text":"is_ket( ψ :: Vector ) :: Bool\n\nReturns true if vector ψ is a valid ket representation of a quamtum state:\n\nψ is a real or complex-valued vector.\nψ is normalized with respect to the bra-ket inner prodcut (ψ'ψ == 0).\n\n\n\n\n\n","category":"function"},{"location":"States/#Ket-Types-1","page":"States","title":"Ket Types","text":"","category":"section"},{"location":"States/#","page":"States","title":"States","text":"AbstractKet","category":"page"},{"location":"States/#QBase.States.AbstractKet","page":"States","title":"QBase.States.AbstractKet","text":"AbstractKet <: AbstractVector{Complex{Float64}}\n\nThe abstract type representing a quantum state ket. Since kets are contained within a complex-valued hilbert space, they are appropriately AbstractVector{Complex{Float64}} subtypes. An AbstractKet cannot be instantiated, it serves as a supertype from which ket types are defined.\n\n\n\n\n\n","category":"type"},{"location":"States/#","page":"States","title":"States","text":"The QBase.States module provides two concrete subtypes of AbstractKet:","category":"page"},{"location":"States/#","page":"States","title":"States","text":"Ket\nQubitKet","category":"page"},{"location":"States/#QBase.States.Ket","page":"States","title":"QBase.States.Ket","text":"Ket( ψ :: Vector{Complex{Float64}} ) <: AbstractKet\n\nA ket representation of a general quantum state. When given invalid input, the constructor, Ket(ψ), throws:\n\nDomainError – If ψ is not normalized.\nMethodError – If ψ is not a column vector ([a,b,c] or [a;b;c])\n\n\n\n\n\n","category":"type"},{"location":"States/#QBase.States.QubitKet","page":"States","title":"QBase.States.QubitKet","text":"QubitKet( ψ :: Vector{Complex{Float64}} ) <: AbstractKet\n\nA ket representation of a 2-dimensional quantum state. When given invalid input, the constructor, QubitKet(ψ), throws:\n\nDomainError – If ψ is not normalized.\nMethodError – If ψ is not a column vector ([a,b] or [a;b]).\n\n\n\n\n\n","category":"type"},{"location":"States/#Ket-Constructors-1","page":"States","title":"Ket Constructors","text":"","category":"section"},{"location":"States/#Singlet-States-1","page":"States","title":"Singlet States","text":"","category":"section"},{"location":"States/#","page":"States","title":"States","text":"bloch_qubit_ket","category":"page"},{"location":"States/#QBase.States.bloch_qubit_ket","page":"States","title":"QBase.States.bloch_qubit_ket","text":"bloch_qubit_ket( θ :: Real, ϕ :: Real ) :: QubitKet\n\nReturns the qubit ket for the specified spherical coordinate on the surface of bloch sphere, (r=1, θ, ϕ).\n\nInputs:\n\nθ ∈ [0, π]: polar angle (w.r.t z-axis)\nϕ ∈ [0, 2π]: azimuthal angle (x-y plane)\n\nA DomainError is thrown if inputs θ and/or ϕ do are not within the valid range.\n\n\n\n\n\n","category":"function"},{"location":"States/#Triplet-States-1","page":"States","title":"Triplet States","text":"","category":"section"},{"location":"States/#","page":"States","title":"States","text":"trine_qubit_kets\nmirror_symmetric_qubit_kets","category":"page"},{"location":"States/#QBase.States.trine_qubit_kets","page":"States","title":"QBase.States.trine_qubit_kets","text":"trine_qubit_kets :: Vector{QubitKet}\n\nThe triplet of kets representing three quantum states separated by equal angles in the z-x plane of bloch sphere.\n\njulia> QBase.trine_qubit_kets == [[1.0, 0], [0.5, sqrt(3)/2], [0.5, -sqrt(3)/2]]\ntrue\n\n\n\n\n\n","category":"constant"},{"location":"States/#QBase.States.mirror_symmetric_qubit_kets","page":"States","title":"QBase.States.mirror_symmetric_qubit_kets","text":"mirror_symmetric_qubit_kets( θ :: Real ) :: Vector{QubitKet}\n\nReturns the triplet of qubit kets in the z-x plane of bloch sphere. The first ket is  0rangle and the other two are symmetric across the z-axis.\n\nInput:\n\nθ ∈ [0,π/2]: the hilbert space angle between |0> and symmetric kets.\n\n\n\n\n\n","category":"function"},{"location":"States/#Density-Matrix-Representation-1","page":"States","title":"Density Matrix Representation","text":"","category":"section"},{"location":"States/#","page":"States","title":"States","text":"Quantum states can be represented in matrix form.","category":"page"},{"location":"States/#","page":"States","title":"States","text":"is_density_matrix\nAbstractDensityMatrix\nDensityMatrix\nQubit","category":"page"},{"location":"States/#QBase.States.is_density_matrix","page":"States","title":"QBase.States.is_density_matrix","text":"is_density_matrix( ρ :: Matrix ) :: Bool\n\nReturns true if input ho is:\n\nHermitian\nPositive Semi-Definite\nTrace[ρ] = 1 (normalization)\n\n\n\n\n\n","category":"function"},{"location":"States/#QBase.States.AbstractDensityMatrix","page":"States","title":"QBase.States.AbstractDensityMatrix","text":"AbstractDensityMatrix <: AbstractMatrix{Complex{Float64}}\n\nThe abstract type representing all density matrices.\n\n\n\n\n\n","category":"type"},{"location":"States/#QBase.States.DensityMatrix","page":"States","title":"QBase.States.DensityMatrix","text":"DensityMatrix( ρ :: Matrix{Complex{Float64}} ) <: AbstractDensityMatrix\n\nThe density matrix representation of a quantum state. The constructor, DensityMatrix(ρ) throws a DomainError if is_density_matrix(ρ) is false.\n\n\n\n\n\n","category":"type"},{"location":"States/#QBase.States.Qubit","page":"States","title":"QBase.States.Qubit","text":"Qubit( ρ :: Matrix{Complex{Float64}} ) <: AbstractDensityMatrix\n\nThe 2x2 density matrix representation of a qubit. The constructor, Qubit(ρ) throws a DomainError if is_density_matrix(ρ) is false or if ρ is not 2x2 in dimension.\n\n\n\n\n\n","category":"type"},{"location":"States/#Density-Matrix-Constructors-1","page":"States","title":"Density Matrix Constructors","text":"","category":"section"},{"location":"States/#","page":"States","title":"States","text":"The density matrix rho can be constructed from ket psirangle by taking the outer-product of the ket with itself, psiranglelanglepsi = rho.","category":"page"},{"location":"States/#","page":"States","title":"States","text":"pure_state\npure_qubit","category":"page"},{"location":"States/#QBase.States.pure_state","page":"States","title":"QBase.States.pure_state","text":"pure_state( ψ :: AbstractKet ) :: Qubit\n\nA state is considered \"pure\" if it is rank-one. A rank-one density matrix is constructed by  taking the outer-product of a ket state.\n\nThe method alternatively accepts a Vector input.\n\npure_state( ψ :: Vector ) :: Qubit\n\nA DomainError is thrown if ψ is not a valid ket.`\n\n\n\n\n\n","category":"function"},{"location":"States/#QBase.States.pure_qubit","page":"States","title":"QBase.States.pure_qubit","text":"pure_qubit( ψ :: AbstractKet ) :: Qubit\n\nA qubit is considered \"pure\" if it is rank-one. A rank-one density matrix is constructed by  taking the outer-product of a ket state.\n\nThe method alternatively accepts a Vector input.\n\npure_qubit( ψ :: Vector ) :: Qubit\n\nA DomainError is thrown if ψ is not a valid ket.`\n\n\n\n\n\n","category":"function"},{"location":"States/#","page":"States","title":"States","text":"The rank of the density matrix can be greater than 1. If a density matrix has rank greater than one, we call the state mixed.","category":"page"},{"location":"States/#","page":"States","title":"States","text":"mixed_state\nmixed_qubit\nbloch_qubit","category":"page"},{"location":"States/#QBase.States.mixed_state","page":"States","title":"QBase.States.mixed_state","text":"mixed_state( priors :: QMath.Marginals, ρ_states :: Vector{<:AbstractDensityMatrix} ) :: DensityMatrix\n\nConstructs the statistical mixture (weighted average) of quantum states. The method accepts states as type DensityMatrix or subbtypes AbstractDensityMatrix.\n\n\n\n\n\n","category":"function"},{"location":"States/#QBase.States.mixed_qubit","page":"States","title":"QBase.States.mixed_qubit","text":"mixed_qubit( priors :: QMath.Marginals, ρ_states :: Vector{Qubit} ) :: Qubit\n\nConstructs the statistical mixture (weighted average) of qubits.\n\n\n\n\n\n","category":"function"},{"location":"States/#QBase.States.bloch_qubit","page":"States","title":"QBase.States.bloch_qubit","text":"Returns the qubit density matrix for quantum state described by a coordinate on bloch sphere.\n\nSpherical Coordinates:\n\nStates on the surface of bloch sphere may be described by spherical coordinates.\n\nbloch_qubit(θ::Real, ϕ::Real) :: Qubit\n\nθ ∈ [0, π]: polar angle (w.r.t z-axis).\nϕ ∈ [0, 2π]: azimuthal angle (x-y plane)\n\nCartesian Coordinates:\n\nStates within the volume of bloch sphere may be described in cartesian coordinates.\n\nbloch_qubit(x::Real, y::Real, z::Real) :: Qubit\n\nwhere x, y, and z are constrained to the unit sphere, 0 <= norm([x,y,z]) <= 1.\n\nA DomainError is thrown if the coordinates (x,y,z) do  not adhere to constraints.\n\n\n\n\n\n","category":"function"},{"location":"States/#Triplets-1","page":"States","title":"Triplets","text":"","category":"section"},{"location":"States/#","page":"States","title":"States","text":"trine_qubits\nmirror_symmetric_qubits","category":"page"},{"location":"States/#QBase.States.trine_qubits","page":"States","title":"QBase.States.trine_qubits","text":"trine_qubits :: Vector{Qubit}\n\nReturns the qubit trine states in density matrix form.\n\n\n\n\n\n","category":"constant"},{"location":"States/#QBase.States.mirror_symmetric_qubits","page":"States","title":"QBase.States.mirror_symmetric_qubits","text":"mirror_symmetric_qubits(θ) :: Vector{Qubit}\n\nReturns a set of 3 mirror symmetric qubit density matrices. The first state is 0ranglelangle 0 the other two are symmetric about the  0rangle axis.\n\nInput:\n\nθ ∈ [0,π/2]: the hilbert space angle between 0rangle and psi_23rangle.\n\n\n\n\n\n","category":"function"},{"location":"States/#Quadruplets-1","page":"States","title":"Quadruplets","text":"","category":"section"},{"location":"States/#","page":"States","title":"States","text":"sic_qubits\nbb84_qubits","category":"page"},{"location":"States/#QBase.States.sic_qubits","page":"States","title":"QBase.States.sic_qubits","text":"sic_qubits :: Vector{Qubit}\n\nThe quadruplet of symmetric informationally complete (SIC) qubits. The qubits are the vertices of a tetrahedron inscribed on bloch sphere.\n\njulia> QBase.sic_qubits\n4-element Array{QBase.States.Qubit,1}:\n [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]\n [0.33333333333333337 + 0.0im 0.47140452079103173 + 0.0im; 0.47140452079103173 + 0.0im 0.6666666666666666 + 0.0im]\n [0.33333333333333337 + 0.0im -0.2357022603955158 - 0.4082482904638631im; -0.2357022603955158 + 0.4082482904638631im 0.6666666666666666 + 0.0im]\n [0.33333333333333337 - 0.0im -0.2357022603955158 + 0.4082482904638631im; -0.2357022603955158 - 0.4082482904638631im 0.6666666666666666 - 0.0im]\n\n\n\n\n\n","category":"constant"},{"location":"States/#QBase.States.bb84_qubits","page":"States","title":"QBase.States.bb84_qubits","text":"bb84_qubits :: Vector{Qubit}\n\nThe quadruplet of qubits used in the BB84 Quantum Key Distribution protocol. The states are 0ranglelangle 0, 1ranglelangle 1, +ranglelangle +, and - ranglelangle -.\n\njulia> QBase.bb84_qubits\n4-element Array{QBase.States.Qubit,1}:\n [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]\n [0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]\n [0.5 + 0.0im 0.5 + 0.0im; 0.5 + 0.0im 0.5 + 0.0im]\n [0.5 + 0.0im -0.5 + 0.0im; -0.5 + 0.0im 0.5 + 0.0im]\n\n\n\n\n\n","category":"constant"},{"location":"overview/#QBase.jl-Overview-1","page":"Overview","title":"QBase.jl - Overview","text":"","category":"section"},{"location":"overview/#","page":"Overview","title":"Overview","text":"A base library for quantum information.","category":"page"},{"location":"overview/#","page":"Overview","title":"Overview","text":"QBase.jl provides:","category":"page"},{"location":"overview/#","page":"Overview","title":"Overview","text":"Types for quantum data structures.\nA catalog of methods for creating quantum states, POVMs, and unitaries.\nMethods for evolving quantum systems and performing quantum measurements.","category":"page"},{"location":"Information/#","page":"Information","title":"Information","text":"CurrentModule = QBase.Information","category":"page"},{"location":"Information/#QBase.Information-1","page":"Information","title":"QBase.Information","text":"","category":"section"},{"location":"Information/#","page":"Information","title":"Information","text":"Information","category":"page"},{"location":"Information/#QBase.Information","page":"Information","title":"QBase.Information","text":"Methods for quantifying information and randomness.\n\nEntropy quantifiers are taken with respect to base-2 logarithms. The entropy is understood as the number of bits {0,1} required to communicate a random result with certainty.\n\n\n\n\n\n","category":"module"},{"location":"Information/#Entropy-1","page":"Information","title":"Entropy","text":"","category":"section"},{"location":"Information/#","page":"Information","title":"Information","text":"shannon_entropy\nvon_neumann_entropy\njoint_entropy\nconditional_entropy","category":"page"},{"location":"Information/#QBase.Information.shannon_entropy","page":"Information","title":"QBase.Information.shannon_entropy","text":"shannon_entropy( probabilities :: QMath.Marginals ) :: Float64\n\nThe classical entropy of a probability distribution.\n\n\n\n\n\n","category":"function"},{"location":"Information/#QBase.Information.von_neumann_entropy","page":"Information","title":"QBase.Information.von_neumann_entropy","text":"von_neumann_entropy( ρ :: States.DensityMatrix ) :: Float64\n\nThe von neumann entropy of a density matrix.\n\n\n\n\n\n","category":"function"},{"location":"Information/#QBase.Information.joint_entropy","page":"Information","title":"QBase.Information.joint_entropy","text":"joint_entropy(priors :: QMath.Marginals, conditionals :: QMath.Conditionals) :: Float64\n\nReturns the entropy for the union of pdf P(xy).\n\n\n\n\n\n","category":"function"},{"location":"Information/#QBase.Information.conditional_entropy","page":"Information","title":"QBase.Information.conditional_entropy","text":"conditional_entropy(priors::QMath.Marginals, conditionals::QMath.Conditionals) :: Float64\n\nReturns the conditional entropy for the system with specified priors and conditionals.\n\n\n\n\n\n","category":"function"},{"location":"Information/#Information-1","page":"Information","title":"Information","text":"","category":"section"},{"location":"Information/#","page":"Information","title":"Information","text":"holevo_bound\nholevo_information\nmutual_information","category":"page"},{"location":"Information/#QBase.Information.holevo_bound","page":"Information","title":"QBase.Information.holevo_bound","text":"holevo_bound( priors :: QMath.Marginals, ρ_states :: Vector{<:AbstractDensityMatrix} ) :: Float64\n\nComputes the upper bound of a quantum channel's information capacity. The information shared through a quantum channel cannot exceed a classical channel of the same dimension.\n\n\n\n\n\n","category":"function"},{"location":"Information/#QBase.Information.holevo_information","page":"Information","title":"QBase.Information.holevo_information","text":"holevo_information( priors :: QMath.Marginals, ρ_states :: Vector{<:AbstractDensityMatrix} ) :: Float64\n\nComputes the holevo (mutual) information shared through a quantum channel.\n\n\n\n\n\n","category":"function"},{"location":"Information/#QBase.Information.mutual_information","page":"Information","title":"QBase.Information.mutual_information","text":"mutual_information( priors :: QMath.Marginals, conditionals :: QMath.Conditionals ) :: Float64\n\nThe entropy of the overlap between p(x) and p(y). The information shared from y to x.\n\n\n\n\n\n","category":"function"},{"location":"Information/#State-Discrimination-1","page":"Information","title":"State Discrimination","text":"","category":"section"},{"location":"Information/#","page":"Information","title":"Information","text":"success_probability\nerror_probability","category":"page"},{"location":"Information/#QBase.Information.success_probability","page":"Information","title":"QBase.Information.success_probability","text":"success_probability(priors::QMath.Marginals, ρ_states::Vector{<:States.AbstractDensityMatrix}, Π::Observables.AbstractPOVM) :: Float64\n\nThe probability of correctly distinguishing quantum states with the specifed POVM.\n\n\n\n\n\n","category":"function"},{"location":"Information/#QBase.Information.error_probability","page":"Information","title":"QBase.Information.error_probability","text":"error_probability(priors::QMath.Marginals, ρ_states::Vector{<:States.AbstractDensityMatrix}, Π::Observables.AbstractPOVM) :: Float64\n\nThe probability of incorrectly distinguishing quantum states with the specifed POVM.\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Home","title":"Home","text":"CurrentModule = QBase","category":"page"},{"location":"#QBase-1","page":"Home","title":"QBase","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"warning: Under Development\nMissing documentation and broken code exists. Breaking changes will be made in future commits","category":"page"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [QBase]","category":"page"},{"location":"Observables/#","page":"Observables","title":"Observables","text":"CurrentModule = QBase.Observables","category":"page"},{"location":"Observables/#QBase.Observables-1","page":"Observables","title":"QBase.Observables","text":"","category":"section"},{"location":"Observables/#","page":"Observables","title":"Observables","text":"Observables","category":"page"},{"location":"Observables/#QBase.Observables","page":"Observables","title":"QBase.Observables","text":"Observables describe measurable quantities of a quantum system. Quantum measurement is a probabilistic process. The measurement outcomes are not definite, but described by conditional probabilities.\n\nThe QBase.Observables submodule provides types and constructors for representing quantum observables.\n\n\n\n\n\n","category":"module"},{"location":"Observables/#Types-1","page":"Observables","title":"Types","text":"","category":"section"},{"location":"Observables/#","page":"Observables","title":"Observables","text":"AbstractPOVM\nis_povm\nPOVM\nQubitPOVM","category":"page"},{"location":"Observables/#QBase.Observables.AbstractPOVM","page":"Observables","title":"QBase.Observables.AbstractPOVM","text":"AbstractPOVM <: AbstractMatrix{Complex{Float64}}\n\nThe abstract type representing positive-operator valued measures (POVMs). An AbstractPOVM cannot be instantiated, it serves as a supertype from which concrete types are derived.\n\n\n\n\n\n","category":"type"},{"location":"Observables/#QBase.Observables.is_povm","page":"Observables","title":"QBase.Observables.is_povm","text":"is_povm( Π :: Vector ) :: Bool\n\nReturns true if Π is a POVM. The following constraints must be satisfied:\n\nEach POVM element is hermitian\nEach POVM element positive semi-definite\nThe POVM is complete: sum(Π) == identity matrix\n\n\n\n\n\n","category":"function"},{"location":"Observables/#QBase.Observables.POVM","page":"Observables","title":"QBase.Observables.POVM","text":"POVM( Π :: Vector{Matrix} ) <: AbstractPOVM\n\nPositve-operator valued measures (POVMs) represent a general quantum measurement. Each POVM-element is a hermitian, positive-semidefinite matrix. The sum of all POVM-elements yields the identity matrix. The constructor, POVM(Π) throws a DomainError if the provided array of matrices, Π is not a valid POVM.\n\n\n\n\n\n","category":"type"},{"location":"Observables/#QBase.Observables.QubitPOVM","page":"Observables","title":"QBase.Observables.QubitPOVM","text":"QubitPOVM( Π :: Vector{Matrix} ) <: AbstractPOVM\n\nA general qubit measurement. A DomainError is thrown if Π does not contain 2x2 elements or if it is not a valid POVM.\n\n\n\n\n\n","category":"type"},{"location":"Observables/#Constructors-1","page":"Observables","title":"Constructors","text":"","category":"section"},{"location":"Observables/#","page":"Observables","title":"Observables","text":"mirror_symmetric_qubit_3povm\nasymmetric_qubit_3povm\ntrine_qubit_povm\nsic_qubit_povm\nsqrt_povm","category":"page"},{"location":"Observables/#QBase.Observables.mirror_symmetric_qubit_3povm","page":"Observables","title":"QBase.Observables.mirror_symmetric_qubit_3povm","text":"mirror_symmetric_qubit_3povm( θ :: Real ) :: QubitPOVM\n\nConstructs a QubitPOVM aligned with the three mirror symmetric qubit states. The first measurement is aligned with the 0rangle state and the remaining two are symmetric across the z-axis.\n\nA DomainError is thrown if argument θ ∉ [π/4,π/2].\n\n\n\n\n\n","category":"function"},{"location":"Observables/#QBase.Observables.asymmetric_qubit_3povm","page":"Observables","title":"QBase.Observables.asymmetric_qubit_3povm","text":"asymmetric_qubit_3povm( θ1::Real, θ2::Real ) :: QubitPOVM\n\nConstructs a general non-orthogonal 3-element QubitPOVM.\n\nInputs:\n\nθ1 ∈ (0,π/2] or [-π/2,0), the angle element 2 makes with 0rangle state.\nθ2 ∈ [-π/2,0) or (0,π/2], the angle element 3 makes with 0rangle state.\n\nA DomainError is thrown if θ1 and θ2 are not in the valid ranges. Note that one angle must be positive and the other negative.\n\n\n\n\n\n","category":"function"},{"location":"Observables/#QBase.Observables.trine_qubit_povm","page":"Observables","title":"QBase.Observables.trine_qubit_povm","text":"trine_qubit_povm :: QubitPOVM\n\nThe POVM with elements parallel to the trine qubit states.\n\n\n\n\n\n","category":"constant"},{"location":"Observables/#QBase.Observables.sic_qubit_povm","page":"Observables","title":"QBase.Observables.sic_qubit_povm","text":"sic_qubit_povm :: QubitPOVM\n\nThe POVM with elements parallel to the symmetric informationally complete (SIC) qubits.\n\n\n\n\n\n","category":"constant"},{"location":"Observables/#QBase.Observables.sqrt_povm","page":"Observables","title":"QBase.Observables.sqrt_povm","text":"sqrt_povm(priors :: QMath.Marginals, states :: Vector{<:States.AbstractDensityMatrix}) :: POVM\n\nReturns the \"pretty good\" square-root povm for the given density operator states and priors.\n\n\n\n\n\n","category":"function"},{"location":"Observables/#Quantum-Measurement-1","page":"Observables","title":"Quantum Measurement","text":"","category":"section"},{"location":"Observables/#","page":"Observables","title":"Observables","text":"kraus_operators\nmeasurement_probabilities\nnaimark_dilation","category":"page"},{"location":"Observables/#QBase.Observables.kraus_operators","page":"Observables","title":"QBase.Observables.kraus_operators","text":"kraus_operators(Π::AbstractPOVM) :: Array{Array{Complex{Float64},2},1}\n\nReturns the Kraus operators for the provided POVM. In general, the Kraus operators form a continuum and are non-unique. In this method, the construction\n\nk_i = sqrtPi_i)otimes irangle\n\n\n\n\n\n","category":"function"},{"location":"Observables/#QBase.Observables.measurement_probabilities","page":"Observables","title":"QBase.Observables.measurement_probabilities","text":"measurement_probabilities(Π::AbstractPOVM, ρ_set::Array{<:AbstractDensityMatrix}) :: QMath.Conditionals\n\nReturns the conditional probability matrix for the set of states and povm. The Probability for a POVM elemtn Π_el and state ρ are determined using the Born rule, tr(Π_el*ρ).\n\n\n\n\n\n","category":"function"},{"location":"Observables/#QBase.Observables.naimark_dilation","page":"Observables","title":"QBase.Observables.naimark_dilation","text":"naimark_dilation( Π::AbstractPOVM )\n\nReturns the dilated projectors which are equivalent to the provided POVM. During measurement, the measured state must be tensored with the ancilla.\n\n\n\n\n\n","category":"function"}]
}
