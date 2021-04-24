export partial_trace, computational_basis_vectors

export n_product_id

# validation methods
export commutes, is_hermitian, is_positive_semidefinite, is_orthonormal_basis, is_complete

"""
    partial_trace(
        ρ :: AbstractMatrix, dims :: Vector{Int64}, id :: Int64
    ) :: Matrix

Performs the partial trace on matrix `ρ` with respect to the subsystem at the specified `id`
and returns resulting reduced matrix.
The returned matrix is square with dimension equal to the product of `dims` with the element
corresponding to `id` removed.

*Inputs:*
* `ρ` : The matrix on which the partial trace is performed. This matrix should be square
        and have dimension equal to the product of the `dims`.
* `dims` :  A vector containing positive integer elements which specify the
        dimension of each subsystem. *E.g.* A 3-qubit system has `dims = [2,2,2]`,
        a system containing a qubit and a qutrit has `dims = [2,3]`.
* `id` : A positive integer indexing the `dims` vector to signify
        the subsystem to be traced out.

The partial trace can be understood as a quantum operation ``\\mathcal{E}`` mapping the input
Hilbert space to the output Hilbert space,
``\\mathcal{E} \\; : \\; \\mathcal{H}_X \\rightarrow \\mathcal{H}_Y``. A quantum
operation admits the operator sum representation, ``\\mathcal{E}(\\rho) = \\sum_i E_i  \\rho E_i^{\\dagger}``.
For example, given a 3-qubit system with density matrix ``\\rho_{ABC}``, the partial trace with
respect to system ``B``, is explicitly,

```math
\\rho_{AC} = \\text{Tr}_B[\\rho_{ABC}] = \\mathcal{E}_B(\\rho_{ABC}) = \\sum_i E_i \\rho_{ABC} E_i^{\\dagger}
```

where ``E_i = \\mathbb{I}_A \\otimes \\langle i |_B \\otimes\\mathbb{I}_C``, represents
the Kraus operators for the quantum operation and ``\\rho_{AC}`` is the reduced density
operator remaining after system ``B`` is traced out.

A `DomainError` is thrown if `ρ` is not square or of proper dimension.
"""
function partial_trace(ρ::AbstractMatrix, dims::Vector{Int64}, id::Int64) :: Matrix
    dim = *(dims...)

    if size(ρ) != (dim, dim)
        throw(DomainError(ρ, "Matrix ρ should be square and have dimension equal to the product of subsystem dimensions."))
    end

    trace_dim = dims[id]

    pre_dim = (id > 1) ? *(dims[1:(id-1)]...) : 1
    pre_id = Matrix(1I, pre_dim, pre_dim)

    post_dim = (id < length(dims)) ? *(dims[(id + 1):end]...) : 1
    post_id = Matrix(1I, post_dim, post_dim)

    # TODO: This method can be sped up significantly by avoiding kraus operator construction
    # and coputing the partial trace by a more direct means.
    basis_vectors = computational_basis_vectors(trace_dim)
    kraus_operators = map(bv -> kron(pre_id, bv', post_id), basis_vectors)

    sum(map( E -> E*ρ*E', kraus_operators))
end

"""
    computational_basis_vectors( dim::Int64 ) ::Vector{ Vector{Int64} }

Returns the set of orthonormal column vectors ``\\{|1\\rangle,\\dots,|i\\rangle,\\dots,|d\\rangle \\}``
spanning a vector space of dimension `dim`.
Note that ``|1\\rangle = ( 1, 0, \\dots, 0 )^T``.
"""
function computational_basis_vectors(dim::Int64) :: Vector{Vector{Int64}}
    _basis_vec(i) = begin
        vec = zeros(Int64, dim)
        vec[i] = 1
        return vec
    end

    map(_basis_vec, 1:dim)
end

"""
    is_hermitian( M :: AbstractMatrix{<:Number}; atol=ATOL :: Float64 ) :: Bool

Returns `true` if the supplied matrix is hermitian (self-adjoint), `M == M'`.
"""
function is_hermitian(M :: AbstractMatrix{<:Number}; atol=ATOL :: Float64) :: Bool
    indsm, indsn = axes(M)
    if indsm != indsn
        return false
    end

    for i = indsn, j = i:last(indsn)
        if M[i,j] != adjoint(M[j,i]) && !isapprox(M[i,j],M[j,i],atol=atol)
            return false
        end
    end

    return true
end

"""
    is_positive_semidefinite(M :: AbstractMatrix) :: Bool

Returns `true` if the supplied matrix is square and positive semidefinite (all eigen values
are real and greater than or equal to 0).
"""
function is_positive_semidefinite(M :: AbstractMatrix{<:Number}; atol=ATOL) :: Bool
    if !isequal(size(M)...)
        return false
    end

    for λ in eigvals(M)
        if (real(λ) < 0) && !isapprox(real(λ), 0, atol=atol)
            return false
        elseif !isapprox(imag(λ), 0, atol=atol)
            return false
        end
    end

    return true
end

"""
    is_orthonormal_basis(basis_vecs ::  Vector;  atol=ATOL ::  Float64) :: Bool

Returns `true` if `basis_vecs` ``\\{|\\psi_j\\rangle\\}_{j=1}^n forms an orthonormal basis:

```math
\\langle \\psi_j |\\psi_j \\rangle = 1 \\quad \\text{and} \\quad \\langle \\psi_j | \\psi_{k} \\rangle = 0
```

for all ``j`` and ``k`` where ``j\\neq k``.
"""
function is_orthonormal_basis(basis_vecs :: Vector{<:AbstractVector}; atol=ATOL :: Float64) :: Bool
    num_vecs = length(basis_vecs)
    if num_vecs != length(basis_vecs[1])
        return false
    end

    for i = 1:num_vecs, j = i:num_vecs
        inner_prod = basis_vecs[j]'basis_vecs[i]
        if i == j &&  !isapprox(inner_prod,1,atol=atol)
            return false
        elseif j > i  && !isapprox(inner_prod,0,atol=atol)
            return false
        end
    end
    return true
end

"""
    n_product_id( tensor_coord :: Vector{Int64}, dims :: Vector{Int64} ) :: Int64

For a given tensor product of subspace dimension `dims` and tensor coordingate `tensor_coord`
returns the corresponding matrix index of the kronecker product.
The tenosr product `dims` is a vector of positive definite integers specifying
the dimension of each of the matrices in the tensor product.
The `tensor_coord` is a vector of positive definite integers specifying the target
index in each matrix in the tensor product.
The `dims` and `tensor_coord` should describe either row or column indices of the
tensor product.

A `DomainError` is thrown if any index in `tensor_coord` is not valid or the number of
elements in `tensor_coord` does not equal that of `dims`.
"""
function n_product_id(tensor_coord::Vector{Int64}, dims::Vector{Int64}) :: Int64
    num_ids = length(tensor_coord)
    if num_ids != length(dims)
        throw(DomainError(tensor_coord, "id does not contain the right number of elements"))
    elseif !all(x -> dims[x] ≥ tensor_coord[x] ≥ 1, 1:num_ids)
        throw(DomainError(dims, "not all num inputs are valid"))
    end

    id = 1
    multiplier = 1
    for i in num_ids:-1:1
        multiplier *= (i+1 > num_ids) ? 1 : dims[i+1]
        id += (tensor_coord[i] - 1) * multiplier
    end

    id
end

"""
    commutes(A :: Matrix, B :: Matrix; atol=ATOL :: Float64) :: Bool

Returns true if matrix `A` commutes with matrix `B`. Two matrices commute if
`A*B == B*A`. The matrices must have compatible dimensions:
* A: Matrix, mxn dimensions
* B: Matrix, nxm dimensions

A `DomainError` is thrown if the matrices are not compatible.
"""
function commutes(A :: AbstractMatrix{<:Number}, B :: Matrix; atol=ATOL :: Float64) :: Bool
    if (size(A)[1] != size(B)[2]) || (size(A)[2] != size(B)[1])
        throw(DomainError((A,B), "size(A) is not compatible with size(B)."))
    end

    all(isapprox.(A*B, B*A, atol=atol))
end

"""
    is_complete(Π :: Vector{<:AbstractMatrix}; atol=ATOL :: Float64) :: Bool

Returns `true` if the provided set of matrices sums to the identity.
"""
function is_complete(Π :: Vector{<:AbstractMatrix}; atol=ATOL :: Float64) :: Bool
    if !isapprox(sum(Π), I, atol=atol)
        return false
    end

    return true
end
