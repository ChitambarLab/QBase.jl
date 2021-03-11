export partial_trace, computational_basis_vectors

export n_product_id

# validation methods
export commutes, is_hermitian, is_positive_semidefinite

"""
    partial_trace(
        ρ :: Matrix, system :: Vector{Int64}, id :: Int64
    ) :: Matrix

Performs the partial trace on matrix `ρ` with respect to the subsystem at the specified `id`.
The partial trace is a quantum operation ``\\mathcal{E}`` mapping the input
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

*Inputs:*
* `ρ` : The matrix on which the partial trace is performed. This matrix should be square
        and have dimension equal to the product of the `system`.
* `system` :  A vector containing positive integer elements which specify the
        dimension of each subsystem. *E.g.* A 3-qubit system has `system = [2,2,2]`,
        a system containing a qubit and a qutrit has `subsystem_dims = [2,3]`.
* `id` : A positive integer indexing the `subsystem_dims` vector to signify
        the subsytem which the partial trace will trace out.

*Output:*
* The reduced matrix which results after the partial trace is performed on `ρ`. This
        matrix is square and has dimension of product of `subsystem_dims` with the element
        corresponding to `id` removed.

A `DomainError` is thrown if `ρ` is not square or of proper dimension.
"""
function partial_trace(ρ::Matrix, system::Vector{Int64}, id::Int64) :: Matrix
    dim = *(system...)

    if size(ρ) != (dim, dim)
        throw(DomainError(ρ, "Matrix ρ should be square and have dimension equal to the product of subsystem dimensions."))
    end

    trace_dim = system[id]

    pre_dim = (id > 1) ? *(system[1:(id-1)]...) : 1
    pre_id = Matrix(1I, pre_dim, pre_dim)

    post_dim = (id < length(system)) ? *(system[(id + 1):end]...) : 1
    post_id = Matrix(1I, post_dim, post_dim)

    basis_vectors = computational_basis_vectors(trace_dim)
    kraus_operators = map(bv -> kron(pre_id, bv', post_id), basis_vectors)

    sum(map( E -> E*ρ*E', kraus_operators))
end

"""
    computational_basis_vectors( dim::Int64 ) ::Vector{ Vector{Int64} }

Returns an orthonormal set of column vectors spanning a vector space of dimension `dim`.
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
    is_hermitian( matrix :: Matrix; atol=ATOL :: Float64 ) :: Bool

Returns `true` if the supplied matrix is hermitian (self-adjoint).
"""
function is_hermitian(A :: Matrix; atol=ATOL :: Float64) :: Bool
    indsm, indsn = axes(A)
    if indsm != indsn
        return false
    end

    for i = indsn, j = i:last(indsn)
        if A[i,j] != adjoint(A[j,i]) && !isapprox(A[i,j],A[j,i],atol=atol)
            return false
        end
    end

    return true
end

"""
    is_positive_semidefinite(matrix :: Matrix) :: Bool

Returns `true` if the supplied matrix is square and positive semidefinite (all eigen values
are real and greater than or equal to 0).
"""
function is_positive_semidefinite(A :: Matrix; atol=ATOL) :: Bool
    if !isequal(size(A)...)
        return false
    end

    for λ in eigvals(A)
        if (real(λ) < 0) && !isapprox(real(λ), 0, atol=atol)
            return false
        elseif !isapprox(imag(λ), 0, atol=atol)
            return false
        end
    end

    return true
end

"""
    n_product_id( system_ids :: Vector{Int64}, system :: Vector{Int64} ) :: Int64

For a given tensor product `system` and coordingate `system_ids` returns the corresponding
index of the tensor product.
The tenosr product `system` is a vector of positive definite integers specifying
the dimension of each of the matrices in the tensor product.
The `system_ids` is a vector of positive definite integers specifying the target
index in each matrix in the tensor product.
The `system` and `system_ids` should describe either row or column indices of the
tensor product.

A `DomainError` is thrown if any index in `system_ids` is not valid or the number of
elements in `system_ids` does not equal that of `system`.
"""
function n_product_id(system_ids :: Vector{Int64}, system::Vector{Int64}) :: Int64
    num_ids = length(system_ids)
    if num_ids != length(system)
        throw(DomainError(system_ids, "id does not contain the right number of elements"))
    elseif !all(x -> system[x] ≥ system_ids[x] ≥ 1, 1:num_ids)
        throw(DomainError(system, "not all num inputs are valid"))
    end

    id = 1
    multiplier = 1
    for i in num_ids:-1:1
        multiplier *= (i+1 > num_ids) ? 1 : system[i+1]
        id += (system_ids[i] - 1) * multiplier
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
function commutes(A :: Matrix, B :: Matrix; atol=ATOL :: Float64) :: Bool
    if (size(A)[1] != size(B)[2]) || (size(A)[2] != size(B)[1])
        throw(DomainError((A,B), "size(A) is not compatible with size(B)."))
    end

    all(isapprox.(A*B, B*A, atol=atol))
end
