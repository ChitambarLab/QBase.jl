export partial_trace, computational_basis_vectors

# validation methods
export commutes, is_hermitian, is_square, is_positive_Semidefinite

"""
    partial_trace(
        ρ::AbstractMatrix, subsystem_dims::Vector{Int64}, subsystem_id::Int64
    ) :: AbstractMatrix

Performs the partial trace on matrix `ρ` with respect to the subsystem at the specified `subsystem_id`.
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
        and have dimension equal to the product of the `subsystem_dims`.
* `subsystem_dims` :  A vector containing positive integer elements which specify the
        dimension of each subsystem. *E.g.* A 3-qubit system has `subsystem_dims = [2,2,2]`,
        a system containing a qubit and a qutrit has `subsystem_dims = [2,3]`.
* `subsystem_id` : A positive integer indexing the `subsystem_dims` vector to signify
        the subsytem which the partial trace will trace out.

*Output:*
* The reduced matrix which results after the partial trace is performed on `ρ`. This
        matrix is square and has dimension of product of `subsystem_dims` with the element
        corresponding to `subsystem_id` removed.

A `DomainError` is thrown if `ρ` is not square or of proper dimension.
"""
function partial_trace(ρ::AbstractMatrix, subsystem_dims::Vector{Int64}, subsystem_id::Int64) :: AbstractMatrix
    system_dim = .*(subsystem_dims...)

    if size(ρ) != (system_dim, system_dim)
        throw(DomainError(ρ, "Matrix ρ should be square and have dimension equal to the product of subsystem dimensions."))
    end

    trace_dim = subsystem_dims[subsystem_id]

    pre_dim = (subsystem_id > 1) ? .*(subsystem_dims[1:(subsystem_id-1)]...) : 1
    pre_id = Matrix(1I, pre_dim, pre_dim)

    post_dim = (subsystem_id < length(subsystem_dims)) ? .*(subsystem_dims[(subsystem_id + 1):end]...) : 1
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
    identity = Matrix(1I, dim, dim)
    map(i -> identity[:,i], 1:dim)
end

"""
    is_hermitian( matrix :: Matrix ) :: Bool

Returns `true` if the supplied matrix is hermitian (self-adjoint).
"""
function is_hermitian(matrix :: Matrix) :: Bool
    isapprox(matrix, adjoint(matrix), atol=10e-6)
end


"""
    is_square( matrix :: Matrix ) :: Bool

Returns `true` if the provided matrix is square.
"""
function is_square(matrix :: Matrix) :: Bool
    dims = size(matrix)
    (dims[1] == dims[2])
end

"""
    is_positive_semidefinite(matrix :: Matrix) :: Bool

Returns `true` if the supplied matrix is positive semidefinite (all eigen values
are real and greater than or equal to 0).

A `DomainError` is thrown if the provided matrix is not square.
"""
function is_positive_semidefinite(matrix :: Matrix) :: Bool
    if !(QMath.is_square(matrix))
        throw(DomainError(matrix, "Provided matrix must be square"))
    end

    evals = eigvals(matrix)

    is_pos_sd = true
    for eval in evals
        if !isapprox(imag(eval), 0, atol=10e-4)
            is_pos_sd = false
            break
        elseif (real(eval) < 0) && !isapprox(real(eval), 0, atol=10e-4)
            is_pos_sd = false
            break
        end
    end

    is_pos_sd
end

"""
    commutes(A :: Matrix, B :: Matrix) :: Bool

Returns true if matrix `A` commutes with matrix `B`. Two matrices commute if
`A*B == B*A`. The matrices must have compatible dimensions:
* A: Matrix, mxn dimensions
* B: Matrix, nxm dimensions

A `DomainError` is thrown if the matrices are not compatible.
"""
function commutes(A :: Matrix, B :: Matrix) :: Bool
    if (size(A)[1] != size(B)[2]) || (size(A)[2] != size(B)[1])
        throw(DomainError((A,B), "size(A) is not compatible with size(B)."))
    end

    all(isapprox.(A*B, B*A, atol=1e-7))
end

# """
# block_diagonal:
#   Constructs a block diagonal matrix from an array of matrices.
#
# Input:
#   matrices: Array of matrices to be block diagonalized
#
# Output:
#   block_matrix: Matrix, the block diagonalized form of matrices
# """
function block_diagonal(matrices)
    num_matrices = size(matrices)[1]

    block_matrix = []
    for row_id in 1:num_matrices

        row = []
        for col_id in 1:num_matrices

            if row_id == col_id
                new_matrix = matrices[col_id]
            else
                zero_cols = size(matrices[col_id])[2]
                zero_rows = size(matrices[row_id])[1]

                new_matrix = zeros(Int64,(zero_rows,zero_cols))
            end

            if col_id == 1
                row = new_matrix
            else
                row = cat(row, new_matrix, dims=2)
            end
        end

        if row_id == 1
            block_matrix = row
        else
            block_matrix = cat(block_matrix, row, dims=1)
        end
    end

    block_matrix
end
