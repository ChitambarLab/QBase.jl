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


# """
# id(dim)
#     Returns the identity matrix for the specified dimension.
# """
function id(dim)
    diagm(0 => fill(1,dim))
end

# """
# partial_trace(ρ, subsystem_dims, index)
#     Traces out the subsystem at the specified index.
#
# Inputs:
#     ρ: Matrix, a square density matrix
#     subsystem_dims: Tuple, contains positive integer elements which specify the
#         dimension of each subsystem.
#     index: Integer, the index of the subsystem list to trace out
#
# Returns:
#     ptrace: Matrix, the partial trace of the specified system.
# """
function partial_trace(ρ, subsystem_dims, index)

    trace_dim = subsystem_dims[index]
    trace_id = id(trace_dim)

    pre_dim = 1
    for sub_dim in subsystem_dims[1:(index-1)]
        pre_dim = pre_dim * sub_dim
    end

    pre_id = id(pre_dim)

    post_dim = 1
    for sub_dim in subsystem_dims[(index + 1):end]
        post_dim = post_dim * sub_dim
    end

    post_id = id(post_dim)

    partial_dim = post_dim*pre_dim
    ptrace = zeros((partial_dim,partial_dim))
    for i in 1:trace_dim

        row_multiplier = kron(pre_id, kron(trace_id[i,:]', post_id))
        col_multiplier = kron(pre_id, kron(trace_id[:,i], post_id))

        ptrace += row_multiplier * ρ * col_multiplier
    end

    ptrace
end


# """
# commutes(A,B):
#
#     Returns true if matrix A commutes with matrix B.
#
# Inputs:
#     A: Matrix, mxn dimensions
#     B: Matrix, nxm dimensions
#
# Output:
#     commutes: Boolean, true if A*B = B*A
# """
function commutes(A, B)
    if (size(A)[1] != size(B)[2]) || (size(A)[2] != size(B)[1])
        throw(ArgumentError("size(A) not compatible with size(B)"))
    end

    all(isapprox.(A*B, B*A, atol=1e-7))
end
