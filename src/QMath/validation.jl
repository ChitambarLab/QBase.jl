# """
# is_hermitian(matrix)
#     Checks if the supplied matrix is hermitian (self-adjoint).
#
# Input:
#     matrix: A matrix to check
#
# Output:
#     Boolean, true if the matrix is hermitian.
# """
function is_hermitian(matrix)
    isapprox(matrix, adjoint(matrix), atol=10e-6)
end

# """
# is_sqaure_matrix(matrix)
#
#     Returns true if the provided matrix is square and 2-dimensional.
# """
function is_square_matrix(matrix)
    dims = size(matrix)
    (dims[1] == dims[2]) & (length(dims) == 2)
end

# """
# is_positive_semidefinite(matrix)
#     Checks if the supplied matrix is positive semoidefinite (all eigen values
#     are real and greater than or equal to 0).
#
# Input:
#     matrix: Matrix to check for validity
#
# Output:
#     is_pos_sd: Boolean, true if matrix is positive semi-definite
# """
function is_positive_semidefinite(matrix)

    if !(QMath.is_square_matrix(matrix))
        throw(ArgumentError("Provided matrix must be square"))
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
