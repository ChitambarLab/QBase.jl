using Combinatorics: permutations, combinations

export stirling2, stirling2_partitions, stirling2_matrices
export permutations, permutation_matrices
export combinations, n_choose_k_matrices
export base_n_val

"""
    stirling2( n :: Int64, k :: Int64  ) :: Int64

Counts the number of ways to partition `n` items into `k` unlabelled groups.
This quantity is known as Stirling's number of the 2nd kind.

Throws a `DomainError` if inputs do not satisfy `n ≥ k ≥ 1`.
"""
function stirling2(n :: Int64, k :: Int64) :: Int64
    if !(n ≥ k ≥ 1)
        throw(DomainError((n,k), "Inputs (n,k) do not satisfy `n ≥ k ≥ 1`."))
    end

    sum(map( i -> (-1)^i*binomial(k, i)*(k-i)^n, 0:k) )/factorial(k)
end

"""
    stirling2_partitions( n :: Int64, k :: Int64 ) :: Vector{Vector{Vector{Int64}}}

Enumerates the unique partitions of `n` items into `k` unlabelled sets.
Each partition is a vector containing a  set of `k` vectors designating each group.

This recursive algorithm was inspired by [this blog](https://devblogs.microsoft.com/oldnewthing/20140324-00/?p=1413).
"""
function stirling2_partitions(n :: Int64, k :: Int64) :: Vector{Vector{Vector{Int64}}}
    if (n ==  0) || (k == 0)
        return Vector{Vector{Vector{Int64}}}(undef,0)
    elseif n <  k
        return Vector{Vector{Vector{Int64}}}(undef,0)
    elseif (k == 1)
        return [[[1:n...]]]
    else
        # partitions include all (n-1, k-1) paritions with [n] as k-th group
        partitions1 = _stirling2_add_partition(n, stirling2_partitions(n-1, k-1))

        # partitions include all (n-1,  k) partitions with [n] appended to each existing group
        partitions2 = _stirling2_extend_partitions(n, k, stirling2_partitions(n-1, k))

        return cat(partitions1, partitions2, dims=1)
    end
end

# helper function for adding a new group to a partition
function _stirling2_add_partition(n, partitions)
    new_partitions = []
    if length(partitions) == 0
        new_partitions = [[[n]]]
    else
        new_partitions = map(s -> cat(s, [[n]], dims=1), partitions)
    end

    new_partitions
end

# helper function for extending existing groups with a new element
function _stirling2_extend_partitions(n, k, partitions)
    new_partitions = []
    map(partition -> begin
        for i in 1:k
            push!(new_partitions, cat(partition[1:k .!= i], [cat(partition[i], [n], dims=1)], dims=1))
        end
    end, partitions)

    new_partitions
end

"""
    stirling2_matrices( n :: Int64, k :: Int64 ) :: Vector{Matrix{Int64}}

Generates the set of matrices with `k` rows and `n` columns where rows correspond
to the groups and columns are the grouped elements. A non-zero element designates
that the column id is grouped into the corresponding row.

A `DomainError` is thrown if `n ≥ k ≥ 1` is not satisfied.
"""
function stirling2_matrices(n :: Int64, k :: Int64) :: Vector{Matrix{Int64}}
    if !(n ≥ k ≥ 1)
        throw(DomainError((n,k), "Inputs (n,k) do not satisfy `n ≥ k ≥ 1`."))
    end

    map(partition -> begin
        m = zeros(Int64, k, n)
        for row_id in 1:k
            col_ids = partition[row_id]
            m[row_id,col_ids] .= 1
        end

        return m
    end, stirling2_partitions(n, k))
end

"""
    permutation_matrices( dim :: Int64 ) :: Vector{Matrix{Int64}}

Generates the set of square permutation matrices of dimension `dim`.
"""
function permutation_matrices(dim :: Int64) :: Vector{Matrix{Int64}}
    map(perm_ids -> Matrix{Int64}(I,dim,dim)[:,perm_ids], permutations(1:dim))
end

"""
    n_choose_k_matrices( n :: Int64, k :: Int64 ) :: Vector{Matrix{Int64}}

Generates a set of `n` by `k` matrices which represent all combinations of selecting
`k` columns from `n` rows.  Each column, contains a single non-zero element and
`k` rows contain a non-zero element.

A `DomainError` is thrown if `n ≥ k ≥ 1` is not satisfied.
"""
function n_choose_k_matrices(n :: Int64, k :: Int64) :: Vector{Matrix{Int64}}
    if !(n ≥ k ≥ 1)
        throw(DomainError((n,k), "Inputs (n,k) must satisfy `n ≥ k ≥ 1`."))
    end

    map(combo -> begin
        m = zeros(Int64,n,k)
        for col_id in 1:k
            m[combo[col_id],col_id] = 1
        end

        return m
    end, combinations(1:n, k))
end

"""
    base_n_val(
        num_array :: Vector{Int64},
        base :: Int64;
        big_endian=true :: Bool
    ) :: Int64

Given an array representing a number in base-n returns the value of that
number in base-10.

Inputs:
* `num_array` - Vector containing semi-positive integers less than base.
* `base` - The base-n number represented by num_array.
* `big_endian` - `true` if most significant place is at index 1, else `false`.
"""
function base_n_val(num_array :: Vector{Int64}, base :: Int64; big_endian = true :: Bool) :: Int64

    if maximum(num_array) >= base
        throw(ArgumentError("max value of $num_array is not valid in $base"))
    elseif minimum(num_array) < 0
        throw(ArgumentError("min value of $num_array must be zero or greater"))
    end

    pad = length(num_array)

    place_vals = big_endian ? map(
        x -> base^(pad - x[1]), enumerate(num_array)
    ) : map(
        x -> base^(x[1] - 1), enumerate(num_array)
    )

    val = place_vals'*num_array

    val
end
