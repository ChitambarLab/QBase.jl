using Combinatorics: permutations

# """
# stirling_2(n, k):
#
#     Counts the number of ways to partition n items into k unlabelled groups.
#     This quantity is known as Stirling's number of the 2nd kind.
#
# Inputs:
#     n: Integer, n > 1
#     k: Integer, n ≧ k ≧ 1
#
# Outputs:
#     Integer, stirling number of the 2nd kind.
# """
function stirling2(n, k)
    sum(map( i -> (-1)^i*binomial(k, i)*(k-i)^n, 0:k) )/factorial(k)
end


# """
# permutation_maps:
#   Generates the set of all square matrices which permute the elements of a vector.
#
# Input:
#   dim: Integer, length of square matrix
#
# Output:
#   maps: Array, the complete set of permutation matrices
#
# """
function permutation_maps(dim)
    ids = 1:dim

    id_permutations = collect( permutations(ids) )

    maps = []
    for perm in id_permutations
        map = zeros(Int64, (dim,dim))
        for id in ids
            map[perm[id], id] = 1
        end

        push!(maps, map)
    end

    maps
end

# """
# base_n_val(num_array, base; endian="big"):
#
#     Given an array representing a number in base-n returns the value of that
#     number in base-10.
#
# Inputs:
#     num_array: Array, contains semi-positive integers less than base.
#     base: Integer, the base-number which is represented by num_array.
#     endian: String, "big" if most significant place is at index 1 else "little"
#
# Output:
#     val: Integer, the base-10 representation the base-n number
# """
function base_n_val(num_array, base; endian="big")

    if maximum(num_array) >= base
        throw(ArgumentError("max value of $num_array is not valid in $base"))
    elseif minimum(num_array) < 0
        throw(ArgumentError("min value of $num_array must be zero or greater"))
    end

    pad = length(num_array)

    place_vals =
    if endian == "big"
        place_vals = map(
            (x) -> base^(pad - x[1]),
            enumerate(num_array)
        )
    else
        place_vals = map(
            (x) -> base^(x[1] - 1),
            enumerate(num_array)
        )
    end

    val = place_vals'*num_array

    val
end
