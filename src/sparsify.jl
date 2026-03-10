"""
This module contains functions to convert the symmetry group into a standard
sparse matrix format. 
This is very useful for using them in standard differentiable programming,
where the inplace operations of this module are not supported.
"""

# We will override the sparse from SparseArrays


@doc raw"""
    struct SparseSymmetryGroup{T}
        symmetries :: Vector{ <: SparseMatrixCSC{T}}
    end

This struct is used to store the symmetry group of a matrix in a sparse format.
"""
struct SparseSymmetryGroup{T}
    symmetries :: Vector{ <: SparseMatrixCSC{T}}
end


@doc raw"""
    function get_sparsified_representation(symmetry_group :: Symmetries{T}) :: SparseSymmetryGroup{T} where T

This function takes a symmetry group and returns its sparse representation.
Note that the matrix are transposed for convenience in matrix multiplication.
So we have

S = symmat' where S is the sparse matrix and symmat is the original matrix.
"""
function SparseArrays.sparse(symmetry_group :: Symmetries{T}) :: SparseSymmetryGroup{T} where T
    # Get the number of symmetries
    n_symmetries = length(symmetry_group.symmetries)

    # Create a vector to store the sparse matrices
    sparse_matrices = Vector{SparseMatrixCSC{T}}(undef, n_symmetries)

    # Loop over each symmetry and convert it to a sparse matrix
    for i in 1:n_symmetries
        sparse_matrices[i] = get_sparse_symmetry(symmetry_group.symmetries[i], symmetry_group.irt[i])
    end

    return SparseSymmetryGroup(sparse_matrices)
end



function get_sparse_symmetry(symmat :: AbstractMatrix{T}, irt :: Vector{Int}; Scache=nothing) where T
    # Get a standard matrix
    natsc = length(irt)
    ndims = size(symmat, 1)
    n_modes = natsc * ndims

    # Get the matrix to be sparsified
    S = Scache

    if S === nothing
        S = zeros(T, n_modes, n_modes)
    end
    S .= 0.0

    for i ∈ 1:natsc 
        j = irt[i]
        i_ind_start = ndims * (i - 1) +  1
        j_ind_start = ndims * (j - 1) +  1

        for k in i_ind_start:i_ind_start + ndims -1
            k_orig = k - i_ind_start + 1
            for h in j_ind_start:j_ind_start + ndims  -1
                h_orig = h - j_ind_start + 1
                S[h, k] = symmat[h_orig, k_orig]
            end
        end
    end

    # Convert the matrix to sparse format
    # Since it is more convenient to apply the matrix by column, we transpose it
    sparse(S')
end

@doc raw"""
    function apply_sparse_symmetry(sparse_s :: SparseMatrixCSC{T}, v :: AbstractArray{U}) where {T, U}

This function applies the sparse symmetry matrix to a displacement vector.

The inplace version should be nonallocating.
"""
function apply_sparse_symmetry(sparse_s :: SparseMatrixCSC, v :: AbstractArray) 
    (v' * sparse_s)'
end
function apply_sparse_symmetry!(output :: AbstractArray{T}, sparse_s :: SparseMatrixCSC, v :: AbstractArray{T}; buffer = default_buffer()) where T
    @no_escape buffer begin
        w = @alloc(T, 1, length(output))
        w .= v' * sparse_s
        output .= w'
        nothing
    end
end


