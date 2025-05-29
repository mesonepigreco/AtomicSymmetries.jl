@doc raw"""
    apply_symmetry_vectorq!(target_vector :: AbstractMatrix{Complex{T}}, original_vector :: AbstractMatrix{Complex{T}}, symmetry_operation :: AbstractMatrix{U}, irt :: Vector{Int}, q_points :: AbstractMatrix{T})


Apply the symmetry on the original vector in q space

## Parameters

- `target_vector` : The result (modified inplace) (3n x nq)
- `original_vector` : The original vector (3n x nq)
- `symmetry_operation` : The 3x3 symmetry 
- `irt` : The atom-atom association by symmetry
- `irt_q` : The q-q association by symmetry
"""
function apply_symmetry_vectorq!(target_vector :: AbstractMatrix{Complex{T}}, original_vector :: AbstractMatrix{Complex{T}},
        symmetry_operation :: AbstractMatrix{U}, irt :: AbstractVector{Int}, irt_q :: AbstractVector{Int})
    # Apply symmetries 
    nq = length(irt_q)
    n_dims = size(symmetry_operation, 1)
    n_atoms = size(target_vector, 1) ÷ n_dims

    for iq in 1:nq
        jq = irt_q[iq]
        for i in 1:n_atoms
            j = irt[i]

            @views mul!(target_vector[n_dims * (j - 1) + 1 : n_dims * j, jq], 
                        symmetry_operation, 
                        original_vector[n_dims * (i - 1) +1 : n_dims * i, iq],
                        T(1.0), T(1.0))
        end
    end
end

@doc raw"""
    apply_symmetry_matrixq!(target_matrix :: AbstractArray{Complex{T}, 3}, original_matrix :: AbstractArray{Complex{T}, 3}, symmetry_operation :: AbstractMatrix{U}, irt :: Vector{Int}, q_points :: AbstractMatrix{T}; buffer = default_buffer())


Apply the symmetry on the matrix in q space
This subroutine assumes the convention that the phase factor is for each supercell, not atoms.
In other words, all the atoms coordinates are computed from the same origin of the supercell they are associated with.

## Parameters

- `target_vector` : The result (modified inplace) (3n x nq)
- `original_vector` : The original vector (3n x nq)
- `symmetry_operation` : The 3x3 symmetry 
- `irt` : The atom-atom association by symmetry
- `irt_q` : The q-q association by symmetry
- `buffer` : The Bumper.jl buffer for caching memory allocations [Optional]
"""
function apply_symmetry_matrixq!(target_matrix :: AbstractArray{Complex{T}, 3},
        original_matrix :: AbstractArray{Complex{T}, 3},
        sym :: AbstractMatrix{U},
        irt :: AbstractVector{Int},
        irt_q :: AbstractVector{Int}; buffer = default_buffer())

    nq = size(target_matrix, 3)
    n_dims = size(sym, 1)
    n_atoms = size(target_matrix, 1) ÷ n_dims

    @no_escape buffer begin
        work = @alloc(T, n_dims, n_dims)
        for iq in 1:nq
            iq_s = irt_q[iq]
            for i ∈ 1:n_atoms
                i_s = irt[i]
                for j in 1:n_atoms 
                    j_s = irt[j]
                    @views mul!(work, 
                                original_matrix[n_dims*(i_s-1) + 1: n_dims*i_s, n_dims*(j_s-1)+1 : n_dims*j_s, iq_s], 
                                sym, T(1.0), T(0.0))
                    @views mul!(result[n_dims*(i-1) + 1: n_dims*i, n_dims*(j - 1) + 1: n_dims*j, iq], 
                        sym', work, 1.0, 1.0)
                end
            end
        nothing
    end
end

@doc raw"""
    get_transformed_q!(irt_q :: AbstractVector{Int}, q_points :: AbstractVector{T}, sym_mat :: AbstractMatrix)

Assume that the q vectors and the symmetry operation are performed in crystal space
"""
function get_transformed_q!(irt_q :: AbstractVector{Int}, q_points :: AbstractMatrix{T}, sym_mat :: AbstractMatrix; buffer = default_buffer())
    nq = size(q_points, 2)
    ndims = size(q_points, 1)
    @no_escape buffer begin
        tmpvector = @alloc(ndims)
        tmp2 = @alloc(ndims)
        distance = @alloc(nq)
        
        for i in 1:nq
            @views mul!(tmpvector, sym_mat, q_points[:, i])

            # Check the closest q point
            min_distance = T(Inf)
            min_index = 0
            for j in 1:nq
                @views tmp2 = q_points[:, i] - tmpvector
                tmp2 .-= floor.(tmp2)
                distance = sum(abs2, tmp2) 
                if distance < min_distance
                    min_index = j
                    min_distance = distance
                end
            end

            # TODO: Check if it is this or the opposite
            irt_q[i] = min_index
        end
        nothing
    end
end
