@doc raw"""
    SymmetriesQSpace(symmetries :: Symmetries{T}, q_points :: AbstractMatrix{T}) :: SymmetriesQSpace{T} where T

    struct SymmetriesQSpace{T} <: GenericSymmetries where T
        symmetries :: Symmetries{T}
        irt_q :: Vector{Vector{Int}}
        minus_q_index :: Vector{Int}
    end

This structure contains the information to perform the symmetrization of a dynamical matrix directly in q space.
Note that the `q_points` needs to be in crystal coordinates.
"""
struct SymmetriesQSpace{T} <: GenericSymmetries 
    symmetries :: Symmetries{T}
    irt_q :: Vector{Vector{Int}}
    minus_q_index :: Vector{Int}
end
function SymmetriesQSpace(symmetries :: Symmetries{T}, q_points :: AbstractMatrix{T}; buffer = default_buffer()) :: SymmetriesQSpace{T} where T
    n_symmetries = length(symmetries)
    n_q = size(q_points, 2)

    irt_q = Vector{Vector{Int}}(undef, n_symmetries)
    for i in 1:n_symmetries
        irt_q[i] = Vector{Int}(undef, n_q)
        get_irt_q!(irt_q[i], q_points, symmetries.symmetries[i]; buffer = buffer)
    end

    minus_q_index = zeros(Int, n_q)
    get_minus_q!(minus_q_index, q_points)

    SymmetriesQSpace(symmetries, irt_q, minus_q_index)
end


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
        symmetry_operation :: AbstractMatrix{U}, irt :: AbstractVector{Int}, irt_q :: AbstractVector{Int}) where {T, U}
    # Apply symmetries 
    nq = length(irt_q)
    n_dims = size(symmetry_operation, 1)
    n_atoms = size(target_vector, 1) ÷ n_dims

    for iq in 1:nq
        jq = irt_q[iq]
        for i in 1:n_atoms
            j = irt[i]

            @views mul!(target_vector[n_dims * (j - 1) + 1: n_dims * j, jq], 
                        symmetry_operation, 
                        original_vector[n_dims * (i - 1) + 1: n_dims * i, iq],
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
        irt_q :: AbstractVector{Int}; buffer = default_buffer()) where {T, U}

    nq = size(target_matrix, 3)
    n_dims = size(sym, 1)
    n_atoms = size(target_matrix, 1) ÷ n_dims

    @no_escape buffer begin
        work = @alloc(Complex{T}, n_dims, n_dims)
        for iq in 1:nq
            iq_s = irt_q[iq]
            for i ∈ 1:n_atoms
                i_s = irt[i]
                for j in 1:n_atoms 
                    j_s = irt[j]
                    @views mul!(work, 
                                original_matrix[n_dims*(i_s-1) + 1: n_dims*i_s, n_dims*(j_s-1)+1 : n_dims*j_s, iq_s], 
                                sym, T(1.0), T(0.0))
                    @views mul!(target_matrix[n_dims*(i-1) + 1: n_dims*i, n_dims*(j - 1) + 1: n_dims*j, iq], 
                        sym', work, 1.0, 1.0)
                end
            end
        end
        nothing
    end
end

@doc raw"""
    get_irt_q!(irt_q :: AbstractVector{Int}, q_points :: AbstractVector{T}, sym_mat :: AbstractMatrix)

Get the correspondance q' = Sq on the provided q grid.
Always assume everything is in crystal coordinates
"""
function get_irt_q!(irt_q :: AbstractVector{Int}, q_points :: AbstractMatrix{T}, sym_mat :: AbstractMatrix; buffer = default_buffer()) where T
    nq = size(q_points, 2)
    ndims = size(q_points, 1)
    @no_escape buffer begin
        tmpvector = @alloc(T, ndims)
        tmp2 = @alloc(T, ndims)
        distance = @alloc(T, nq)
        
        for i in 1:nq
            @views mul!(tmpvector, sym_mat, q_points[:, i])

            # Check the closest q point
            min_distance = T(Inf)
            min_index = 0
            for j in 1:nq
                @views tmp2 .= q_points[:, j] - tmpvector
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

@doc raw"""
    get_minus_q!(minus_q_index :: AbstractVector{Int}, q_points :: AbstractMatrix{T}; buffer = default_buffer()) where T

Identify for each q point what is the corresponding -q:

``
\vec q \longrightarrow -\vec q + \vec G
``

where ``\vec G`` is a reciprocal vector. Since this is done in crystal coordinates``\vec G`` are all possible integers.
"""
function get_minus_q!(minus_q_index :: AbstractVector{Int}, q_points :: AbstractMatrix{T}; buffer = default_buffer()) where T
    @no_escape buffer begin
        ndims = size(q_points, 1)
        nq = size(q_points, 2)
        tmpvector = @alloc(T, ndims)
        tmp2 = @alloc(T, ndims)
        distance = @alloc(T, nq)
        
        for i in 1:nq
            @views tmpvector .= -q_points[:, i]

            # Check the closest q point
            min_distance = T(Inf)
            min_index = 0
            for j in 1:nq
                @views tmp2 .= q_points[:, j] - tmpvector
                tmp2 .-= floor.(tmp2)
                distance = sum(abs2, tmp2) 
                if distance < min_distance
                    min_index = j
                    min_distance = distance
                end
            end

            # TODO: Check if it is this or the opposite
            minus_q_index[i] = min_index
        end
        nothing
    end
end

@doc raw"""
    symmetrize_vector_q!(target_gamma :: AbstractVector{T}, original_q :: AbstractArray{Complex{T}, 2}, symmetries :: Symmetries, irt_q :: Vector{Vector{Int}}; buffer = default_buffer() where T

Impose the symmetrization of a vector in q space.
Since the symmetrization also imposes translational symmetries, the result is always a vector only at gamma.

The symmetrized vector is supposed to be a displacement (so no translations are applied)


## Parameters

- `target_gamma` : The `n_at * n_dims` output symmetrized vector at ``\Gamma``
- `original_q` : The original vector in q-space of size `nat*n_dims, nq`
- `symmetries` : The symmetry group
- `irt_q` : A vector (one for each symmetry) of the correspondances of q points. For each symmetry can be obtained from `get_irt_q!`
- `gamma_index` : Specify which q vector is ``\Gamma``. If not specified, it is assumed to be the first one
"""
function symmetrize_vector_q!(target_gamma :: AbstractVector{T}, original_q :: AbstractArray{Complex{T}, 2}, symmetries :: Symmetries, irt_q :: Vector{Vector{Int}}; gamma_index=1, buffer = default_buffer()) where {T}

    n_modes = size(original_q, 1)
    n_q = size(original_q, 2)

    @assert gamma_index <= n_q, "Error, the number of q points ($n_q) cannot be lower than the index of Γ ($gamma_index)"

    @no_escape buffer begin
        tmp_vector = @alloc(Complex{T}, n_modes, n_q)

        for i in 1:length(symmetries)
            sym_mat = sym.symmetries[i]
            irt = sym.irt[i]
            q_irt = irt_q[i]

            apply_symmetry_vectorq!(tmp_vector, original_q, sym_mat, irt)
        end

        tmp_vector ./= length(symmetries)
        target_gamma .= real(tmp_vector[:, gamma_index])
        nothing
    end
end

@doc raw"""
    symmetrize_matrix_q!(target_q :: AbstractArray{Complex{T}, 3}, original_q :: AbstractArray{Complex{T}, 3}, symmetries :: Symmetries, irt_q :: Vector{Vector{Int}}; buffer = default_buffer() where T
    symmetrize_matrix_q!(matrix_q :: AbstractArray{Complex{T}, 3}, q_symmetries :: SymmetriesQSpace; buffer = default_buffer())  where T
    symmetrize_matrix_q!(target_q :: AbstractArray{T, 3}, original_q :: AbstractArray{Complex{T}, 3}, q_symmetries :: SymmetriesQSpace; buffer = default_buffer())  where T


Impose the symmetrization of a dynamical matrix in q space.
The matrix must be in crystal coordinates.


## Parameters

- `target_q` : The symmetrized matrix of size `n_modes, n_modes, nq` (modified in-place).
- `original_q` : The original matrix in q-space of size `n_modes, n_modes, nq`. It could be the same as target_q
- `symmetries` : The symmetry group
- `irt_q` : A vector (one for each symmetry) of the correspondances of q points. For each symmetry can be obtained from `get_irt_q!`
- `minus_q_index` : A vector containing for each `q` the corresponding ``\vec {q'} = -\vec q + \vec G``, where ``\vec G`` is a generic reciprocal lattice vector.
"""
function symmetrize_matrix_q!(target_q :: AbstractArray{Complex{T}, 3}, original_q :: AbstractArray{Complex{T}, 3}, symmetries :: Symmetries, irt_q :: Vector{Vector{Int}}, minus_q_index::Vector{Int}; buffer = default_buffer())  where T

    n_modes = size(original_q, 1)
    n_q = size(original_q, 3)

    @no_escape buffer begin
        tmp_matrix = @alloc(Complex{T}, n_modes, n_modes, n_q)
        target_q .= Complex{T}(0.0)

        for i in 1:length(symmetries)
            sym_mat = symmetries.symmetries[i]
            irt = symmetries.irt[i]
            q_irt = irt_q[i]

            apply_symmetry_matrixq!(tmp_matrix, original_q, sym_mat, irt, q_irt; buffer=buffer)
            println("After symmetrization: sym $i")
            @show original_q
            @show tmp_matrix
        end

        tmp_matrix ./= length(symmetries)
        nothing

        # Apply the hermitianity
        for iq in 1:n_q
            @views tmp_matrix[:,:, iq] .+= tmp_matrix[:, :, iq]'
        end
        tmp_matrix ./= T(2)

        # Apply the time-reversal symmetry
        target_q .= tmp_matrix
        for iq in 1:n_q
            @views target_q[:, :, iq] .+= conj.(tmp_matrix[:, :, minus_q_index[iq]]')
        end
        target_q ./= T(2)
    end
end
function symmetrize_matrix_q!(target_q :: AbstractArray{Complex{T}, 3}, original_q :: AbstractArray{Complex{T}, 3}, q_symmetries :: SymmetriesQSpace; buffer = default_buffer())  where T
    symmetrize_matrix_q!(target_q, original_q, q_symmetries.symmetries, q_symmetries.irt_q, q_symmetries.minus_q_index; buffer=buffer)
end
function symmetrize_matrix_q!(matrix_q :: AbstractArray{Complex{T}, 3}, q_symmetries :: SymmetriesQSpace; buffer = default_buffer())  where T
    symmetrize_matrix_q!(matrix_q, matrix_q, q_symmetries; buffer=buffer)
end

