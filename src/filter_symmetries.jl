@doc raw"""
    symmetry_violation_mask!(mask :: Vector{Bool}, 
                          symmetry_group :: Symmetries,
                          vector :: AbstractVector,
                          cell :: AbstractMatrix;
                          buffer=default_buffer())


Identify which are the indices of the symmetries that 
a perturbation with vector given by `vector` is violating.

## Parameters

- `mask` : The vector of Bool. True if the `vector` violates the corresponding symmetry, false otherwise
- `symmetry_group` : The group with all the symmetry operations
- `vector` : A vector in cartesian space whose invariance under symmetry operation is checked.
- `cell` : The primitive cell of the structure. Required to correctly apply symmetries in cartesian space.
- `buffer` : Optional, Bumper.jl buffer for stack allocation.
"""
function symmetry_violation_mask!(mask :: Vector{Bool}, 
                      symmetry_group :: Symmetries,
                      vector :: AbstractVector,
                      cell :: AbstractMatrix;
                      buffer=default_buffer())
    n_dims = get_dimensions(symmetry_group)
    nat = length(vector) ÷ n_dims

    @assert length(vector) % n_dims == 0 "Error, length(vector) must be a multiple of the dimension ($(symmetry_group.dimension))."
    @assert nat == 1 || nat == get_n_atoms(symmetry_group) "Error, wrong number of atoms in the provided vector: ($(n_dims))."

    mask .= false

    @no_escape buffer begin
        irt = @alloc(Int, nat)
        tmp_vector = @alloc(T, size(vector)...)
        crystal_coordinates = @alloc(T, length(vector))
        irt .= 0
        irt[1] = 1

        # Convert to crystalline coordinates
        get_crystal_coords!(reshape(crystal_coordinates, n_dims, :), 
                           reshape(vector, n_dims, :), 
                           cell; buffer=buffer)

        for i in length(symmetry_group) :-1: 1
            if nat > 1
                irt .= symmetry_group.irt[i]
            end
            
            # Check the invariance
            tmp_vector .= 0
            apply_sym_centroid!(tmp_vector, crystal_coordinates, symmetry_group.symmetries[i], symmetry_group.dimension, irt)
            # println("Symmetry $i: ", symmetry_group.symmetries[i])
            # println("Vector: ", vector)
            # println("Transformed vector: ", tmp_vector)

            if norm(tmp_vector - crystal_coordinates) > 1e-6
                mask[i] = true
            end
        end
    end
end


@doc raw"""
    filter_invariant_symmetries!(symmetry_group :: Symmetries, vector :: AbstractVector, cell :: AbstractMatrix;
                                 buffer=default_buffer())

This subroutine filters the symmetries in `symmetry_group` that does 
not leave the vector `vector` invariant under their transformation.

Note that `vector` is assumed in Cartesian coordinates (not crystalline) since version 0.6.0

For example, inversion symmetry is not compatible with a vector that is not the null vector,
so it will be removed from the symmetry group.
While a reflection symmetry is compatible with any vector lying on the reflection plane, 
so it will be kept in the symmetry group.

This is userful if the symmetries are evaluated from a set of atomic positions,
But we then need symmetries that also are invariant under an external perturbation along a certain direction.
"""
function filter_invariant_symmetries!(symmetry_group :: Symmetries, vector :: AbstractVector{T}, cell :: AbstractMatrix; buffer=default_buffer()) where T


    @no_escape buffer begin
        mask = @alloc(Bool, length(symmetry_group))

        symmetry_mask_violation!(mask, symmetry_group, vector, cell; buffer=buffer)
        # Delete from reverse to avoid bugs
        for i in length(mask):-1:1
            if mask[i]
                deleteat!(symmetry_group.symmetries, i)
                deleteat!(symmetry_group.irt, i)
            end
        end
    end
    
    # Update the symmetry group
    update_symmetry_functions!(symmetry_group)
end

@doc raw"""
    filter_invariant_symmetries!(symmetry_group :: SymmetriesQSpace, vector :: AbstractVector, cell :: AbstractMatrix;
                                 buffer=default_buffer())

This subroutine filters the symmetries in `symmetry_group` that does 
not leave the vector `vector` invariant under their transformation.

It removes the symmetries violated by the vector.
This works for symmetry groups expressed in q space

"""
function filter_invariant_symmetries!(q_symmetry_group :: SymmetriesQSpace, vector :: AbstractVector{T}, cell :: AbstractMatrix; buffer=default_buffer()) where T

    @no_escape buffer begin
        mask = @alloc(Bool, length(symmetry_group))

        symmetry_mask_violation!(mask, q_symmetry_group.symmetries, vector, cell; buffer=buffer)
        # Delete from reverse to avoid bugs
        for i in length(mask):-1:1
            if mask[i]
                deleteat!(q_symmetry_group.symmetries.symmetries, i)
                deleteat!(q_symmetry_group.symmetries.irt, i)
                deleteat!(q_symmetry_group.irt_q, i)
                deleteat!(q_symmetry_group.minus_q_index, i)
            end
        end
    end
    
    # Update the symmetry group
    update_symmetry_functions!(q_symmetry_group.symmetries)
end

