@doc raw"""
    filter_invariant_symmetries!(symmetry_group :: Symmetries, vector :: AbstractVector)

This subroutine filters the symmetries in `symmetry_group` that does 
not leave the vector `vector` invariant under their transformation.

For example, inversion symmetry is not compatible with a vector that is not the null vector,
so it will be removed from the symmetry group.
While a reflection symmetry is compatible with any vector lying on the reflection plane, 
so it will be kept in the symmetry group.

This is userful if the symmetries are evaluated from a set of atomic positions,
But we then need symmetries that also are invariant under an external perturbation along a certain direction.
"""
function filter_invariant_symmetries!(symmetry_group :: Symmetries, vector :: AbstractVector)
    tmp_vector = similar(vector)
    nat = length(vector) ÷ symmetry_group.dimension

    @assert length(vector) % symmetry_group.dimension == 0 "Error, length(vector) must be a multiple of the dimension ($(symmetry_group.dimension))."
    @assert nat == 1 || nat == get_n_atoms(symmetry_group) "Error, wrong number of atoms in the provided vector: ($(symmetry_group.dimension))."

    irt = zeros(Int, nat)
    irt[1] = 1

    for i in length(symmetry_group) :-1: 1
        if nat > 1
            irt .= symmetry_group.irt[i]
        end
        
        # Check the invariance
        tmp_vector .= 0
        apply_sym_centroid!(tmp_vector, vector, symmetry_group.symmetries[i], symmetry_group.dimension, irt)
        # println("Symmetry $i: ", symmetry_group.symmetries[i])
        # println("Vector: ", vector)
        # println("Transformed vector: ", tmp_vector)

        if norm(tmp_vector - vector) > 1e-6
            deleteat!(symmetry_group.symmetries, i)
            deleteat!(symmetry_group.irt, i)
        end
    end

    # Update the symmetry group
    update_symmetry_functions!(symmetry_group)
end
