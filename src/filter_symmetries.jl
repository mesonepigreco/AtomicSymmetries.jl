@doc raw"""
    CosetSymmetryGroup{T <: GenericSymmetries} <: GenericSymmetries

This is a way to represent a symmetry group as a subgroup times its coset representatives.
This is useful to apply symmetries to a perturbation that breaks the group.

The full original group G is obtained by combining every coset representative ``g_k`` with the subgroup ``H``.
Full symmetrization on ``G`` is equivalent to first symmetrize with respect to the cosets,
then impose the symmetries of the subgroup ``H``. For this reason,
the cosets are intended as right cosets, i.e. applying ``g_k`` to the element ``h_j`` of
the representation of ``H`` is done to the right-side.

```math
h_j g_k
```

The subgroup ``H`` is the largest subgroup of ``G`` that leaves the perturbation vector invariant.
The coset representatives ``g_k`` (including the identity) satisfy ``G = \bigcup_k H g_k``.

The constructor can be called as

```julia
CosetSymmetryGroup(original_group :: Symmetries{T}, perturbation_vector :: AbstractVector, cell :: AbstractMatrix; buffer=default_buffer())
CosetSymmetryGroup(original_group :: SymmetriesQSpace{T}, perturbation_vector :: AbstractVector, cell :: AbstractMatrix; buffer=default_buffer())
```

"""
struct CosetSymmetryGroup{T <: GenericSymmetries} <: GenericSymmetries
    subgroup :: T
    cosets :: T
end

get_dimensions(csg :: CosetSymmetryGroup) = get_dimensions(csg.subgroup)
get_n_atoms(csg :: CosetSymmetryGroup) = get_n_atoms(csg.subgroup)
Base.length(csg :: CosetSymmetryGroup) = length(csg.subgroup) * length(csg.cosets)



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
function symmetry_violation_mask!(mask :: AbstractVector{Bool}, 
                      symmetry_group :: Symmetries,
                      vector :: AbstractVector{T},
                      cell :: AbstractMatrix{T};
                      buffer=default_buffer()) where T
    n_dims = get_dimensions(symmetry_group)
    nat = length(vector) ÷ n_dims

    @assert length(vector) % n_dims == 0 "Error, length(vector) must be a multiple of the dimension ($(symmetry_group.dimension))."
    @assert nat == 1 || nat == get_n_atoms(symmetry_group) "Error, wrong number of atoms in the provided vector: ($(n_dims))."

    mask .= false


    irt = zeros(Int, nat)
    tmp_vector = zeros(T, size(vector)...)
    crystal_coordinates = zeros(T, length(vector))
    irt .= 0
    irt[1] = 1

    # Convert to crystalline coordinates
    get_crystal_coords!(reshape(crystal_coordinates, n_dims, :), 
                       reshape(vector, n_dims, :), 
                       cell; buffer=buffer)

    for i in 1:length(symmetry_group)
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

        symmetry_violation_mask!(mask, symmetry_group, vector, cell; buffer=buffer)
        has_translations = length(symmetry_group.translations) == length(mask)
        # Delete from reverse to avoid bugs
        for i in length(mask):-1:1
            if mask[i]
                deleteat!(symmetry_group.symmetries, i)
                deleteat!(symmetry_group.irt, i)
                deleteat!(symmetry_group.unit_cell_translations, i)
                if has_translations
                    deleteat!(symmetry_group.translations, i)
                end
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
        mask = @alloc(Bool, length(q_symmetry_group))

        symmetry_violation_mask!(mask, q_symmetry_group.symmetries, vector, cell; buffer=buffer)
        has_translations = length(q_symmetry_group.symmetries.translations) == length(mask)
        # Delete from reverse to avoid bugs
        for i in length(mask):-1:1
            if mask[i]
                deleteat!(q_symmetry_group.symmetries.symmetries, i)
                deleteat!(q_symmetry_group.symmetries.irt, i)
                deleteat!(q_symmetry_group.symmetries.unit_cell_translations, i)
                deleteat!(q_symmetry_group.irt_q, i)
                if has_translations
                    deleteat!(q_symmetry_group.symmetries.translations, i)
                end
            end
        end
    end
    
    # Update the symmetry group
    update_symmetry_functions!(q_symmetry_group.symmetries)
end


@doc raw"""
    CosetSymmetryGroup(original_group :: Symmetries{T}, perturbation_vector :: AbstractVector, cell :: AbstractMatrix; buffer=default_buffer()) where T

Construct a `CosetSymmetryGroup` by decomposing `original_group` into a subgroup invariant
under `perturbation_vector` and its right coset representatives.

The perturbation vector is assumed in Cartesian coordinates.
"""
function CosetSymmetryGroup(original_group :: Symmetries{T}, perturbation_vector :: AbstractVector, cell :: AbstractMatrix; buffer=default_buffer()) where T
    n_symmetries_original = length(original_group)
    n_dims = get_dimensions(original_group)
    nat = get_n_atoms(original_group)

    # Compute which symmetries are violated by the perturbation
    mask_violation = Vector{Bool}(undef, n_symmetries_original)
    symmetry_violation_mask!(mask_violation, original_group, perturbation_vector, cell; buffer=buffer)

    # Build the subgroup: deepcopy and filter out violated symmetries
    subgroup = deepcopy(original_group)
    filter_invariant_symmetries!(subgroup, perturbation_vector, cell; buffer=buffer)
    n_subgroup = length(subgroup)

    # Track which elements of G are covered by existing cosets
    covered = falses(n_symmetries_original)

    # Mark subgroup elements as covered (they belong to the trivial coset H*e)
    for i in 1:n_symmetries_original
        if !mask_violation[i]
            covered[i] = true
        end
    end

    # Build the cosets Symmetries object, starting with identity
    cosets_group = get_identity_symmetry_group(T; dims=n_dims, n_atoms=nat, translations=true)

    # Find additional coset representatives
    for i in 1:n_symmetries_original
        if covered[i]
            continue
        end

        # i is a new coset representative
        g_i = original_group.symmetries[i]
        g_i_inv = inv(g_i)

        # Add to cosets group
        add_symmetry!(cosets_group, g_i; update=false, irt=copy(original_group.irt[i]))
        if length(original_group.translations) >= i
            push!(cosets_group.translations, copy(original_group.translations[i]))
        end
        push!(cosets_group.unit_cell_translations, copy(original_group.unit_cell_translations[i]))

        # Mark all elements in the right coset H * g_i as covered
        for j in 1:n_symmetries_original
            if covered[j]
                continue
            end
            # Check if original_group[j] * g_i^{-1} ∈ H
            product = original_group.symmetries[j] * g_i_inv
            for k in 1:n_subgroup
                if maximum(abs.(product .- subgroup.symmetries[k])) < 1e-6
                    covered[j] = true
                    break
                end
            end
        end
    end

    @assert all(covered) "Not all elements of G were covered by cosets. This indicates a bug in the coset decomposition."
    @assert length(cosets_group) * n_subgroup == n_symmetries_original "Coset decomposition inconsistent: $(length(cosets_group)) cosets × $(n_subgroup) subgroup ≠ $(n_symmetries_original) original"

    update_symmetry_functions!(cosets_group)

    return CosetSymmetryGroup{Symmetries{T}}(subgroup, cosets_group)
end

@doc raw"""
    CosetSymmetryGroup(original_group :: SymmetriesQSpace{T}, perturbation_vector :: AbstractVector, cell :: AbstractMatrix; buffer=default_buffer()) where T

Construct a `CosetSymmetryGroup` for a q-space symmetry group. The subgroup and cosets
are wrapped as `SymmetriesQSpace` with properly initialized `irt_q`.
"""
function CosetSymmetryGroup(original_group :: SymmetriesQSpace{T}, perturbation_vector :: AbstractVector, cell :: AbstractMatrix; buffer=default_buffer()) where T
    # Build the real-space coset decomposition
    real_coset = CosetSymmetryGroup(original_group.symmetries, perturbation_vector, cell; buffer=buffer)

    # Wrap both subgroup and cosets in SymmetriesQSpace to compute irt_q
    q_points = original_group.q_points
    subgroup_q = SymmetriesQSpace(real_coset.subgroup, q_points; buffer=buffer)
    cosets_q = SymmetriesQSpace(real_coset.cosets, q_points; buffer=buffer)

    return CosetSymmetryGroup{SymmetriesQSpace{T}}(subgroup_q, cosets_q)

end




