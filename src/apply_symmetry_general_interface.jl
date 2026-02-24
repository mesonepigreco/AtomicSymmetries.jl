@doc raw"""
    rotate_vector!(new_vector :: AbstractVector{T}, old_vector :: AbstractVector{T}, cell :: Matrix{T},
            reciprocal_vectors :: Matrix{T},
            symmetry_group :: Symmetries, sym_index :: Int; buffer=default_buffer())
    rotate_vector!(new_vector :: AbstractMatrix{Complex{T}}, old_vector :: AbstractMatrix{Complex{T}},
            cell :: Matrix{T}, reciprocal_vectors :: Matrix{T},
            symmetry_group :: SymmetriesQSpace, sym_index :: Int; buffer=default_buffer())

Apply a **single** symmetry operation on a vector (e.g. a displacement or force).
The vector must be provided in Cartesian coordinates; the conversion to crystal
coordinates and back is handled internally.

This function works both with `SymmetriesQSpace` and with standard real-space
`Symmetries`, thanks to multiple dispatch.

**Real-space version** — the vector has length ``n_\text{dims} \times n_\text{atoms}``
and the operation applies the symmetry rotation plus the atom permutation:

```math
\vec v'_{\text{irt}[a]} = S\, \vec v_{a}
```

**Q-space version** — the vector has size ``(n_\text{modes},\, n_q)`` and the
operation additionally permutes q-points according to the symmetry:

```math
\vec v'_{\text{irt}[a]}(q') = S\, \vec v_{a}(q), \qquad q' = S^{-T} q
```

To symmetrize a vector (average over all symmetries), call this function for
each symmetry and average the result. That is equivalent to what
`symmetrize_vector!` (real space) or `symmetrize_vector_cartesian_q!` (q space)
do internally.

## Parameters

- `new_vector` : Output rotated vector (modified in-place)
- `old_vector` : Input vector to be rotated
- `cell` : Primitive cell matrix (lattice vectors as columns)
- `reciprocal_vectors` : Reciprocal lattice vectors (column-wise)
- `symmetry_group` : The symmetry group (`Symmetries` or `SymmetriesQSpace`)
- `sym_index` : Index of the symmetry operation to apply (``1 \le \text{sym\_index} \le N_\text{sym}``)
- `buffer` : Optional Bumper.jl buffer for stack allocations

## Example

```julia
n_sym = get_nsymmetries(sym_group)
avg = zeros(length(vector))
rotated = zeros(length(vector))
for i in 1:n_sym
    rotated .= 0
    rotate_vector!(rotated, vector, cell, reciprocal_vectors, sym_group, i)
    avg .+= rotated
end
avg ./= n_sym  # same as symmetrize_vector!(vector, cell, sym_group)
```

## See also

- [`symmetrize_vector!`](@ref) — symmetrize a vector in real space (Cartesian)
- [`symmetrize_vector_cartesian_q!`](@ref) — symmetrize a vector in q space (Cartesian)
"""
function rotate_vector!(new_vector :: AbstractVector{T}, old_vector :: AbstractVector{T}, cell :: Matrix{T},
        reciprocal_vectors :: Matrix{T}, symmetry_group :: Symmetries, sym_index :: Int; buffer=default_buffer()) where T

    # Get the dimensions
    n_dims = get_dimensions(symmetry_group)
    n_modes = length(new_vector)
    n_atoms = n_modes ÷ n_dims
    new_vector .= zero(T)

    @no_escape buffer begin
        tmp_vector = @alloc(T, n_modes)

        # Convert to crystal
        cryst_cart_conv!(reshape(tmp_vector, n_dims, :),
                         reshape(old_vector, n_dims, :),
                         cell, reciprocal_vectors, false; q_space=false)

        # Apply the matrix
        apply_sym_centroid!(new_vector, tmp_vector, symmetry_group.symmetries[sym_index], n_dims, symmetry_group.irt[sym_index]; buffer)
        tmp_vector .= new_vector

        # Convert back to cartesian
        cryst_cart_conv!(reshape(new_vector, n_dims, :),
                         reshape(tmp_vector, n_dims, :),
                         cell, reciprocal_vectors, true; q_space=false)
        nothing
    end
end
function rotate_vector!(new_vector :: AbstractMatrix{Complex{T}}, old_vector :: AbstractMatrix{Complex{T}},
        cell :: Matrix{T}, reciprocal_vectors :: Matrix{T},
        symmetry_group :: SymmetriesQSpace, sym_index :: Int; buffer=default_buffer()) where T

    # Get the dimensions
    n_q = size(new_vector, 2)
    n_dims = get_dimensions(symmetry_group)
    n_modes = size(new_vector, 1)
    n_atoms = n_modes ÷ n_dims
    new_vector .= zero(Complex{T})

    @no_escape buffer begin
        tmp_vector = @alloc(Complex{T}, n_modes, n_q)

        # Convert to crystal
        cryst_cart_conv!(reshape(tmp_vector, n_dims, :),
                         reshape(old_vector, n_dims, :),
                         cell, reciprocal_vectors, false; q_space=false)

        # Apply the matrix
        apply_symmetry_vectorq!(new_vector, tmp_vector, symmetry_group[sym_index], symmetry_group.symmetries.irt[sym_index],
                                symmetry_group.irt_q[sym_index])
        tmp_vector .= new_vector

        # Convert back to cartesian
        cryst_cart_conv!(reshape(new_vector, n_dims, :),
                         reshape(tmp_vector, n_dims, :),
                         cell, reciprocal_vectors, true; q_space=false)
        nothing
    end
end


@doc raw"""
    rotate_matrix!(new_matrix :: AbstractMatrix{T}, old_matrix :: AbstractMatrix{T}, cell :: Matrix{T},
            reciprocal_vectors :: Matrix{T},
            symmetry_group :: Symmetries, sym_index :: Int; buffer=default_buffer()) where T
    rotate_matrix!(new_matrix :: AbstractMatrix{T}, old_matrix :: AbstractMatrix{T}, cell :: Matrix{T},
            reciprocal_vectors :: Matrix{T},
            symmetry_group :: SymmetriesQSpace, sym_index :: Int; buffer=default_buffer()) where T

Apply a **single** symmetry operation on an intensive ``n_\text{dims} \times n_\text{dims}``
matrix (e.g. a stress tensor or a dielectric tensor).
The matrix must be provided in Cartesian coordinates; the conversion to crystal
coordinates and back is handled internally.

The operation performed is a second-rank tensor rotation:

```math
M' = S^\top\, M\, S
```

where ``S`` is the symmetry rotation matrix (in crystal coordinates).

Since an intensive matrix is q-independent, both the `Symmetries` and
`SymmetriesQSpace` dispatches perform exactly the same operation
(the `SymmetriesQSpace` version delegates to the `Symmetries` version).

To symmetrize a matrix, call this function for each symmetry and average. For
a cubic crystal, the averaged stress tensor will be proportional to the
identity matrix.

## Parameters

- `new_matrix` : Output rotated matrix (modified in-place)
- `old_matrix` : Input matrix to be rotated
- `cell` : Primitive cell matrix (lattice vectors as columns)
- `reciprocal_vectors` : Reciprocal lattice vectors (column-wise)
- `symmetry_group` : The symmetry group (`Symmetries` or `SymmetriesQSpace`)
- `sym_index` : Index of the symmetry operation to apply (``1 \le \text{sym\_index} \le N_\text{sym}``)
- `buffer` : Optional Bumper.jl buffer for stack allocations

## Example

```julia
n_sym = get_nsymmetries(sym_group)
avg = zeros(3, 3)
rotated = zeros(3, 3)
for i in 1:n_sym
    rotate_matrix!(rotated, stress_tensor, cell, reciprocal_vectors, sym_group, i)
    avg .+= rotated
end
avg ./= n_sym  # symmetrized stress tensor
```

## See also

- [`rotate_dynamical_matrix!`](@ref) — for extensive ``n_\text{modes} \times n_\text{modes`` matrices with atom-block structure
- [`symmetrize_fc!`](@ref) — symmetrize a force constant matrix in real space (Cartesian)
"""
function rotate_matrix!(new_matrix :: AbstractMatrix{T}, old_matrix :: AbstractMatrix{T}, cell :: Matrix{T},
        reciprocal_vectors :: Matrix{T}, symmetry_group :: Symmetries, sym_index :: Int; buffer=default_buffer()) where T

    n_dims = get_dimensions(symmetry_group)

    @no_escape buffer begin
        tmp_matrix = @alloc(T, n_dims, n_dims)
        work = @alloc(T, n_dims, n_dims)

        # Convert to crystal: treat as a single-atom (n_atoms=1) block
        cart_cryst_matrix_conversion!(tmp_matrix, old_matrix, cell; cart_to_cryst=true, buffer=buffer)

        # Apply rotation: new = S' * old_cryst * S
        sym_mat = symmetry_group.symmetries[sym_index]
        mul!(work, tmp_matrix, sym_mat)
        mul!(new_matrix, sym_mat', work)

        # new_matrix is now in crystal coords, convert back to cartesian
        tmp_matrix .= new_matrix
        cart_cryst_matrix_conversion!(new_matrix, tmp_matrix, cell; cart_to_cryst=false, buffer=buffer)

        nothing
    end
end
function rotate_matrix!(new_matrix :: AbstractMatrix{T}, old_matrix :: AbstractMatrix{T}, cell :: Matrix{T},
        reciprocal_vectors :: Matrix{T}, symmetry_group :: SymmetriesQSpace, sym_index :: Int; buffer=default_buffer()) where T
    # Delegate to the Symmetries version — stress tensor is intensive (q-independent)
    rotate_matrix!(new_matrix, old_matrix, cell, reciprocal_vectors, symmetry_group.symmetries, sym_index; buffer=buffer)
end


@doc raw"""
    rotate_dynamical_matrix!(new_matrix :: AbstractMatrix{T}, old_matrix :: AbstractMatrix{T}, cell :: Matrix{T},
            reciprocal_vectors :: Matrix{T},
            symmetry_group :: Symmetries, sym_index :: Int; buffer=default_buffer()) where T
    rotate_dynamical_matrix!(new_matrix :: AbstractArray{Complex{T}, 3}, old_matrix :: AbstractArray{Complex{T}, 3}, cell :: Matrix{T},
            reciprocal_vectors :: Matrix{T},
            symmetry_group :: SymmetriesQSpace, sym_index :: Int; buffer=default_buffer()) where T

Apply a **single** symmetry operation on a dynamical matrix (force-constant-like
matrix with atom blocks). The matrix must be provided in Cartesian coordinates;
the conversion to crystal coordinates and back is handled internally.

**Real-space version** — the matrix has size
``(n_\text{modes}, n_\text{modes})`` where ``n_\text{modes} = n_\text{dims} \times n_\text{atoms}``.
Each ``n_\text{dims} \times n_\text{dims}`` block ``(a, b)`` is rotated by the
symmetry and the atom indices are permuted according to `irt`:

```math
\Phi'_{\text{irt}[a],\, \text{irt}[b]} = S^\top\, \Phi_{a b}\, S
```

**Q-space version** — the matrix has size ``(n_\text{modes}, n_\text{modes}, n_q)``.
In addition to the block-wise rotation and atom permutation, the q-point is
also permuted and phase factors from fractional translations are included:

```math
D'_{\text{irt}[a],\, \text{irt}[b]}(q')
= e^{2\pi i\, q \cdot (\vec t_a - \vec t_b)}\,
  S^\top\, D_{a b}(q)\, S
```

where ``q' = S^{-T} q`` and ``\vec t_a`` are the unit-cell translations
that bring the symmetry-transformed atom back into the primitive cell.

To symmetrize a dynamical matrix, call this function for each symmetry and
average the result. That is equivalent to what `symmetrize_fc!` (real space)
or `symmetrize_matrix_cartesian_q!` (q space) do internally.

## Parameters

- `new_matrix` : Output rotated matrix (modified in-place)
- `old_matrix` : Input matrix to be rotated
- `cell` : Primitive cell matrix (lattice vectors as columns)
- `reciprocal_vectors` : Reciprocal lattice vectors (column-wise)
- `symmetry_group` : The symmetry group (`Symmetries` or `SymmetriesQSpace`)
- `sym_index` : Index of the symmetry operation to apply (``1 \le \text{sym\_index} \le N_\text{sym}``)
- `buffer` : Optional Bumper.jl buffer for stack allocations

## Example

```julia
n_sym = get_nsymmetries(sym_group)
avg = zeros(n_modes, n_modes)
rotated = zeros(n_modes, n_modes)
for i in 1:n_sym
    rotated .= 0
    rotate_dynamical_matrix!(rotated, fc, cell, reciprocal_vectors, sym_group, i)
    avg .+= rotated
end
avg ./= n_sym  # same as symmetrize_fc!(fc, cell, sym_group)
```

## See also

- [`symmetrize_fc!`](@ref) — symmetrize a force constant matrix in real space (Cartesian)
- [`symmetrize_matrix_cartesian_q!`](@ref) — symmetrize a dynamical matrix in q space (Cartesian)
- [`rotate_matrix!`](@ref) — for intensive ``n_\text{dims} \times n_\text{dims}`` matrices (stress, dielectric)
"""
function rotate_dynamical_matrix!(new_matrix :: AbstractMatrix{T}, old_matrix :: AbstractMatrix{T}, cell :: Matrix{T},
        reciprocal_vectors :: Matrix{T}, symmetry_group :: Symmetries, sym_index :: Int; buffer=default_buffer()) where T

    n_dims = get_dimensions(symmetry_group)
    n_modes = size(old_matrix, 1)
    n_atoms = n_modes ÷ n_dims
    new_matrix .= zero(T)

    @no_escape buffer begin
        tmp_matrix = @alloc(T, n_modes, n_modes)

        # Convert Cartesian -> crystal
        cart_cryst_matrix_conversion!(tmp_matrix, old_matrix, cell; cart_to_cryst=true, buffer=buffer)

        # Apply symmetry (handles atom permutation and block-wise rotation)
        apply_sym_fc!(new_matrix, tmp_matrix, symmetry_group.symmetries[sym_index], n_dims, symmetry_group.irt[sym_index]; buffer=buffer)

        # Convert crystal -> Cartesian
        tmp_matrix .= new_matrix
        cart_cryst_matrix_conversion!(new_matrix, tmp_matrix, cell; cart_to_cryst=false, buffer=buffer)

        nothing
    end
end
function rotate_dynamical_matrix!(new_matrix :: AbstractArray{Complex{T}, 3}, old_matrix :: AbstractArray{Complex{T}, 3}, cell :: Matrix{T},
        reciprocal_vectors :: Matrix{T}, symmetry_group :: SymmetriesQSpace, sym_index :: Int; buffer=default_buffer()) where T

    n_dims = get_dimensions(symmetry_group)
    n_modes = size(old_matrix, 1)
    n_q = size(old_matrix, 3)
    new_matrix .= zero(Complex{T})

    @no_escape buffer begin
        tmp_matrix = @alloc(Complex{T}, n_modes, n_modes, n_q)

        # Convert Cartesian -> crystal (3D version loops over q-slices)
        cart_cryst_matrix_conversion!(tmp_matrix, old_matrix, cell; cart_to_cryst=true, buffer=buffer)

        # Apply symmetry (handles atom + q-point permutation + phase factors)
        apply_symmetry_matrixq!(new_matrix, tmp_matrix,
                                symmetry_group[sym_index],
                                symmetry_group.symmetries.irt[sym_index],
                                symmetry_group.irt_q[sym_index],
                                symmetry_group.symmetries.unit_cell_translations[sym_index],
                                symmetry_group.q_points;
                                buffer=buffer)

        # Convert crystal -> Cartesian
        tmp_matrix .= new_matrix
        cart_cryst_matrix_conversion!(new_matrix, tmp_matrix, cell; cart_to_cryst=false, buffer=buffer)

        nothing
    end
end


@doc raw"""
    rotate_centroid!(new_centroid :: AbstractVector{T}, old_centroid :: AbstractVector{T}, cell :: Matrix{T},
            reciprocal_vectors :: Matrix{T},
            symmetry_group :: Symmetries, sym_index :: Int; buffer=default_buffer()) where T

Apply a **single** symmetry operation on a centroid vector (e.g. atomic positions).
The centroid must be provided in Cartesian coordinates; the conversion to crystal
coordinates and back is handled internally.

This function works only in real-space, as in q space a centroid is, by definition, only a ``\Gamma`` vector.

**Real-space version** — the centroid has length ``n_\text{dims} \times n_\text{atoms}``
and the operation applies the symmetry rotation plus the atom permutation,
and includes translations if present in the symmetry group:

```math
\vec r'_{\text{irt}[a]} = S\, \vec r_{a} + \vec t
```

where ``S`` is the symmetry rotation matrix (in crystal coordinates) and ``\vec t``
is the fractional translation associated with the symmetry operation.

To symmetrize a centroid (average over all symmetries), call this function for
each symmetry and average the result. That is equivalent to what
`symmetrize_positions!` (real space) does internally.

## Parameters

- `new_centroid` : Output rotated centroid (modified in-place)
- `old_centroid` : Input centroid to be rotated
- `cell` : Primitive cell matrix (lattice vectors as columns)
- `reciprocal_vectors` : Reciprocal lattice vectors (column-wise)
- `symmetry_group` : The symmetry group (`Symmetries` or `SymmetriesQSpace`)
- `sym_index` : Index of the symmetry operation to apply (``1 \le \text{sym\_index} \le N_\text{sym}``)
- `buffer` : Optional Bumper.jl buffer for stack allocations

## Example

```julia
n_sym = get_nsymmetries(sym_group)
avg = zeros(length(centroid))
rotated = zeros(length(centroid))
for i in 1:n_sym
    rotated .= 0
    rotate_centroid!(rotated, centroid, cell, reciprocal_vectors, sym_group, i)
    avg .+= rotated
end
avg ./= n_sym  # same as symmetrize_positions!(centroid, cell, sym_group)
```

## See also

- [`symmetrize_positions!`](@ref) — symmetrize atomic positions in real space (Cartesian)
- [`rotate_vector!`](@ref) — for translation-invariant vectors (forces, displacements)
"""
function rotate_centroid!(new_centroid :: AbstractVector{T}, old_centroid :: AbstractVector{T}, cell :: Matrix{T},
        reciprocal_vectors :: Matrix{T}, symmetry_group :: Symmetries, sym_index :: Int; buffer=default_buffer()) where T

    # Get the dimensions
    n_dims = get_dimensions(symmetry_group)
    n_modes = length(new_centroid)
    n_atoms = n_modes ÷ n_dims
    new_centroid .= zero(T)

    # Get translation if available
    translation = nothing
    if sym_index <= length(symmetry_group.translations)
        translation = symmetry_group.translations[sym_index]
    end

    @no_escape buffer begin
        tmp_centroid = @alloc(T, n_modes)
        transformed_centroid = @alloc(T, n_modes)
        transformed_centroid .= zero(T)  # Initialize to zero

        # Convert to crystal coordinates
        get_crystal_coords!(reshape(tmp_centroid, n_dims, :),
                           reshape(old_centroid, n_dims, :),
                           cell; buffer=buffer)

        # Apply the symmetry with translation
        apply_sym_centroid!(transformed_centroid, tmp_centroid, 
                           symmetry_group.symmetries[sym_index], 
                           n_dims, 
                           symmetry_group.irt[sym_index];
                           translation=translation,
                           buffer=buffer)
        
        # Apply unit cell translations to bring transformed positions into the primitive cell
        for i in 1:n_atoms
            start_index = n_dims * (i - 1) + 1
            end_index = n_dims * i  
            @views transformed_centroid[start_index : end_index] .-= symmetry_group.unit_cell_translations[sym_index][:, i]
        end
       
        # Convert back to cartesian coordinates
        get_cartesian_coords!(reshape(new_centroid, n_dims, :),
                             reshape(transformed_centroid, n_dims, :),
                             cell)

        nothing
    end
end
