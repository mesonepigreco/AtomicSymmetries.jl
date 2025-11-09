@doc raw"""
    get_crystal_coords!(crystal :: AbstractMatrix{T}, cartesian :: AbstractMatrix{T}, cell :: AbstractMatrix{T}; buffer = default_buffer())

Convert cartesian coordinates into crystal ones. 
Optionally, a Bumper buffer can be provided to avoid memory allocations.

"""
function get_crystal_coords!(crystal :: AbstractMatrix{T}, cartesian :: AbstractMatrix{T}, cell :: AbstractMatrix{T}; buffer = default_buffer()) where T
    n_atoms = size(cartesian, 2)
    @no_escape buffer begin
        metric_tensor = @alloc(T, 3, 3)
        inv_metric_tensor = @alloc(T, 3, 3)
        tmp_vectors = @alloc(T, 3, n_atoms)

        mul!(metric_tensor, cell', cell)
        inv_metric_tensor .= inv(metric_tensor)

        mul!(tmp_vectors, cell', cartesian)
        mul!(crystal, inv_metric_tensor, tmp_vectors)
        nothing # <-- Avoid returning crystal
    end
end

function get_cartesian_coords!(cartesian :: AbstractMatrix{T}, crystal :: AbstractMatrix{T}, cell :: AbstractMatrix{T}) where T
    mul!(cartesian, cell, crystal)
end

@doc raw"""
    to_primitive_cell_cryst!(cryst_coords :: AbstractVector{T}, closest_vector :: AbstractVector{T})

Convert the crystal coordinates into the primitive cell closest
to the provided vector
"""
function to_primitive_cell_cryst!(cryst_coords :: AbstractVector{T}, closest_vector :: AbstractVector{T}) where {T}
    ndim = length(cryst_coords)
    for i in 1:ndim
        cryst_coords[i] = cryst_coords[i] - round(cryst_coords[i] - closest_vector[i])
    end
end



@doc raw"""
    to_primitive_cell_cart!(cartesian_coords :: AbstractVector{T}, cell :: AbstractMatrix{T}; buffer = default_buffer())
    to_primitive_cell_cart!(cartesian_coords :: AbstractMatrix{T}, cell :: AbstractMatrix{T}; buffer = default_buffer())

Put the atoms as closest as possible to the origin of the cell.
"""
function to_primitive_cell_cart!(cartesian_coords :: AbstractVector{T}, cell :: AbstractMatrix{T}; buffer = default_buffer()) where T
    @no_escape buffer begin
        crystal_coords = @alloc(T, size(cartesian_coords)...)
        get_crystal_coords!(crystal_coords, cartesian_coords, cell)
        to_primitive_cell_cryst!(crystal_coords, [0.0, 0.0, 0.0])
        get_cartesian_coords!(cartesian_coords, crystal_coords, cell)
        nothing
    end
end
function to_primitive_cell_cart!(cartesian_coords :: AbstractMatrix{T}, cell :: AbstractMatrix{T}; buffer = default_buffer()) where T
    nat = size(cartesian_coords, 2)
    @no_escape buffer begin
        crystal_coords = @alloc(T, size(cartesian_coords)...)
        get_crystal_coords!(crystal_coords, cartesian_coords, cell; buffer=buffer)
        zvect = @alloc(T, 3)
        zvect .= 0.0 

        for i in 1:nat
            @views to_primitive_cell_cryst!(crystal_coords[:, i], zvect)
        end
        get_cartesian_coords!(cartesian_coords, crystal_coords, cell)
        nothing
    end
end




@doc """
    cart_cryst_matrix_conversion!(dest :: AbstractMatrix{T}, matrix :: AbstractMatrix{T}, cell :: AbstractMatrix{T}; cart_to_cryst :: Bool = true, buffer=default_buffer()) where T
    cart_cryst_matrix_conversion!(dest :: AbstractArray{Complex{T}, 3}, matrix :: AbstractArray{Complex{T}, 3}, cell :: AbstractMatrix{T}; cart_to_cryst = true, buffer=default_buffer()) where T


Convert a matrix `matrix` from cartesian to crystal coordinates.
The result are stored in `dest`.
The primitive vectors are stored as columns of the `cell` matrix.

If `cart_to_cryst` is `true`, the conversion is from cartesian to crystal coordinates.
Otherwise, the conversion is from crystal to cartesian coordinates.

This function exploits Bumper stack memory allocations. 
It is possible to pass the stack as an argument with the buffer keyword

It works either with supercell real force constant matrices, and with force constant matrices directly written in q space.
In the latter case, the dimension of the matrix is expected to be (n_modes, n_modes, nq)
"""
function cart_cryst_matrix_conversion!(dest :: AbstractMatrix{U}, matrix :: AbstractMatrix{U}, cell :: AbstractMatrix{T}; cart_to_cryst = true, buffer=default_buffer()) where {T, U <: Union{T, Complex{T}}}
    dim = size(cell, 1)
    nmodes = size(matrix, 1)
    n_atoms = nmodes ÷ dim

    dest .= 0.0 

    @no_escape buffer begin
        metric_tensor = @alloc(T, dim, dim)
        inv_metric_tensor = @alloc(T, dim, dim)
        transform_matrix = @alloc(T, dim, dim)
        tmp_matrix = @alloc(U, dim, dim)

        mul!(metric_tensor, cell', cell)
        inv_metric_tensor .= inv(metric_tensor) #TODO: Allocating
        mul!(transform_matrix, inv_metric_tensor, cell')

        if cart_to_cryst
            transform_matrix .= inv(transform_matrix)
        end

        for i in 1:n_atoms
            for j in 1:n_atoms
                @views mul!(tmp_matrix, matrix[dim*(i-1) + 1 : dim*i, dim*(j-1)+1 : dim*j], 
                            transform_matrix)
                @views mul!(dest[dim*(i-1) + 1 : dim*i, dim*(j-1)+1 : dim*j], transform_matrix', tmp_matrix)
            end
        end
        nothing
    end
end
function cart_cryst_matrix_conversion!(dest :: AbstractArray{Complex{T}, 3}, 
        matrix :: AbstractArray{Complex{T}, 3}, cell :: AbstractMatrix{T}; cart_to_cryst = true, buffer=default_buffer()) where T
    nq = size(dest, 3)
    @assert nq == size(matrix, 3)

    for iq in 1:nq
        @views cart_cryst_matrix_conversion!(dest[:, :, iq], matrix[:, :, iq], cell; 
                                             cart_to_cryst = cart_to_cryst, buffer=buffer)
    end
end


@doc raw"""
    get_reciprocal_lattice!(reciprocal_vectors :: Matrix{T}, cell :: Matrix{T})

Compute the reciprocal lattice vectors from the primitive cell.

Reciprocal lattice vectors ``\boldsymbol{B}`` (columns of the matrix)
satisfy the property

```math
\boldsymbol{B}^\dagger \boldsymbol{A} = 2\pi \boldsymbol{I}
```

where ``\boldsymbol{A}`` is the matrix whose columns are the direct lattice
vectors, and ``\boldsymbol{I}`` is the identity matrix.


## Parameters

- `reciprocal_vectors` : The ``\boldsymbol{B}`` matrix, whose columns are
    the reciprocal lattice vectors (modified in-place).
- `cell` : The matrix whose columns are the primitive direct lattice vectors.
"""
function get_reciprocal_lattice!(reciprocal_vectors :: Matrix{T}, cell :: Matrix{T}) where T
    n_dims = size(cell, 1)
    tmp_cell_t = SMatrix{n_dims, n_dims}(cell')
    reciprocal_vectors .= inv(tmp_cell_t)
    reciprocal_vectors .*= (2π)
end





@doc raw"""
    cryst_cart_conv!(target, source, primitive_cell, reciprocal_vectors, cryst_to_cart; q_space=false)

In-place conversion of coordinates between crystallographic and Cartesian systems,
supporting both real and reciprocal (q) space.

The function performs one of four transformations based on the boolean flags
`cryst_to_cart` and `q_space`. It computes `target = α * T * source`, where
`T` is the transformation matrix and `α` is a scaling factor.

# Arguments
- `target::AbstractArray{T}`: The destination array, which is modified in-place.
- `source::AbstractArray{T}`: The source array containing the coordinates to be transformed.
- `primitive_cell::AbstractMatrix{U}`: The matrix whose columns represent the
  primitive lattice vectors (e.g., $A = [a₁, a₂, a₃]$).
- `reciprocal_vectors::AbstractMatrix{U}`: The matrix whose columns represent the
  reciprocal lattice vectors (e.g., $B = [b₁, b₂, b₃]$).
- `cryst_to_cart::Bool`: The direction of the transformation.
  - `true`: Crystallographic coordinates (unitless) to Cartesian (units of length or 1/length).
  - `false`: Cartesian to Crystallographic.

# Keyword Arguments
- `q_space::Bool = false`: Toggles between real space and reciprocal (q) space.
  - `false`: Real-space transformation.
  - `true`: Reciprocal-space (q-space) transformation.

# Operations Performed

Let `A = primitive_cell` and `B = reciprocal_vectors`. The function assumes the
standard physics definition where ``A^T B = 2\pi I``.



The function calculates `target = α * T * source` based on the following cases:

1.  **`cryst_to_cart=true`, `q_space=false`**: (Cryst → Cart, Real Space)
    - `T = A`
    - `α = 1.0`
    - `target = A * source`

2.  **`cryst_to_cart=true`, `q_space=true`**: (Cryst → Cart, Q-Space)
    - `T = B`
    - `α = 1.0`
    - `target = B * source`

3.  **`cryst_to_cart=false`, `q_space=false`**: (Cart → Cryst, Real Space)
    - `T = B'` (Transpose of `reciprocal_vectors`)
    - `α = 1.0`
    - `target = B' * source`
4.  **`cryst_to_cart=false`, `q_space=true`**: (Cart → Cryst, Q-Space)
    - `T = A'` (Transpose of `primitive_cell`)
    - `α = 1 / (2π)`
    - `target = (1 / (2π)) * A' * source` (This is correct, as ``B^{-1} = \frac{1}{2\pi} A^T``)
"""
function cryst_cart_conv!(target :: AbstractArray{T}, source :: AbstractArray{T},
        primitive_cell :: AbstractMatrix{U},
        reciprocal_vectors :: AbstractMatrix{U},
        cryst_to_cart :: Bool;
        q_space :: Bool = false) where {T, U}

    if cryst_to_cart && q_space
        mul!(target, reciprocal_vectors, source)
    elseif !cryst_to_cart && !q_space
        mul!(target, reciprocal_vectors', source, U(1/(2π)), zero(U))
    elseif !cryst_to_cart && q_space
        mul!(target, primitive_cell', source, U(1/(2π)), zero(U))
    else
        mul!(target, primitive_cell, source)
    end
end
