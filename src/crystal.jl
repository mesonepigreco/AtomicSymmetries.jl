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
function cart_cryst_matrix_conversion!(dest :: AbstractMatrix{U}, matrix :: AbstractMatrix{U}, cell :: AbstractMatrix{T}; cart_to_cryst = true, buffer=default_buffer()) where {U <: Union{T, Complex{T}}, T}
    dim = size(cell, 1)
    nmodes = size(matrix, 1)
    n_atoms = nmodes ÷ dim

    dest .= 0.0 

    @no_escape buffer begin
        metric_tensor = @alloc(T, dim, dim)
        inv_metric_tensor = @alloc(T, dim, dim)
        transform_matrix = @alloc(T, dim, dim)
        tmp_matrix = @alloc(T, dim, dim)

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

