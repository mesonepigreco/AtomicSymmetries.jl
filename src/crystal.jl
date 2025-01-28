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
